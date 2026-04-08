library(splines)

#' Fit a Lifeplus Model
#'
#' @description
#' Fits a Bayesian hierarchical model to life expectancy data using
#' CmdStan. The model is composed of three components: a transition model,
#' data model, and shock model, each of which may be configured independently.
#'
#' @details
#' The function constructs a model from the chosen specifications, prepares
#' the data, runs MCMC sampling via \code{\link[cmdstanr]{CmdStanModel}}, and returns
#' an object of class \code{lifeplus} containing posterior draws, diagnostics, and
#' metadata.
#'
#' All areas are assumed to have observations at every time point in the (non-held-out)
#' data.
#'
#' @param data a data frame containing at least the columns specified by \code{y}, \code{time}, and \code{area}.
#' @param y A string giving the name of the numeric outcome column in \code{data} (typically life expectancy in years).
#' @param time A string giving the name of the integer-valued time column in \code{data} (e.g., calendar year).
#' @param area A string giving the name of the column in \code{data} that identifies the geographic or administrative unit for each observation.
#' @param start_time An integer specifying the first time point of the estimation window, or \code{NULL} (default)
#'   to use \code{min(data[[time]])}.
#' @param end_time An integer specifying the last time point of the projection window, or \code{NULL} (default)
#'   to use \code{max(data[[time]])}.
#' @param transition_model A transition model specification object created by a
#'   \code{transition_model_*} constructor (default: \code{\link{transition_model_double_logistic}()}).
#' @param data_model A data model specification object created by a
#'   \code{data_model_*} constructor (default: \code{\link{data_model_normal}()}).
#' @param shock_model A shock model specification object created by a
#'   \code{shock_model_*} constructor (default: \code{\link{shock_model_none}()}).
#' @param held_out A logical vector of length \code{nrow(data)} indicating
#'   which observations to hold out from model fitting, or \code{NULL} (default)
#'   to include all observations.
#' @param ... additional arguments passed to \code{CmdStanModel::sample}.
#'
#' @return An object of class \code{lifeplus}, which is a named list containing:
#' \describe{
#'   \item{\code{samples}}{A \code{\link[cmdstanr]{CmdStanMCMC}} fit object.}
#'   \item{\code{data}}{The original input data frame.}
#'   \item{\code{stan_data}}{The named list passed to Stan.}
#'   \item{\code{times}}{Integer vector of time-points in the estimation and projection window (\code{start_time:end_time})}
#'   \item{\code{areas}}{Vector of unique area identifiers.}
#'   \item{\code{grid}}{Grid of transition function inputs.}
#'   \item{\code{elapsed}}{A \code{link{difftime}} object for wall-clock sampling time.}
#'   \item{\code{posteriors}}{Posterior summaries.}
#'   \item{\code{diagnose}}{MCMC diagnostic sumary from CmdStan.}
#'   \item{\code{transition_model}, \code{data_model}, \code{shock_model}}{Model specifications used for this fit.}
#' }
#'
#' @examples
#' \dontrun{
#' example_life_data <- simulate_lifeplus(num_areas = 1, num_times = 40)
#'
#' fit <- lifeplus(
#'   data = example_life_data,
#'   y = "y",
#'   time = "time",
#'   area = "area",
#'   transition_model = transition_model_double_logistic(hierarchical = FALSE),
#'   data_model = data_model_normal(),
#'   chains = 4,
#'   parallel_chains = 4,
#'   iter_warmup = 1000,
#'   iter_sampling = 1000
#' )
#' }
#'
#' @seealso
#' Transition models: \code{\link{transition_model_double_logistic}()},
#'    \code{\link{transition_model_spline}()}, \code{\link{transition_model_gaussian_process}()},
#'
#' Data models: \code{\link{data_model_normal}()}, \code{\link{data_model_outlier}()}
#'
#' Shock models: \code{\link{shock_model_none}()}, \code{\link{shock_model_regularized_horseshoe}()}
#'
#' @export
lifeplus <- function(
  data,
  y,
  time,
  area,
  start_time = NULL,
  end_time = NULL,
  transition_model = transition_model_double_logistic(),
  data_model = data_model_normal(),
  shock_model = shock_model_none(),
  held_out = NULL,
  ...
) {
  ###### Initial argument checks
  checkmate::assert_class(transition_model, "lifeplus_transition_model")
  checkmate::assert_class(data_model, "lifeplus_data_model")
  checkmate::assert_class(shock_model, "lifeplus_shock_model")
  checkmate::assert_data_frame(data, min.rows = 1)
  checkmate::assert_subset(y, colnames(data))
  checkmate::assert_subset(time, colnames(data))
  checkmate::assert_subset(area, colnames(data))
  checkmate::assert_integer(start_time, null.ok = TRUE)
  checkmate::assert_integer(end_time, null.ok = TRUE)
  checkmate::assert_logical(held_out, len = nrow(data), null.ok = TRUE)

  args <- list(...)

  # Save original dataset
  original_data <- data

  # Initialize start and end time if necessary
  if (is.null(start_time)) {
    start_time <- min(data[[time]])
  }
  if (is.null(end_time)) {
    end_time <- max(data[[time]])
  }

  if (is.null(held_out)) {
    held_out <- rep(FALSE, nrow(data))
  }

  # Make sure the observed data are within the estimation period
  checkmate::assert_integer(data[[time]], lower = start_time, upper = end_time)

  check_model_compatibility(transition_model, data_model, shock_model)

  # Load model
  stan_model <- load_model(transition_model, data_model, shock_model)
  # Setup data for Stan
  areas <- unique(data[[area]])
  times <- seq(start_time, end_time, 1)

  data$c <- match(data[[area]], areas)
  data$t <- match(data[[time]], times)

  if (length(held_out) == 1 && held_out == FALSE) {
    held_out = rep(0, nrow(data))
  } else {
    held_out = as.numeric(held_out)
  }

  t_last <- max(data$t)

  # Set up grid
  grid <- c(seq(from = 0, to = 110, by = 1))
  num_grid <- length(grid)

  obs <- matrix(
    NA,
    nrow = length(areas),
    ncol = max(data$t[held_out == 0])
  )
  for (c in seq_along(areas)) {
    o <- order(data$t[data$c == c & held_out == 0])
    obs[c, ] <- data[[y]][data$c == c & held_out == 0][o]
  }

  stan_data <- list(
    C = nrow(obs),
    T = ncol(obs),
    Tpred = length(times),

    y = obs,

    num_grid = num_grid,
    grid = grid,

    shock_diff_mode = 0,
    include_prior = 0
  )

  # Augment stan_data with additional elements from data model
  if (!is.null(data_model$stan_data)) {
    if (is.function(data_model$stan_data)) {
      data_model$stan_data <- data_model$stan_data(data[[y]], grid)
    }
    stan_data <- c(stan_data, data_model$stan_data)
  }

  # Augment stan_data with additional elements from transition model
  if (!is.null(transition_model$stan_data)) {
    if (is.function(transition_model$stan_data)) {
      transition_model$stan_data <- transition_model$stan_data(data[[y]], grid)
    }
    stan_data <- c(stan_data, transition_model$stan_data)
  }

  # Augment stan_data with additional elements from shock model
  if (!is.null(shock_model$stan_data)) {
    if (is.function(shock_model$stan_data)) {
      shock_model$stan_data <- shock_model$stan_data(data[[y]], grid)
    }
    stan_data <- c(stan_data, shock_model$stan_data)
  }

  # Sampling
  start <- Sys.time()
  fit <- stan_model$sample(
    stan_data,
    ...
  )
  elapsed <- Sys.time() - start

  result <- list(
    samples = fit,

    data = original_data,
    stan_data = stan_data,

    times = times,
    areas = areas,

    grid = grid,

    elapsed = elapsed,

    # Save arguments
    y = y,
    time = time,
    area = area,
    held_out = held_out,

    transition_model = transition_model,
    data_model = data_model,
    shock_model = shock_model
  )

  cat("Extracting posteriors...\n")

  result$posteriors <- lifeplus_posteriors(
    result,
    ifelse(is.null(args$parallel_chains), 1, args$parallel_chains)
  )

  result$diagnose <- fit$diagnostic_summary()

  attr(result, "class") <- "lifeplus"

  result
}

#' Print a Lifeplus Model Fit
#'
#' @description
#' Displays a concise summary of a fitted \code{lifeplus} model.
#'
#' @param x An object of class \code{lifeplus}, as returned by \code{\link{lifeplus}()}.
#' @param ... Additional arguments passed to other methods (currently ignored).
#'
#' @return Invisibly returns \code{x}.
#'
#' @seealso \code{\link{lifeplus}()}
#'
#' @export
print.lifeplus <- function(x, ...) {
  cli::cli_rule(left = "{.strong lifeplus} model fit")

  ##### Inputs #####
  n_areas <- length(x$areas)
  n_obs <- nrow(x$data)
  n_held_out <- sum(x$held_out)
  start_time <- min(x$times)
  end_time <- max(x$times)
  n_time <- length(x$times)

  cli::cli_h3("Data")
  cli::cli_ul()
  cli::cli_li("Outcome variable: {.field {x$y}}")
  cli::cli_li(
    "Time variable: {.field {x$time}} ({.val {start_time}} to {.val {end_time}})"
  )
  cli::cli_li(
    "Area variable: {.field {x$area}} ({.val {n_areas}} areas)"
  )
  cli::cli_li("Observations: {.val {n_obs}} ({.val {n_held_out}} held out)")
  cli::cli_end()

  ##### Model components ######
  cli::cli_h3("Model components")
  cli::cli_ul()
  cli::cli_li("Transition: {.emph {x$transition_model$name}}")
  cli::cli_li("Data: {.emph {x$data_model$name}}")
  cli::cli_li("Shock: {.emph {x$shock_model$name}}")
  cli::cli_end()

  ##### Sampling #####
  cli::cli_h3("Sampling")
  n_chains <- length(x$samples$metadata()$id)
  n_iter <- x$samples$metadata()$iter_sampling
  n_warmup <- x$samples$metadata()$iter_warmup

  cli::cli_ul()
  cli::cli_li("Chains: {.val {n_chains}}")
  cli::cli_li("Iterations: {.val {n_warmup}} warmup + {.val {n_iter}} sampling")
  cli::cli_li("Elapsed: {format(x$elapsed, digits = 3)}")
  cli::cli_end()

  ##### Diagnostics #####
  diag <- x$diagnose
  n_divergent <- sum(diag$num_divergent)
  n_max_treedepth <- sum(diag$num_max_treedepth)
  any_low_ebfmi <- any(diag$ebfmi < 0.3)

  cli::cli_h3("Diagnostics")
  if (n_divergent > 0) {
    cli::cli_alert_danger("{.val {n_divergent}} divergent transition{?s}")
  } else {
    cli::cli_alert_success("No divergent transitions")
  }
  if (n_max_treedepth > 0) {
    cli::cli_alert_danger(
      "{.val {n_max_treedepth}} transition{?s} hit max treedepth"
    )
  } else {
    cli::cli_alert_success("No max treedepth warnings")
  }
  if (any_low_ebfmi) {
    cli::cli_alert_danger("Low E-BFMI detected (< 0.3)")
  } else {
    cli::cli_alert_success("No E-BFMI warnings")
  }

  invisible(x)
}
