#' @importFrom generics tidy
#' @export
generics::tidy

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
#' @param posterior_quantiles A numeric vector of probabilities in \code{(0, 1)}
#'   at which to summarise posterior draws. Defaults to \code{\link{lifeplus_default_quantiles}()}, which provides
#'   a fine grid from the 0.1% to 99.9% quantile.
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
#'   \item{\code{posterior_quantiles}}{The numeric vector of quantile probabilities used for posterior summaries.}
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
  posterior_quantiles = lifeplus_default_posterior_quantiles(),
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
  checkmate::assert_numeric(start_time, null.ok = TRUE)
  checkmate::assert_numeric(end_time, null.ok = TRUE)
  checkmate::assert_logical(held_out, len = nrow(data), null.ok = TRUE)

  checkmate::assert_numeric(
    posterior_quantiles,
    lower = 0,
    upper = 1,
    any.missing = FALSE,
    min.len = 1,
    sorted = TRUE,
    .var.name = "posterior_quantiles"
  )

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
  start_time <- as.integer(start_time)
  end_time <- as.integer(end_time)

  if (is.null(held_out)) {
    held_out <- rep(FALSE, nrow(data))
  }

  # Make sure the observed data are within the estimation period
  checkmate::assert_numeric(data[[time]], lower = start_time, upper = end_time)

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
    ifelse(is.null(args$parallel_chains), 1, args$parallel_chains),
    posterior_quantiles = posterior_quantiles
  )

  result$diagnose <- fit$diagnostic_summary()
  result$posterior_quantiles <- posterior_quantiles

  attr(result, "class") <- "lifeplus"

  result
}
