library(splines)

#' Fit Lifeplus model
#'
#' @param data a data frame.
#' @param y column name of outcome variable..
#' @param time column name of time variable.
#' @param area column name of the area of each observation
#' @param start_time start time of estimates.
#' @param end_time end time of estimates.
#' @param transition_model transition model specification (default: `transition_model_logistic()`)
#' @param data_model data model specification (default: `data_model_normal()`)
#' @param shock_model shock term specification (default: `shock_model_none()`)
#' @param held_out binary vector indicating which observations are held out. Set to FALSE to hold out no observations.
#' @param ... additional arguments passed to CmdStanModel::sample.
#'
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
  held_out = FALSE,
  ...
) {
  ###### Initial argument checks
  checkmate::check_class(transition_model, "lifeplus_transition_model")
  checkmate::check_class(data_model, "lifeplus_data_model")
  checkmate::check_class(shock_model, "lifeplus_shock_model")
  checkmate::check_data_frame(data, min.rows = 1)
  checkmate::check_subset(y, colnames(data))
  checkmate::check_subset(time, colnames(data))
  checkmate::check_subset(area, colnames(data))
  checkmate::check_integer(start_time, upper = end_time, null.ok = TRUE)
  checkmate::check_integer(end_time, lower = start_time, null.ok = TRUE)
  checkmate::check_logical(held_out, len = nrow(data))

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

  # Make sure the observed data are within the estimation period
  checkmate::check_integer(data[[time]], lower = start_time, upper = end_time)

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
      stan_data <- c(stan_data, data_model$stan_data(data[[y]], grid))
    } else {
      stan_data <- c(stan_data, data_model$stan_data)
    }
  }

  # Augment stan_data with additional elements from transition model
  if (!is.null(transition_model$stan_data)) {
    if (is.function(transition_model$stan_data)) {
      stan_data <- c(stan_data, transition_model$stan_data(data[[y]], grid))
    } else {
      stan_data <- c(stan_data, transition_model$stan_data)
    }
  }

  # Augment stan_data with additional elements from shock model
  if (!is.null(shock_model$stan_data)) {
    if (is.function(shock_model$stan_data)) {
      stan_data <- c(stan_data, shock_model$stan_data(data[[y]], grid))
    } else {
      stan_data <- c(stan_data, shock_model$stan_data)
    }
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
