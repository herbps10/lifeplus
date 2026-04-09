#' Parse Stan variable indices from CmdStan summary output
#' @param variables Character vector of Stan-style index names (e.g. \code{"eta[1,2]"})
#'
#' @importFrom stringr str_match str_split
#' @return A list with \code{name} (character) and \code{indices} (integer matrix)
#' @keywords internal
#' @noRd
parse_stan_indices <- function(variables) {
  base_match <- stringr::str_match(
    variables,
    "^([a-zA-Z_][a-zA-Z0-9_]*)\\[(.+)\\]$"
  )

  index_strings <- stringr::str_split(base_match[, 3], ",")
  indices <- do.call(rbind, lapply(index_strings, as.integer))
  list(name = base_match[, 2], indices = indices)
}

#' Summarise draws from CmdStan fit
#' @keywords internal
#' @noRd
summarise_draws <- function(fit, variables, n_cores, posterior_quantiles) {
  fit$samples$summary(
    variables,
    ~ stats::quantile(.x, probs = posterior_quantiles),
    posterior::default_convergence_measures(),
    .cores = n_cores
  )
}

#' Check if transition model was estimated hierarchically
#' @keywords internal
#' @noRd
is_hierarchical_transition <- function(fit) {
  h <- fit$transition_model$hierarchical
  !is.null(h) && isTRUE(as.logical(h))
}

#' Default summary posterior quantile probabilities
#' @details
#' The default quantiles are 0.1%, 1%, 2.5%, 25%, 50%, 75%, 90%, 97.5%, 99%, 99.9%.
#'
#' @return numeric vector of quantile probabilities
#' @export
lifeplus_default_posterior_quantiles <- function() {
  c(0.001, 0.01, 0.025, 0.1, 0.25, 0.5, 0.75, 0.9, 0.975, 0.99, 0.999)
}

#' Extract posterior summaries for temporal (area x time) parameters
#'
#' @param fit A \code{lifeplus} object.
#' @param n_cores Integer number of cores.
#' @param posterior_quantiles Numeric vector of quantile probabilities.
#'
#' @return a data frame with posterior quantiles, convergence diagnostics,
#'   and columns \code{area} and \code{time}
#'
#' @noRd
#' @keywords internal
process_temporal <- function(fit, var, n_cores, posterior_quantiles) {
  post <- summarise_draws(fit, var, n_cores, posterior_quantiles)

  parsed <- parse_stan_indices(post$variable)
  post$variable <- parsed$name
  post$area <- fit$areas[parsed$indices[, 1]]
  post$time <- fit$times[parsed$indices[, 2]]

  return(post)
}

#' Extract posterior summaries for area-level transition functions
#'
#' Processes \code{transition_function_pred[c, g]} over the age grid
#'
#' @param fit A \code{lifeplus} object.
#' @param n_cores Integer number of cores.
#' @param posterior_quantiles Numeric vector of quantile probabilities.
#'
#' @return a data frame with posterior quantiles, convergence diagnostics,
#'   and columns \code{area} and \code{x} (age grid value)
#'
#' @noRd
#' @keywords internal
process_transition_functions <- function(fit, n_cores, posterior_quantiles) {
  post <- summarise_draws(
    fit,
    "transition_function_pred",
    n_cores,
    posterior_quantiles
  )

  parsed <- parse_stan_indices(post$variable)
  post$variable <- parsed$name
  post$area <- fit$areas[parsed$indices[, 1]]
  post$x <- fit$grid[parsed$indices[, 2]]

  return(post)
}

#' Extract posterior summaries for the hierarchical mean transition function
#'
#' Only applicable when the transition model uses a hierarchical specification.
#'
#' @param fit A \code{lifeplus} object.
#' @param n_cores Integer number of cores.
#' @param posterior_quantiles Numeric vector of quantile probabilities.
#'
#' @return a data frame with posterior quantiles, convergence diagnostics,
#'   and column \code{x} (age grid value), or \code{NULL}
#'
#' @noRd
#' @keywords internal
process_transition_function_mean <- function(
  fit,
  n_cores,
  posterior_quantiles
) {
  if (!is_hierarchical_transition(fit)) {
    return(NULL)
  }

  post <- summarise_draws(
    fit,
    "transition_function_pred_mean",
    n_cores,
    posterior_quantiles
  )

  parsed <- parse_stan_indices(post$variable)
  post$variable <- parsed$name
  post$x <- fit$grid[parsed$indices[, 1]]

  return(post)
}

#' Extract posterior summaries for transition model parameters
#' Delegates to the transition model's \code{extract_params} method if one
#' is defined. Returns \code{NULL} for models that do not expose parameters.
#'
#' @param fit A \code{lifeplus} object.
#' @param n_cores Integer number of cores.
#' @param posterior_quantiles Numeric vector of quantile probabilities.
#'
#' @return A list with element \code{params} (data frame) and optionally
#' \code{correlation}, or \code{NULL}
#'
#' @noRd
#' @keywords internal
process_transition_params <- function(fit, n_cores, posterior_quantiles) {
  extractor <- fit$transition_model$extract_params
  if (is.null(extractor)) {
    return(NULL)
  }
  extractor(fit, n_cores, posterior_quantiles)
}


#' Extract posterior summaries for shock model parameters
#' Delegates to the shock model's \code{extract_params} method if one
#' is defined. Returns \code{NULL} for models that do not expose parameters.
#'
#' @param fit A \code{lifeplus} object.
#' @param n_cores Integer number of cores.
#' @param posterior_quantiles Numeric vector of quantile probabilities.
#'
#' @return A list with element \code{params} (data frame) or \code{NULL}
#'
#' @noRd
#' @keywords internal
process_shock_params <- function(fit, n_cores, posterior_quantiles) {
  extractor <- fit$shock_model$extract_params
  if (is.null(extractor)) {
    return(NULL)
  }
  extractor(fit, n_cores, posterior_quantiles)
}


#' Extract posterior summaries for data model parameters
#' Delegates to the data model's \code{extract_params} method if one
#' is defined. Returns \code{NULL} for models that do not expose parameters.
#'
#' @param fit A \code{lifeplus} object.
#' @param n_cores Integer number of cores.
#' @param posterior_quantiles Numeric vector of quantile probabilities.
#'
#' @return A list with element \code{params} (data frame) or \code{NULL}
#'
#' @noRd
#' @keywords internal
process_shock_params <- function(fit, n_cores, posterior_quantiles) {
  extractor <- fit$data_model$extract_params
  if (is.null(extractor)) {
    return(NULL)
  }
  extractor(fit, n_cores, posterior_quantiles)
}

#' Extract and summarise all posteriors from a lifeplus fit
#'
#' Called internally by \code{\link{lifeplus}()} after sampling. Computes
#' quantile summaries and convergence diagnostics for all major parameter groups.
#'
#' @param fit A \code{lifeplus} object (must constain \code{samples}, \code{areas},
#'   \code{times}, \code{grid}, \code{shock_model}, and \code{transition_model}).
#' @param n_cores Integer number of cores for parallel summary computation (default: 1).
#' @param posterior_quantiles Summary quantiles (defaults to \code{\link{lifeplus_default_posterior_quantiles}()}).
#'
#' @return A named list with elements:
#' \describe{
#'   \item{\code{life}}{Area x time posterior life expectancy summaries.}
#'   \item{\code{shocks}}{Area x time posterior shock summaries.}
#'   \item{\code{shockfree}}{Area x time posterior shock-free life expectancy summaries.}
#'   \item{\code{transition_functions}}{Area x age-grid transition function summaries.}
#'   \item{\code{transition_function_mean}}{Hierarchical mean transition function summaries (or \code{NULL})}
#'   \item{\code{transition_params}}{Transition model parameter summaries (or \code{NULL})}
#' }
#'
#' @noRd
#' @keywords internal
lifeplus_posteriors <- function(
  fit,
  n_cores = 1L,
  posterior_quantiles = lifeplus_default_posterior_quantiles()
) {
  value <- list(
    life = process_temporal(fit, "eta", n_cores, posterior_quantiles),
    transition_functions = process_transition_functions(
      fit,
      n_cores,
      posterior_quantiles
    ),
    transition_function_mean = process_transition_function_mean(
      fit,
      n_cores,
      posterior_quantiles
    ),
    transition_params = process_transition_params(
      fit,
      n_cores,
      posterior_quantiles
    ),
  )

  if (fit$shock_model$name != "none") {
    value$shockfree <- process_temporal(
      fit,
      "eta_shockfree",
      n_cores,
      posterior_quantiles
    )

    value$shocks <- process_temporal(
      fit,
      "shock2",
      n_cores,
      posterior_quantiles
    )

    value_shock_params <- process_shock_params(
      fit,
      n_cores,
      posterior_quantiles
    )
  }
}
