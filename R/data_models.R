#' Normal Data Model
#'
#' @description
#' White noise data model that assumes the residuals follow independent mean-zero normal distributions.
#'
#' @details
#' The residuals are modeled as \eqn{\epsilon_{c,t} \sim N(0, \sigma^2)}.
#' The prior for the residual scale parameter is \eqn{\sigma \sim N(m, s^2)},
#' where \eqn{m} and \eqn{s} are user-specified via the \code{prior_mean} and \code{prior_sd} arguments, respectively.
#'
#' @param prior_mean Prior mean for white noise standard deviation
#' @param prior_sd Prior standard deviation for white noise standard deviation
#' @param fixed_sd if \code{NULL} (default), then residual standard deviation parameter is estimated from the data; if set to
#'   a positive number, then residual standard deviation is fixed to that value.
#'
#' @return lifeplus_data_model
#' @export
data_model_normal <- function(prior_mean = 0, prior_sd = 1, fixed_sd = NULL) {
  checkmate::assert_number(fixed_sd, null.ok = TRUE, lower = 0)
  structure(
    list(
      name = "normal",
      prior_mean = prior_mean,
      prior_sd = prior_sd,
      fixed_sd = fixed_sd,
      stan_data = list(
        fix_epsilon_sigma = !is.null(fixed_sd),
        epsilon_sigma_fixed = ifelse(is.null(fixed_sd), 0, fixed_sd),
        epsilon_sigma_prior_mu = prior_mean,
        epsilon_sigma_prior_sd = prior_sd
      ),
      print_info = normal_print_info,
      extract_params = normal_extract_params
    ),
    class = "lifeplus_data_model"
  )
}

#' Extract posterior summaries for normal data model parameters
#' @noRd
#' @keywords internal
normal_extract_params <- function(fit, n_cores, posterior_quantiles) {
  vars <- c("epsilon_sigma")
  post <- summarise_draws(fit, vars, n_cores, posterior_quantiles)

  list(params = post)
}

#' Print details for the normal data model
#'
#' @param x A \code{lifeplus_data_model} object with \code{name = "data_model"}.
#'
#' @noRd
#' @keywords internal
normal_print_info <- function(x) {
  cli::cli_h3("Settings")
  cli::cli_ul()
  if (is.null(x$fixed_sd)) {
    cli::cli_li("Prior mean: {.val {x$prior_mean}}")
    cli::cli_li("Prior standard deviation: {.val {x$prior_sd}}")
  } else {
    cli::cli("Fixed standard deviation: {.val {x$fixed_sd}}")
  }
  cli::cli_end()
}

#' Outlier Data Model
#'
#'
#' @details
#' Any observed absolute differences in life expectancy greater than the user-specified outlier threshold (`outlier_threshold`) are ignored.
#' The remaining residuals are modeled as in the normal data model:  \eqn{\epsilon_{c,t} \sim N(0, \sigma^2)}.
#' The prior for the residual scale parameter is \eqn{\sigma \sim N(m, s^2)},
#' where \eqn{m} and \eqn{s} are user-specified via the \code{prior_mean} and \code{prior_sd} arguments, respectively.
#'
#' @param outlier_threshold threshold for ignoring observed life expectancy differences.
#' @param prior_mean Prior mean for white noise standard deviation
#' @param prior_sd Prior standard deviation for white noise standard deviation
#' @param fixed_sd if \code{NULL} (default), then residual standard deviation parameter is estimated from the data; if set to
#'   a positive number, then residual standard deviation is fixed to that value.
#'
#' @return lifeplus_data_model
#' @export
data_model_outlier <- function(
  outlier_threshold = 5,
  prior_mean = 0,
  prior_sd = 1,
  fixed_sd = NULL
) {
  structure(
    list(
      name = "outlier",
      outlier_threshold = outlier_threshold,
      prior_mean = prior_mean,
      prior_sd = prior_sd,
      fixed_sd = fixed_sd,
      stan_data = list(
        fix_epsilon_sigma = !is.null(fixed_sd),
        epsilon_sigma_fixed = ifelse(is.null(fixed_sd), 0, fixed_sd),
        outlier_threshold = outlier_threshold,
        epsilon_sigma_prior_mu = prior_mean,
        epsilon_sigma_prior_sd = prior_sd
      ),
      print_info = outlier_print_info,
      extract_params = outlier_extract_params
    ),
    class = "lifeplus_data_model"
  )
}

#' Extract posterior summaries for outlier data model parameters
#' @noRd
#' @keywords internal
outlier_extract_params <- function(fit, n_cores, posterior_quantiles) {
  vars <- c("epsilon_sigma")
  post <- summarise_draws(fit, vars, n_cores, posterior_quantiles)

  list(params = post)
}

#' Print details for the outlier data model
#'
#' @param x A \code{lifeplus_data_model} object with \code{name = "data_model"}.
#'
#' @noRd
#' @keywords internal
outlier_print_info <- function(x) {
  cli::cli_h3("Settings")
  cli::cli_ul()
  cli::cli_li("Outlier threshold: {.val {x$outlier_threshold}}")
  if (is.null(x$fixed_sd)) {
    cli::cli_li("Prior mean: {.val {x$prior_mean}}")
    cli::cli_li("Prior standard deviation: {.val {x$prior_sd}}")
  } else {
    cli::cli("Fixed standard deviation: {.val {x$fixed_sd}}")
  }
  cli::cli_end()
}

#' Print a Lifeplus Data Model Specification
#'
#' @description
#' Displays a human-readable summary of a data model specification.
#'
#' @param x An object of class \code{lifeplus_data_model}.
#' @param ... Additional arguments (currently ignored).
#'
#' @return Invisibly returns \code{x}
#'
#' @export
print.lifeplus_data_model <- function(x, ...) {
  cli::cli_rule(left = "{.strong lifeplus} data model")
  cli::cli_text("Type: {.emph {x$name}}")

  if (!is.null(x$print_info)) {
    x$print_info(x)
  }

  cli::cli_rule()
  invisible(x)
}
