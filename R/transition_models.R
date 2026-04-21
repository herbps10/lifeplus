##### Double Logistic ######

#' Double Logistic Transition Model
#'
#' @description
#' Specifies a double logistic transition model for use in \code{\link{lifeplus}()}.
#' The double logistic function models the relationship between age and the rate of
#' change of mortality with six parameters.
#'
#' @details
#' The six parameters of the double logistic function are:
#' \describe{
#'   \item{\code{Delta1}}{Shape parameter constrained in \code{(0, 50)}.}
#'   \item{\code{Delta2}}{Shape parameters constrained in \code{(0, 50)}.}
#'   \item{\code{Delta3}}{Shape parameters constrained in \code{(0, 50)}.}
#'   \item{\code{Delta4}}{Shape parameters constrained in \code{(5, 50)}.}
#'   \item{\code{k}}{Shape parameter constrained to \code{(0, 10)}.}
#'   \item{\code{z}}{Asymptotic growth rate parameter constrained to \code{(0, 1.15/5)}.}
#' }
#'
#' When \code{hierarchical = TRUE} (the default), the parameters are modeled with
#' a multivariate normal hierarchical prior across areas, allowing partial pooling.
#' When \code{FALSE}, each area's parameters are estimated independently.
#'
#' @param hierarchical Logical; whether to estimate the double logistic parameters
#'   hierarchically across areas (default: \code{TRUE}). Setting this to \code{TRUE}
#'   enables partial pooling, which can improve estimation for areas with sparse data.
#' @param multivariate_prior Logical; whether to apply a multivariate normal prior with
#'   hierarchical estimation.
#' @param log_scale Logical; whether parameters are estimated on log scale (default: \code{FALSE}).
#' @param prior_mean Prior mean for Delta1--Delta4, k, and z parameters (default: 0).
#' @param prior_sd Prior standard deviation for Delta1--Delta4, k, and z parameters (default: 0).
#'
#' @return An object of class \code{lifeplus_transition_model}, a named list containing:
#' \describe{
#'   \item{\code{name}}{Character string \code{"double_logistic"}.}
#'   \item{\code{hierarchical}}{Logical; whether hierarchical estimation is enabled.}
#'   \item{\code{multivariate_prior}}{Logical; whether multivariate normal prior is enabled.}
#'   \item{\code{log_scale}}{Logical; whether to model hierarchical standard deviations on log-scale.}
#'   \item{\code{prior_mean}}{Prior mean.}
#'   \item{\code{prior_sd}}{Prior standard deviation.}
#'   \item{\code{param_names}}{Character vector of parameter names for labelling posterior summaries.}
#'   \item{\code{D_phi}}{Dimension of shared parameters for expectation propagation.}
#'   \item{\code{stan_data}}{A named list of data passed to the Stan model.}
#'   \item{\code{extract_params}}{A function for extracting transition parameters from a fitted model
#'     (used internally)}
#'   \item{\code{print_info}}{A function for printing the transition model.}
#' }
#'
#' @examples
#' # Default: hierarchical estimation
#' transition_model_double_logistic()
#'
#' # Independent estimation per area
#' transition_model_double_logistic(hierarchical = FALSE)
#'
#' @seealso
#' \code{\link{lifeplus}()} for model fitting.
#'
#' @export
transition_model_double_logistic <- function(
  hierarchical = TRUE,
  multivariate_prior = TRUE,
  log_scale = FALSE,
  prior_mean = 0,
  prior_sd = 1
) {
  checkmate::assert_flag(hierarchical)

  n_params <- 6L

  param_names <- c("Delta1", "Delta2", "Delta3", "Delta4", "k", "z")

  if (log_scale == FALSE) {
    Delta_sigma_prior_mean <- rep(0, n_params)
    Delta_sigma_prior_sd <- rep(1, n_params)
  } else {
    Delta_sigma_prior_mean <- rep(-2, n_params)
    Delta_sigma_prior_sd <- rep(1, n_params)
  }

  structure(
    list(
      name = "double_logistic",
      hierarchical = hierarchical,
      param_names = param_names,
      multivariate_prior = multivariate_prior,
      log_scale = log_scale,
      D_phi = ifelse(hierarchical == 1, n_params * 2, 0),
      prior_mean = Delta_prior_mean,
      prior_sd = Delta_prior_sd,
      stan_data = list(
        D = n_params,
        hierarchical = as.integer(hierarchical),
        Delta_constrain = rep(1, n_params),
        Delta_lower = c(0, 0, 0, 5, 0, 0),
        Delta_upper = c(50, 50, 50, 50, 10, 1.15 / 5),
        Delta_prior_mean = rep(prior_mean, n_params),
        Delta_prior_sd = rep(prior_sd, n_params),
        Delta_sigma_prior_mean = Delta_sigma_prior_mean,
        Delta_sigma_prior_sd = Delta_sigma_prior_sd,
        Delta_sigma_lower = rep(0, n_params),
        Delta_multi = as.numeric(multivariate_prior),
        Delta_log_scale = as.integer(log_scale)
      ),
      extract_params = double_logistic_extract_params,
      print_info = double_logistic_print_info
    ),
    class = "lifeplus_transition_model"
  )
}

#' Extract posterior summaries for double logistic parameters
#'
#' @noRd
#' @keywords internal
double_logistic_extract_params <- function(fit, n_cores, posterior_quantiles) {
  param_names <- fit$transition_model$param_names

  post <- summarise_draws(fit, "Delta", n_cores, posterior_quantiles)

  parsed <- parse_stan_indices(post$variable)
  post$variable <- param_names[parsed$indices[, 2]]
  post$area <- fit$areas[parsed$indices[, 1]]

  correlation <- NULL
  if (
    is_hierarchical_transition(fit) &&
      isTRUE(as.logical(fit$transition_model$multivariate_prior == TRUE))
  ) {
    correlation <- summarise_draws(
      fit,
      "Omega_Delta",
      n_cores,
      posterior_quantiles
    )

    parsed_corr <- parse_stan_indices(correlation$variable)
    correlation$variable1 <- param_names[parsed_corr$indices[, 1]]
    correlation$variable2 <- param_names[parsed_corr$indices[, 2]]
  }

  list(params = post, correlation = correlation)
}

#' Print details for the double logistic transition model
#'
#' @param x A \code{lifeplus_transition_model} object with \code{name = "double_logistic"}.
#'
#' @noRd
#' @keywords internal
double_logistic_print_info <- function(x) {
  if (x$hierarchical) {
    cli::cli_alert_success("Hierarchical estimation: enabled")
  } else {
    cli::cli_alert_info("Hierarchical estimation: disabled")
  }
}

##### Spline #####

#' Spline Transition Model
#'
#' @param hierarchical Whether to estimate hierarchical Gaussian Process
#' @param degree Degree of B-spline basis functions
#' @param knots Number of B-spline knots (evenly spaced)

#' @return An object of class \code{lifeplus_transition_model} containing:
#' \describe{
#'   \item{\code{name}}{Character string \code{"gaussian_process"}}
#'   \item{\code{hierarchical}}{Logical; whether hierarchical estimation is enabled.}
#'   \item{\code{degree}}{Integer; degree of B-spline basis functions.}
#'   \item{\code{knots}}{Integer; number of B-spline knots}
#'   \item{\code{stan_data}}{A function used internally to form the data passed to Stan for the spline model}
#'   \item{\code{extract_params}}{A function for extracting spline parameters from a fitted model (used internally).}
#'   \item{\code{print_info}}{A function for printing the transition model.}
#' }
#'
#' @export
transition_model_spline <- function(
  hierarchical = TRUE,
  degree = 2,
  knots = 7
) {
  checkmate::assert_flag(hierarchical)
  checkmate::assert_count(degree, positive = TRUE)
  checkmate::assert_count(knots, positive = TRUE)

  degree <- as.integer(degree)
  knots <- as.integer(knots)

  stan_data_fn <- function(y, grid) {
    build_spline_stan_data(y, grid, degree, knots, hierarchical)
  }

  structure(
    list(
      name = "spline",
      hierarchical = hierarchical,
      degree = degree,
      knots = knots,
      stan_data = stan_data_fn,
      extract_params = spline_extract_params,
      print_info = spline_print_info
    ),
    class = "lifeplus_transition_model"
  )
}

build_spline_stan_data <- function(y, grid, degree, knots, hierarchical) {
  knots_vec <- sort(c(seq(0, max(y) / 110, length.out = knots), 1, 2))
  B <- t(splines::bs(
    grid,
    knots = knots_vec,
    degree = degree,
    intercept = FALSE
  ))
  B <- B[1:(nrow(B) - 1), ]
  num_grid <- length(grid)
  num_basis <- nrow(B)
  ext_knots <- c(
    rep(knots_vec[1], degree),
    knots_vec,
    rep(knots[length(knots)], degree)
  )

  list(
    spline_degree = degree,
    num_basis = num_basis,
    num_knots = length(knots_vec),
    B = B,
    knots = knots_vec,
    ext_knots = ext_knots,

    hierarchical = as.numeric(hierarchical),
    alpha_constrain = rep(c(1), num_basis),
    alpha_lower = rep(c(0), num_basis),
    alpha_upper = rep(c(10), num_basis),
    alpha_prior_mean = rep(c(-2), num_basis),
    alpha_prior_sd = rep(c(2), num_basis)
  )
}

#' Extract posterior summaries for spline process model parameters
#' @noRd
#' @keywords internal
spline_extract_params <- function(fit, n_cores, posterior_quantiles) {
  vars <- c("alpha")
  post <- summarise_draws(fit, vars, n_cores, posterior_quantiles)

  parsed <- parse_stan_indices(post$variable)
  post$variable <- parsed$name
  post$area <- fit$areas[parsed$indices[, 1]]
  post$basis <- parsed$indices[, 2]

  list(params = post, correlation = NULL)
}

#' Print details for the spline transition model
#'
#' @param x A \code{lifeplus_transition_model} object with \code{name = "spline"}.
#'
#' @noRd
#' @keywords internal
spline_print_info <- function(x) {
  if (x$hierarchical) {
    cli::cli_alert_success("Hierarchical estimation: enabled")
  } else {
    cli::cli_alert_info("Hierarchical estimation: disabled")
  }

  cli::cli_h3("Approximation settings")
  cli::cli_ul()
  cli::cli_li("Spline degree: {.val {x$degree}}")
  cli::cli_li("Spline knots: {.val {x$knots}}")
  cli::cli_end()
}

##### Gaussian Process ######

#' Approximate Gaussian Process Transition Model
#'
#' @description
#' Specifies an approximate Gaussian Process (GP) transition model for use in
#' \code{\link{lifeplus}()}. The transition function relating age to
#' the rate of change of mortality is given a GP prior, providing a flexible nonparametric
#' alternative to parameteric transition models (e.g., \code{\link{transition_model_double_logistic}()}).
#'
#' @details
#' For computational efficiency, the GP is approximated using a Hilbert space
#' basis function expansion (Riutort-Mayol et al. 2023). The two tuning parameters
#' control the fidelity of this approximation:
#' \describe{
#'   \item{\code{basis_functions}}{The number of basis fucntions \eqn{M}.
#'     Higher values give a better approximation of the true GP but increase
#'     the number of parameters. Values in the range 15-30 generally suffice.}
#'   \item{\code{boundary_factor}}{A multiplier \eqn{L} that defines the boundary
#'     of the approximation domain as \eqn{L} times the range of the input (age) grid.
#'     Values in the range 1.2-2.0 are typical.}
#' }
#'
#' When \code{hierarchical = TRUE} (default), the GP basis weights are modeled hierarchically
#' across areas, enabling partial pooling. When \code{FALSE}, basis weights are estimated independently
#' across areas.
#'
#' @param hierarchical Whether to estimate hierarchical Gaussian Process
#' @param basis_functions Number of basis functions to use in approximation
#' @param boundary_factor Approximation boundary factor tuning parameter
#'
#' @references
#' Riutort-Mayol, G., Bürkner, PC., Andersen, M.R., Solin. A, and Vehtari A. (2023)
#' Practical Hilbert space approximate Bayesian Gaussian processes for probabilistic programming.
#' \emph{Statistics and Computing} 33, 17. \doi{10.1007/s11222-022-10167-2}
#'
#' @return An object of class \code{lifeplus_transition_model} containing:
#' \describe{
#'   \item{\code{name}}{Character string \code{"gaussian_process"}}
#'   \item{\code{hierarchical}}{Logical; whether hierarchical estimation is enabled.}
#'   \item{\code{basis_functions}}{Integer; number of basis functions used.}
#'   \item{\code{boundary_factor}}{Numeric; boundary expansion factor used.}
#'   \item{\code{stan_data}}{A named list of data passed to the Stan model.}
#'   \item{\code{extract_params}}{A function for extracting GP parameters from a fitted model (used internally).}
#'   \item{\code{print_info}}{A function for printing the transition model.}
#' }
#'
#' @examples
#' # Default settings
#' transition_model_gaussian_process()
#'
#' @seealso
#' \code{\link{lifeplus}()} for model fitting.
#' \code{\link{transition_model_double_logistic}()} for a parametric alternative.
#'
#' @export
transition_model_gaussian_process <- function(
  hierarchical = TRUE,
  basis_functions = 20,
  boundary_factor = 1.5
) {
  checkmate::assert_flag(hierarchical)
  checkmate::assert_count(basis_functions, positive = TRUE)
  checkmate::assert_number(boundary_factor, lower = 1)

  basis_functions <- as.integer(basis_functions)

  structure(
    list(
      name = "gaussian_process",
      hierarchical = hierarchical,
      basis_functions = basis_functions,
      boundary_factor = boundary_factor,
      stan_data = list(
        hierarchical = as.integer(hierarchical),
        L = boundary_factor,
        M = basis_functions,
        beta_constrain = rep(c(0), basis_functions),
        beta_lower = rep(c(0), basis_functions),
        beta_upper = rep(c(1), basis_functions),
        beta_prior_mean = rep(c(0), basis_functions),
        beta_prior_sd = rep(c(1), basis_functions)
      ),
      extract_params = gaussian_process_extract_params,
      print_info = gaussian_process_print_info
    ),
    class = "lifeplus_transition_model"
  )
}

#' Extract posterior summaries for Gaussian Process parameters
#' @noRd
#' @keywords internal
gaussian_process_extract_params <- function(fit, n_cores, posterior_quantiles) {
  vars <- c("gp_sigma", "gp_lengthscale")
  post <- summarise_draws(fit, vars, n_cores, posterior_quantiles)

  list(params = post, correlation = NULL)
}

#' Print details for the Gaussian Process transition model
#'
#' @param x A \code{lifeplus_transition_model} object with \code{name = "gaussian_process"}.
#'
#' @noRd
#' @keywords internal
gaussian_process_print_info <- function(x) {
  if (x$hierarchical) {
    cli::cli_alert_success("Hierarchical estimation: enabled")
  } else {
    cli::cli_alert_info("Hierarchical estimation: disabled")
  }

  cli::cli_h3("Approximation settings")
  cli::cli_ul()
  cli::cli_li("Basis functions: {.val {x$basis_functions}}")
  cli::cli_li("Boundary factor: {.val {x$boundary_factor}}")
  cli::cli_end()
}


#' Print a Lifeplus Transition Model Specification
#'
#' @description
#' Displays a human-readable summary of a transition model specification.
#'
#' @param x An object of class \code{lifeplus_transition_model}.
#' @param ... Additional arguments (currently ignored).
#'
#' @return Invisibly returns \code{x}
#'
#' @export
print.lifeplus_transition_model <- function(x, ...) {
  cli::cli_rule(left = "{.strong lifeplus} transition model")
  cli::cli_text("Type: {.emph {x$name}}")

  if (!is.null(x$print_info)) {
    x$print_info(x)
  }

  cli::cli_rule()
  invisible(x)
}
