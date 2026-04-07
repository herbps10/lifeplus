#' Double logistic transition model
#'
#' @param hierarchical whether to estimate double logistic parameters hierarchically (default: TRUE)
#'
#' @return Object of class lifeplus_transition_model
#' @export
transition_model_double_logistic <- function(hierarchical = TRUE) {
  x <- list(
    name = "double_logistic",
    stan_data = list(
      D = 6,
      hierarchical = as.numeric(hierarchical),
      Delta_constrain = c(1, 1, 1, 1, 1, 1),
      Delta_lower = c(0, 0, 0, 5, 0, 0),
      Delta_upper = c(50, 50, 50, 50, 10, 1.15 / 5),
      Delta_prior_mean = c(0, 0, 0, 0, 0, 0),
      Delta_prior_sd = c(1, 1, 1, 1, 1, 1),
      Delta_sigma_lower = c(0, 0, 0, 0, 0, 0)
    )
  )
  class(x) <- "lifeplus_transition_model"
  x
}

#' Spline transition model
#'
#' @param hierarchical Whether to apply hierarchical models to spline coefficients
#' @param degree Spline degree
#' @param knots Number of spline knots
#'
#' @return lifeplus_transition_model
#' @export
transition_model_spline <- function(
  hierarchical = TRUE,
  degree = 2,
  knots = 7
) {
  checkmate::check_flag(hierarchical)
  checkmate::check_integer(degree, lower = 1)
  checkmate::check_integer(knots, lower = 1)

  f <- function(y, grid) {
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
      knots,
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

  x <- list(
    name = "spline",
    stan_data = f
  )

  class(x) <- "lifeplus_transition_model"
  x
}

#' Approximate Gaussian Process transition model
#'
#' @param hierarchical Whether to estimate hierarchical Gaussian Process
#' @param basis_functions Number of basis functions to use in approximation
#' @param boundary_factor Approximation boundary factor tuning parameter
#'
#' @details
#' The transition function is modeled by a Gaussian Process prior.
#' For computational reasons, the Gaussian Process is approximated with a basis function approach
#' (see Riutort-Mayol et al. 2020).
#'
#' References: Riutort-Mayol, G., Bürkner, PC., Andersen, M.R. et al. (2020) Practical Hilbert space approximate Bayesian Gaussian processes for probabilistic programming. Stat Comput 33, 17
#'
#' @return lifeplus_transition_model
#' @export
transition_model_gaussian_process <- function(
  hierarchical = TRUE,
  basis_functions = 20,
  boundary_factor = 1.5
) {
  checkmate::check_flag(hierarchical)
  checkmate::check_integer(basis_functions, lower = 1)
  checkmate::check_integer(boundary_factor, lower = 1)

  x <- list(
    name = "gaussian_process",
    stan_data = list(
      hierarchical = as.integer(hierarchical),
      L = boundary_factor,
      M = basis_functions,
      beta_constrain = rep(c(0), basis_functions),
      beta_lower = rep(c(0), basis_functions),
      beta_upper = rep(c(1), basis_functions),
      beta_prior_mean = rep(c(0), basis_functions),
      beta_prior_sd = rep(c(1), basis_functions)
    )
  )

  class(x) <- "lifeplus_transition_model"
  x
}
