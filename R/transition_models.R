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
#' @param degree Spline degree
#' @param num_knots Number of spline knots
#'
#' @return lifeplus_transition_model
#' @export
transition_model_spline <- function(degree = 2, num_knots = 7) {
  checkmate::check_integer(degree, lower = 1)
  checkmate::check_integer(num_knots, lower = 1)

  f <- function(y, grid) {
    knots <- sort(c(seq(0, max(y) / 110, length.out = num_knots), 1, 2))
    B <- t(splines::bs(grid, knots = knots, degree = degree, intercept = FALSE))
    B <- B[1:(nrow(B) - 1), ]
    num_grid <- length(grid)
    num_basis <- nrow(B)
    ext_knots <- c(
      rep(knots[1], degree),
      knots,
      rep(knots[length(knots)], degree)
    )

    list(
      degree = degree,
      num_knots = num_knots,
      B = B,
      knots = knots,
      ext_knots = ext_knots
    )
  }

  x <- list(
    name = "spline",
    stan_data = f
  )

  class(x) <- "lifeplus_transition_model"
  x
}
