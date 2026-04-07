load_model <- function(transition_model, data_model, shock_model) {
  model_name <- paste0(paste0(
    c(transition_model$name, data_model$name, shock_model$name),
    collapse = "_"
  ))

  return(instantiate::stan_package_model(
    model_name,
    package = "lifeplus"
  ))
}


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

#' Normal data model
#'
#' @param prior_mean Prior mean for white noise standard deviation
#' @param prior_sd Prior standard deviation for white noise standard deviation
#' @details
#' The residuals are modeled as $\epsilon_{c,t} \sim N(0, \sigma^2)$.
#' The prior for the residual scale parameter is $\sigma \sim N(m, s^2)$,
#' where $m$ and $s$ are user-specified via the `prior_mean` and `prior_sd` arguments, respectively.
#'
#' @return lifeplus_data_model
#' @export
data_model_normal <- function(prior_mean = 0, prior_sd = 1) {
  x <- list(
    name = "normal",
    stan_data = list(
      epsilon_sigma_prior_mu = prior_mean,
      epsilon_sigma_prior_sd = prior_sd
    )
  )
  class(x) <- "lifeplus_data_model"
  x
}


#' Outlier data model
#'
#' @param outlier_threshold threshold for ignoring observed life expectancy differences.
#' @param prior_mean Prior mean for white noise standard deviation
#' @param prior_sd Prior standard deviation for white noise standard deviation
#' @details
#' Any observed absolute differences in life expectancy greater than the user-specified outlier threshold (`outlier_threshold`) are ignored.
#' The remaining residuals are modeled as in the normal data model:  $\epsilon_{c,t} \sim N(0, \sigma^2)$.
#' The prior for the residual scale parameter is $\sigma \sim N(m, s^2)$,
#' where $m$ and $s$ are user-specified via the `prior_mean` and `prior_sd` arguments, respectively.
#'
#' @return lifeplus_data_model
#' @export
data_model_outlier <- function(
  outlier_threshold = 5,
  prior_mean = 0,
  prior_sd = 1
) {
  x <- list(
    name = "outlier",
    stan_data = list(
      outlier_threshold = outlier_threshold,
      epsilon_sigma_prior_mu = prior_mean,
      epsilon_sigma_prior_sd = prior_sd
    )
  )
  class(x) <- "lifeplus_data_model"
  x
}

#' Regularized horseshoe prior model for shocks
#' @param scale_global Global scale parameter
#' @param slab_scale Slab scale parameter
#' @param slab_df Slab degrees of freedom parameter
#' @param constrain_negative Whether to constrain shocks to be negative (default: FALSE)
#' @details
#' The regularized horseshoe prior aggressively shrinks the shock term to zero.
#'
#'
#'
#' Reference: Juho Piironen, Aki Vehtari "Sparsity information and regularization in the horseshoe and other shrinkage priors,"
#' Electronic Journal of Statistics, Electron. J. Statist. 11(2), 5018-5051, (2017)
#' @export
shock_model_regularized_horseshoe <- function(
  scale_global = 0.1,
  slab_scale = 10,
  slab_df = 6,
  constrain_negative = FALSE
) {
  checkmate::check_numeric(scale_global, lower = 0)
  checkmate::check_numeric(slab_scale, lower = 0)
  checkmate::check_numeric(slab_df, lower = 0)
  checkmate::check_flag(constrain_negative)

  x <- list(
    name = "regularized_horseshoe",
    stan_data = list(
      scale_global = scale_global,
      slab_scale = slab_scale,
      slab_df = slab_df,
      constrain_negative = as.numeric(constrain_negative)
    )
  )
  class(x) <- "lifeplus_shock_model"
  x
}

#' Fix all shock terms to zero.
#'
#' @export
shock_model_none <- function() {
  x <- list(
    name = "none"
  )
  class(x) <- "lifeplus_shock_model"
  x
}
