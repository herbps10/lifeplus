#' Normal data model
#'
#' @param prior_mean Prior mean for white noise standard deviation
#' @param prior_sd Prior standard deviation for white noise standard deviation
#' @details
#' The residuals are modeled as \eqn{\epsilon_{c,t} \sim N(0, \sigma^2)}.
#' The prior for the residual scale parameter is \eqn{\sigma \sim N(m, s^2)},
#' where \eqn{m} and \eqn{s} are user-specified via the \code{prior_mean} and \code{prior_sd} arguments, respectively.
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
#' The remaining residuals are modeled as in the normal data model:  \eqn{\epsilon_{c,t} \sim N(0, \sigma^2)}.
#' The prior for the residual scale parameter is \eqn{\sigma \sim N(m, s^2)},
#' where \eqn{m} and \eqn{s} are user-specified via the \code{prior_mean} and \code{prior_sd} arguments, respectively.
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
  checkmate::assert_numeric(scale_global, lower = 0)
  checkmate::assert_numeric(slab_scale, lower = 0)
  checkmate::assert_numeric(slab_df, lower = 0)
  checkmate::assert_flag(constrain_negative)

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
