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
