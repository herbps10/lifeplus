#' Regularized horseshoe prior model for shock terms.
#'
#' @param scale_global Global scale parameter
#' @param slab_scale Slab scale parameter
#' @param slab_df Slab degrees of freedom parameter
#' @param constrain_negative Whether to constrain shocks to be negative (default: FALSE)
#'
#' @details
#' The regularized horseshoe prior aggressively shrinks the shock term to zero.
#'
#' @references
#' Juho Piironen, Aki Vehtari "Sparsity information and regularization in the horseshoe and other shrinkage priors,"
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

  structure(
    list(
      name = "regularized_horseshoe",
      scale_global = scale_global,
      slab_scale = slab_scale,
      slab_df = slab_df,
      stan_data = list(
        scale_global = scale_global,
        slab_scale = slab_scale,
        slab_df = slab_df,
        constrain_negative = as.numeric(constrain_negative)
      ),
      print_info = regularized_horseshoe_print_info
    ),
    class = "lifeplus_shock_model"
  )
}

#' Print details for the regularized horseshoe shock model
#'
#' @param x A \code{lifeplus_shock_model} object with \code{name = "regularized_horseshoe"}.
#'
#' @noRd
#' @keywords internal
regularized_horseshoe_print_info <- function(x) {
  cli::cli_h3("Prior settings")
  cli::cli_ul()
  cli::cli_li("Global scale: {.val {x$scale_global}}")
  cli::cli_li("Slab scale: {.val {x$slab_scale}}")
  cli::cli_li("Slab degrees of freedom: {.val {x$slab_df}}")
  cli::cli_end()
}

#' No shock model
#'
#' @description
#' Fixes all shock terms to zero.
#'
#' @export
shock_model_none <- function() {
  structure(
    list(
      name = "none"
    ),
    class = "lifeplus_shock_model"
  )
}


#' Print a Lifeplus Shock Model Specification
#'
#' @description
#' Displays a human-readable summary of a shock model specification.
#'
#' @param x An object of class \code{lifeplus_transition_model}.
#' @param ... Additional arguments (currently ignored).
#'
#' @return Invisibly returns \code{x}
#'
#' @export
print.lifeplus_shock_model <- function(x, ...) {
  cli::cli_rule(left = "{.strong lifeplus} shock model")
  cli::cli_text("Type: {.emph {x$name}}")

  if (!is.null(x$print_info)) {
    x$print_info(x)
  }

  cli::cli_rule()
  invisible(x)
}
