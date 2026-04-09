#' Regularized Horseshoe Shock Model
#'
#' @description
#' Specifies a regularized horseshoe prior for the shock terms in a \code{\link{lifeplus}()} model.
#' This prior aggressively shrinks most shocks toward zero while allowing a small number of
#' large shocks to escape shrinkage.
#'
#' @details
#' The regularized horseshoe prior (Piironen and Vehtari, 2017) is parameterized by several
#' prior settings:
#' \describe{
#'   \item{scale_global}{Controls the overall level of shrinkage.
#'     Smaller values impose stronger sparsity (fewer shocks escape shrinkage). To improve
#'     MCMC efficiency, the global is fixed to this value rather than estimated from the data.}
#'   \item{slab_scale}{Regularizes the magnitude of shocks that escape the global shrinkage.
#'     Larger values allow larger shocks, and smaller values pull even the ``active'' shocks toward zero.}
#'   \item{slab_df}{Degrees of freedom for the slab (regularization) component.
#'     Smaller values allow for heavier tails, and larger values converge toward a Gaussian prior for the slab.}
#'   \item{local_df}{Degrees of freedom for the half-\eqn{t} prior on the local scale parameters.}
#' }
#'
#' Setting \code{constrain_negative = TRUE} constrains all shock terms to be non-positive, which is appropriate
#' when shocks are expected to represent mortality crises (reductions in life expectancy) only.
#'
#' @param scale_global A positive number giving the fixed global scale parameter (default: \code{0.1}).
#' @param slab_scale A positive number giving the slab scale parameter (default: \code{10}).
#' @param slab_df A positive number (>= 1) giving the slab degrees of freedom (default: \code{6}).
#' @param local_df A positive number (>= 1) giving the degrees of freedom for the half-\eqn{t} prior on the local
#'   scale parameters (default: \code{3}).
#' @param constrain_negative Logical; whether to constrain shocks to be non-positive (default: \code{FALSE})
#'
#' @return object of type \code{lifeplus_shock_model}, a named list constaining:
#' \describe{
#'   \item{\code{name}}{Character string \code{"regularized_horseshoe"}}
#'   \item{\code{scale_global}}{The global scale value used.}
#'   \item{\code{slab_scale}}{The slab scale value used.}
#'   \item{\code{slab_df}}{The slab degrees of freedom used.}
#'   \item{\code{local_df}}{The local degrees of freedom used.}
#'   \item{\code{constrain_negative}}{Logical; whether shocks are constrained to be non-positive.}
#'   \item{\code{stan_data}}{A named list of data passed to the Stan model.}
#'   \item{\code{print_info}}{A function for printing model-specific details.}
#' }
#'
#' @examples
#' # Default settings
#' shock_model_regularized_horseshoe()
#'
#' # Stronger sparsity with constrained shocks
#' shock_model_regularized_horseshoe(
#'   scale_global = 0.05,
#'   constrain_negative = TRUE
#' )
#'
#'
#' @references
#' Juho Piironen, Aki Vehtari "Sparsity information and regularization in the horseshoe and other shrinkage priors,"
#' Electronic Journal of Statistics, Electron. J. Statist. 11(2), 5018-5051, (2017)
#' @export
shock_model_regularized_horseshoe <- function(
  scale_global = 0.1,
  slab_scale = 10,
  slab_df = 6,
  local_df = 3,
  constrain_negative = FALSE
) {
  checkmate::assert_number(scale_global, lower = 0)
  checkmate::assert_number(slab_scale, lower = 0)
  checkmate::assert_number(slab_df, lower = 0)
  checkmate::assert_number(local_df, lower = 0)
  checkmate::assert_flag(constrain_negative)

  structure(
    list(
      name = "regularized_horseshoe",
      scale_global = scale_global,
      slab_scale = slab_scale,
      slab_df = slab_df,
      local_df = local_df,
      constrain_negative = constrain_negative,
      stan_data = list(
        scale_global = scale_global,
        slab_scale = slab_scale,
        slab_df = slab_df,
        nu_local = local_df,
        constrain_negative = as.integer(constrain_negative)
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
  cli::cli_li("Local parameter degrees of freedom: {.val {x$local_df}}")
  cli::cli_end()
}

#' No Shock Model
#'
#' @description
#' Fixes all shock terms to zero.
#'
#' #' @return object of type \code{lifeplus_shock_model}, a named list constaining:
#' \describe{
#'   \item{\code{name}}{Character string \code{"none"}}
#' }
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
