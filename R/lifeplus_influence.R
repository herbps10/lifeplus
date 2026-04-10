#' Compute posterior covariance influence measure
#'
#' @description
#' Computes the posterior covariance between the log likelihood and the given parameter.
#'
#' @param x An object of class \code{lifeplus}, as returned by \code{\link{lifeplus}()}.
#' @param parameters Vector of model parameter names
#'
#' @return vector of influence measures
#' @importFrom posterior as_draws_matrix
#' @export
lifeplus_influence <- function(x, parameters) {
  checkmate::assert_class(x, "lifeplus")

  log_lik_mat <- posterior::as_draws_matrix(x$log_lik)
  param_draws <- posterior::as_draws_matrix(x$samples$draws(parameters))

  influence <- cov(log_lik_mat, param_draws)

  rownames(influence) <- paste0("obs_", 1:ncol(log_lik_mat))

  influence
}
