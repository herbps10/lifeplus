#' Cavity prior specification for expectation propagation.
#'
#' @param mu numeric vector specifying cavity prior mean
#' @param Sigma numeric matrix specifying cavity prior covariance matrix
#'
#' @return object of class \code{lifeplus_cavity_prior}
#' @export
lifeplus_cavity_prior <- function(mu, Sigma) {
  structure(
    list(mu = mu, Sigma = Sigma),
    class = "lifeplus_cavity_prior"
  )
}
