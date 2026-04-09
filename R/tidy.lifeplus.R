#' Tidy a Lifeplus Model Fit
#'
#' @description
#' Extracts posterior summaries from a fitted \code{lifeplus} model as a tidy tibble,
#' following the conventions of the \href{https://broom.tidymodels.org/}{broom} package.
#' Each row represents one parameter (area x time combination) with its posterior median and
#' credible interval.
#'
#' @param x An object of class \code{lifeplus}, as returned by \code{\link{lifeplus}()}.
#' @param conf.level A scalar in \code{(0, 1)} specifying the credible interval width.
#'   Defaults to \code{0.95} (i.e., a 95\% interval).
#' @param component A character string specifying which variable to extract.
#'   One of \code{"eta"} (default), \code{"shocks"}, or \code{"eta_shockfree"}. The latter two
#'   require the presence of a shock model.
#' @param ... Additional arguments (currently ignored)
#'
#' @return A \code{\link[tibble]{tibble}} with columns:
#' \describe{
#'   \item{\code{term}}{Parameter name.}
#'   \item{\code{area}}{Area identifier.}
#'   \item{\code{time}}{Time point.}
#'   \item{\code{estimate}}{Posterior median.}
#'   \item{\code{conf.low}}{Lower bound of the credible interval.}
#'   \item{\code{conf.high}}{Upper bound of the credible interval.}
#'   \item{\code{conf.level}}{Credible interval level.}
#' }
#'
#' @seealso \code{\link{lifeplus}()}, \code{\link{summary.lifeplus}()}
#'
#' @importFrom generics tidy
#' @importFrom tibble tibble
#' @export
tidy.lifeplus <- function(
  x,
  conf.level = 0.95,
  component = c("life", "shocks", "shockfree")
) {
  component <- match.arg(component)
  checkmate::assert_number(conf.level, lower = 0, upper = 1)

  if (component %in% c("shocks", "shockfree") && x$shock_model$name == "none") {
    cli::cli_abort(
      "Model was fit without a shocks model; no posterior shocks to retrieve."
    )
  }

  # Extract the posterior
  post <- x$posteriors[[component]]

  quantile_cols <- resolve_conf_columns(x$posterior_quantiles, conf.level)

  result <- tibble::tibble(
    term = component,
    area = post$area,
    time = post$time,
    estimate = post[["50%"]],
    conf.low = post[[quantile_cols$low]],
    conf.high = post[[quantile_cols$high]],
    conf.level = conf.level
  )

  result
}

resolve_conf_columns <- function(posterior_quantiles, conf.level) {
  low_p <- round((1 - conf.level) / 2, 5)
  high_p <- round(1 - low_p, 5)

  if (!(low_p %in% posterior_quantiles) || !(high_p %in% posterior_quantiles)) {
    cli::cli_abort(
      "Quantiles for {.val {conf.level * 100}}% interval not available."
    )
  }

  list(
    low = paste0(low_p * 100, "%"),
    high = paste0(high_p * 100, "%")
  )
}
