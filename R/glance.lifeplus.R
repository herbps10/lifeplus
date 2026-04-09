#' Glance at a Lifeplus Model Fit
#'
#' @description
#' Extracts principle diagnostic and validation results from a Lifeplus model fit
#' following the conventions of the \code{broom} package.
#' Each row represents one fit. If the model fit includes held-out observations,
#' model validation summaries are extracted based on the most recent held-out observations in each area.
#'
#' @param x An object of class \code{lifeplus}, as returned by \code{\link{lifeplus}()}.
#' @param component A character string specifying which variable to use for validations.
#'   One of \code{"life"} (default), \code{"shockfree"}. The latter
#'   requires the presence of a shock model.
#' @param conf.level A scalar in \code{(0, 1)} specifying the credible interval width.
#'   Defaults to \code{0.95} (i.e., a 95% interval).
#' @param diagnostics Logical; whether to include MCMC diagnostic summary measures (default: \code{TRUE})
#' @param validations Logical; whether to include held-out data validation summary measures (default: \code{TRUE})
#' @param ... Additional arguments (currently ignored)
#'
#' @return A \code{\link[tibble]{tibble}} with columns:
#' \describe{
#'   \item{\code{transition.model}}{Transition model used.}
#'   \item{\code{shock.model}}{Shock model used.}
#'   \item{\code{data.model}}{Data model used.}
#'   \item{\code{num_divergent}}{Number of divergent transitions (total across all chains).}
#'   \item{\code{num_max_treedepth}}{Number of transitions that hit max treedepth (total across all chains).}
#'   \item{\code{min_ebfmi}}{Minimum E-BFMI across all chains.}
#'   \item{\code{mean.error}}{Mean error.}
#'   \item{\code{median.error}}{Median error.}
#'   \item{\code{mean.absolute.error}}{Mean absolute error.}
#'   \item{\code{median.absolute.error}}{Median absolute error.}
#'   \item{\code{ci.width}}{Width of credible interval of level specified by \code{conf.level}.}
#'   \item{\code{coverage}}{Empirical coverage of the credible interval of level specified by \code{conf.level}.}
#'   \item{\code{conf.level}}{Credible interval level specified by \code{conf.level}}
#' }
#'
#' @seealso \code{\link{lifeplus}()}, \code{\link{summary.lifeplus}()}
#'
#' @importFrom generics glance
#' @importFrom tibble tibble
#' @importFrom stats median
#'
#' @export
glance.lifeplus <- function(
  x,
  component = c("life", "shockfree"),
  conf.level = 0.95,
  diagnostics = TRUE,
  validations = TRUE,
  ...
) {
  component <- match.arg(component)
  checkmate::assert_number(conf.level, lower = 0, upper = 1)

  if (component %in% c("shockfree") && x$shock_model$name == "none") {
    cli::cli_abort(
      "Model was fit without a shocks model; no posterior shocks to retrieve."
    )
  }

  # Extract the posterior
  cutoff_time <- x$cutoff_time
  n_held_out <- sum(x$data[[x$time]] > cutoff_time)

  result <- tibble::tibble(
    transition.model = x$transition_model$name,
    shock.model = x$shock_model$name,
    data.model = x$data_model$name,
  )

  if (diagnostics) {
    result$num_divergent <- sum(x$diagnose$num_divergent)
    result$num_max_treedepth <- sum(x$diagnose$num_max_treedepth)
    result$min_ebfmi <- min(x$diagnose$ebfmi)
  }

  if (n_held_out > 0 && validations) {
    valid <- x$validation[[component]]
    valid <- valid[valid$time == max(valid$time), ]

    result$mean_error <- mean(valid$error)
    result$median_error <- stats::median(valid$error)
    result$mean_absolute_error <- mean(abs(valid$error))
    result$median_absolute_error <- stats::median(abs(valid$error))
    result$ci.width <- mean(valid[[paste0("ci_width", conf.level * 100, "%")]])
    result$coverage <- mean(valid[[paste0("covered", conf.level * 100, "%")]])
    result$conf.level <- conf.level
  }

  result
}
