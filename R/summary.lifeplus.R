#' Compile diagnostics for a lifeplus summary
#' @param object A \code{lifeplus} object.
#' @param estimates Named list of posterior data frames.
#'
#' @return A named list with diagnostic information.
#'
#' @noRd
#' @keywords internal
compile_diagnostics <- function(object, estimates) {
  diag <- object$diagnose

  mcmc <- list(
    n_divergent = sum(diag$num_divergent),
    n_max_treedepth = sum(diag$num_max_treedepth),
    any_low_ebfmi = any(diag$ebfmi < 0.3)
  )

  rhat_threshold <- 1.1

  high_rhat <- lapply(estimates, function(est) {
    if (is.null(est) || !("rhat" %in% colnames(est))) {
      return(NULL)
    }
    est[!is.na(est$rhat) & est$rhat > rhat_threshold, , drop = FALSE]
  })
  high_rhat <- Filter(Negate(is.null), high_rhat)

  list(mcmc = mcmc, high_rhat = high_rhat, rhat_threshold = rhat_threshold)
}

#' Summarise a Lifeplus Model Fit
#' @description
#' Produces a structured summary of posterior estimates from a fitted \code{lifeplus} model.
#'
#' @param object An object of class \code{lifeplus}, as returned by \code{\link{lifeplus}()}.
#' @param component A character string specifying which posterior component to summarise: one of
#'   "life", "shocks", "shockfree", transition_functions", "transition_function_mean", "transition_params", or "all".
#' @param areas An optional character vector of areas names to include in summary.
#' @param times An optional integer vector of time points to include in summary.
#' @param ... Additional arguments (currently ignored).
#'
#' @return An object of class \code{summary.lifeplus}.
#'
#' @seealso \code{link{lifeplus}()}, \code{\link{tidy.lifeplus}()}, \code{\link{print.lifeplus}()}
#'
#' @export
summary.lifeplus <- function(
  object,
  component = c(
    "life",
    "shocks",
    "shockfree",
    "transition_functions",
    "transition_function_mean",
    "transition_params",
    "data_params",
    "data_params",
    "all"
  ),
  areas = NULL,
  times = NULL,
  ...
) {
  acceptable_components <- c(
    "life",
    "shocks",
    "shockfree",
    "transition_functions",
    "transition_function_mean",
    "transition_params",
    "data_params",
    "shock_params",
    "all"
  )
  checkmate::assert_subset(
    component,
    acceptable_components
  )
  component <- match.arg(component)

  components <- component
  if (component == "all") {
    components <- setdiff(acceptable_components, "all")
  }

  estimates <- lapply(components, function(component) {
    est <- object$posteriors[[component]]
    if (is.null(est)) {
      return(NULL)
    }

    if (component %in% c("transition_params", "data_params", "shock_params")) {
      est <- est$params
    }

    if (!is.null(areas) && "area" %in% colnames(est)) {
      unknown <- setdiff(areas, object$areas)
      if (length(unknown) > 0) {
        cli::cli_warn(
          "{length(unknown)} unknown area{?s} will be ignored: {unknown}"
        )
      }
      est <- est[est$area %in% areas, , drop = FALSE]
    }

    if (!is.null(times) && "time" %in% colnames(est)) {
      unknown <- setdiff(times, object$times)
      if (length(unknown) > 0) {
        cli::cli_warn(
          "{length(unknown)} unknown time{?s} will be ignored: {unknown}"
        )
      }
      est <- est[est$time %in% times, , drop = FALSE]
    }
    est
  })
  names(estimates) <- components

  estimates <- Filter(Negate(is.null), estimates)

  if (length(estimates) == 0) {
    cli::cli_warn("No posterior summaries available.")
  }

  diagnostics <- compile_diagnostics(object, estimates)

  if (length(estimates) == 1 && component != "all") {
    estimates <- estimates[[1]]
  }

  structure(
    list(
      component = component,
      estimates = estimates,
      n_areas = length(if (is.null(areas)) object$areas else areas),
      n_times = length(if (is.null(times)) object$times else times),
      posterior_quantiles = object$posterior_quantiles,
      diagnostics = diagnostics
    ),
    class = "summary.lifeplus"
  )
}

#' Print a Lifeplus Model Summary
#'
#' @param x An object of class \code{summary.lifeplus}.
#' @param n Maximum number of rows to print per component (default: \code{10}).
#' @param ... Additional arguments (currently ignored).
#'
#' @return Invisibly returns \code{x}.
#'
#' @seealso \code{\link{summary.lifeplus}()}
#'
#' @export
print.summary.lifeplus <- function(x, n = 10, ...) {
  cli::cli_rule(left = "{.strong lifeplus} summary")

  cli::cli_text(
    "{.val {x$n_areas}} area{?s}, {.val {x$n_times}} time point{?s}"
  )

  if (is.data.frame(x$estimates)) {
    cli::cli_h3("Posterior summaries ({.field {x$component}})")
    print_truncated(x$estimates, n)
  } else {
    for (component in names(x$estimates)) {
      cli::cli_h3("Posterior summaries ({.field {component}})")
      print_truncated(x$estimates[[component]], n)
    }
  }

  cli::cli_h3("Diagnostics")
  mcmc <- x$diagnostics$mcmc
  if (mcmc$n_divergent > 0) {
    cli::cli_alert_danger("{.val {mcmc$n_divergent}} divergent transition{?s}")
  }
  if (mcmc$n_max_treedepth > 0) {
    cli::cli_alert_danger(
      "{.val {mcmc$n_max_treedepth}} transition{?s} hit max treedepth"
    )
  }
  if (mcmc$any_low_ebfmi) {
    cli::cli_alert_danger("Low E-BFMI detected (< 0.3)")
  }

  high_rhat <- x$diagnostics$high_rhat
  n_high_rhat <- sum(vapply(high_rhat, nrow, 1))
  rhat_threshold <- x$diagnostics$rhat_threshold
  if (n_high_rhat > 0) {
    cli::cli_alert_danger(
      "{.val {n_high_rhat}} parameter{?s} with Rhat > {.val {rhat_threshold}}"
    )
  }

  if (
    mcmc$n_divergent == 0 &&
      mcmc$n_max_treedepth == 0 &&
      !mcmc$any_low_ebfmi &&
      n_high_rhat == 0
  ) {
    cli::cli_alert_success("All diagnostics passed")
  }

  cli::cli_rule()
  invisible(x)
}

#' Print a data frame, truncating to n rows
#' @noRd
#' @keywords internal
print_truncated <- function(df, n) {
  if (nrow(df) <= n) {
    print(df)
  } else {
    print(df[seq_len(n), , drop = FALSE])
    cli::cli_text("{.emph ... and {.val {nrow(df) - n}} more row{?s}}")
  }
}
