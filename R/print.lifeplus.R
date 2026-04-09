#' Print a Lifeplus Model Fit
#'
#' @description
#' Displays a concise summary of a fitted \code{lifeplus} model.
#'
#' @param x An object of class \code{lifeplus}, as returned by \code{\link{lifeplus}()}.
#' @param ... Additional arguments passed to other methods (currently ignored).
#'
#' @return Invisibly returns \code{x}.
#'
#' @seealso \code{\link{lifeplus}()}
#'
#' @export
print.lifeplus <- function(x, ...) {
  cli::cli_rule(left = "{.strong lifeplus} model fit")

  ##### Inputs #####
  n_areas <- length(x$areas)
  n_obs <- nrow(x$data)
  cutoff_time <- x$cutoff_time
  n_held_out <- sum(x$data[[x$time]] > cutoff_time)
  start_time <- min(x$times)
  end_time <- max(x$times)
  n_time <- length(x$times)

  cli::cli_h3("Data")
  cli::cli_ul()
  cli::cli_li("Outcome variable: {.field {x$y}}")
  if (n_held_out == 0) {
    cli::cli_li(
      "Time variable: {.field {x$time}} ({.val {start_time}} to {.val {end_time}})"
    )
  } else {
    cli::cli_li(
      "Time variable: {.field {x$time}} ({.val {start_time}} to {.val {end_time}}, all data after {.val {cutoff_time}} held out)"
    )
  }
  cli::cli_li(
    "Area variable: {.field {x$area}} ({.val {n_areas}} areas)"
  )
  cli::cli_li("Observations: {.val {n_obs}} ({.val {n_held_out}} held out)")
  cli::cli_end()

  ##### Model components ######
  cli::cli_h3("Model components")
  cli::cli_ul()
  cli::cli_li("Transition: {.emph {x$transition_model$name}}")
  cli::cli_li("Data: {.emph {x$data_model$name}}")
  cli::cli_li("Shock: {.emph {x$shock_model$name}}")
  cli::cli_end()

  ##### Sampling #####
  cli::cli_h3("Sampling")
  n_chains <- length(x$samples$metadata()$id)
  n_iter <- x$samples$metadata()$iter_sampling
  n_warmup <- x$samples$metadata()$iter_warmup

  cli::cli_ul()
  cli::cli_li("Chains: {.val {n_chains}}")
  cli::cli_li("Iterations: {.val {n_warmup}} warmup + {.val {n_iter}} sampling")
  cli::cli_li("Elapsed: {format(x$elapsed, digits = 3)}")
  cli::cli_end()

  ##### Diagnostics #####
  diag <- x$diagnose
  n_divergent <- sum(diag$num_divergent)
  n_max_treedepth <- sum(diag$num_max_treedepth)
  any_low_ebfmi <- any(diag$ebfmi < 0.3)

  cli::cli_h3("Diagnostics")
  if (n_divergent > 0) {
    cli::cli_alert_danger("{.val {n_divergent}} divergent transition{?s}")
  } else {
    cli::cli_alert_success("No divergent transitions")
  }
  if (n_max_treedepth > 0) {
    cli::cli_alert_danger(
      "{.val {n_max_treedepth}} transition{?s} hit max treedepth"
    )
  } else {
    cli::cli_alert_success("No max treedepth warnings")
  }
  if (any_low_ebfmi) {
    cli::cli_alert_danger("Low E-BFMI detected (< 0.3)")
  } else {
    cli::cli_alert_success("No E-BFMI warnings")
  }

  invisible(x)
}
