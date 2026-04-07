#' Plot posterior projections
#'
#' @param fit lifeplus fit object
#' @param areas vector of areas to plot (default: all areas in fit)
#' @param start_time first time point to plot (defaults to first data point)
#' @param endart_time final time point to plot (defaults to final projection time point)
#' @param data logical indicating whether to plot raw data (default: TRUE)
#'
#' @export
plot_projections <- function(
  fit,
  areas = fit$areas,
  start_time = min(fit$times),
  end_time = max(fit$times),
  data = TRUE
) {
  checkmate::check_class(fit, "lifeplus")
  checkmate::check_subset(areas, fit$areas)
  checkmate::check_flag(data)

  plot_data <- fit$posteriors$temporal[
    fit$posteriors$temporal$variable == "eta" &
      fit$posteriors$temporal$area %in%
        areas &
      fit$posteriors$temporal$time >= start_time &
      fit$posteriors$temporal$time <= end_time,
  ]

  alpha <- 0.8
  p <- ggplot2::ggplot(plot_data, ggplot2::aes(x = time, y = `50%`)) +
    ggplot2::geom_ribbon(
      ggplot2::aes(
        ymin = `2.5%`,
        ymax = `97.5%`,
        fill = "95%"
      ),
      alpha = alpha
    ) +
    ggplot2::geom_ribbon(
      ggplot2::aes(
        ymin = `10%`,
        ymax = `90%`,
        fill = "80%"
      ),
      alpha = alpha
    ) +
    ggplot2::geom_ribbon(
      ggplot2::aes(
        ymin = `25%`,
        ymax = `75%`,
        fill = "50%"
      ),
      alpha = alpha
    ) +
    ggplot2::geom_line() +
    ggplot2::scale_fill_brewer(direction = -1) +
    ggplot2::facet_wrap(~area)

  if (data == TRUE) {
    point_data <- data.frame(
      x = fit$data[[fit$time]],
      y = fit$data[[fit$y]],
      area = fit$data[[fit$area]]
    )
    p <- p + ggplot2::geom_point(data = point_data, ggplot2::aes(x, y))
  }

  p
}

#' Plot posterior shocks
#'
#' @param fit lifeplus fit object
#' @param areas vector of areas to plot (default: all areas in fit)
#' @param start_time first time point to plot (defaults to first data point)
#' @param endart_time final time point to plot (defaults to final projection time point)
#'
#' @export
plot_shocks <- function(
  fit,
  areas = fit$areas,
  start_time = min(fit$times),
  end_time = max(fit$times)
) {
  checkmate::check_class(fit, "lifeplus")
  checkmate::check_subset(areas, fit$areas)

  plot_data <- fit$posteriors$temporal[
    fit$posteriors$temporal$variable == "shock2" &
      fit$posteriors$temporal$area %in%
        areas &
      fit$posteriors$temporal$time >= start_time &
      fit$posteriors$temporal$time <= end_time,
  ]

  alpha <- 0.8
  p <- ggplot2::ggplot(plot_data, ggplot2::aes(x = time, y = `50%`)) +
    ggplot2::geom_ribbon(
      ggplot2::aes(
        ymin = `2.5%`,
        ymax = `97.5%`,
        fill = "95%"
      ),
      alpha = alpha
    ) +
    ggplot2::geom_ribbon(
      ggplot2::aes(
        ymin = `10%`,
        ymax = `90%`,
        fill = "80%"
      ),
      alpha = alpha
    ) +
    ggplot2::geom_ribbon(
      ggplot2::aes(
        ymin = `25%`,
        ymax = `75%`,
        fill = "50%"
      ),
      alpha = alpha
    ) +
    ggplot2::geom_line(size = 0.5) +
    ggplot2::scale_fill_brewer(direction = -1) +
    ggplot2::facet_wrap(~area)

  p
}


#' Plot posterior transition functions
#'
#' @param fit lifeplus fit object
#' @param areas vector of areas to plot (default: all areas in fit)
#' @param data logical indicating whether to plot raw data (default: TRUE)
#'
#' @return ggplot2 object
plot_transitions <- function(fit, areas = fit$areas, data = TRUE) {
  checkmate::check_class(fit, "lifeplus")
  checkmate::check_subset(areas, fit$areas)
  checkmate::check_flag(data)

  plot_data <- fit$posteriors$transition_functions[
    fit$posteriors$transition_functions$area %in% areas,
  ]

  alpha <- 0.8
  p <- ggplot2::ggplot(
    plot_data,
    ggplot2::aes(x = x, y = `50%`)
  ) +
    ggplot2::geom_ribbon(
      ggplot2::aes(
        ymin = `2.5%`,
        ymax = `97.5%`,
        fill = "95%"
      ),
      alpha = alpha
    ) +
    ggplot2::geom_ribbon(
      ggplot2::aes(
        ymin = `10%`,
        ymax = `90%`,
        fill = "80%"
      ),
      alpha = alpha
    ) +
    ggplot2::geom_ribbon(
      ggplot2::aes(
        ymin = `25%`,
        ymax = `75%`,
        fill = "50%"
      ),
      alpha = alpha
    ) +
    ggplot2::geom_line() +
    ggplot2::scale_fill_brewer(direction = -1) +
    ggplot2::facet_wrap(~area)

  if (data == TRUE) {
    diff <- t(apply(fit$stan_data$y, 1, diff))
    n <- ncol(fit$stan_data$y) - 1
    point_data <- data.frame(
      area = rep(fit$areas[match(areas, fit$areas)], times = n),
      px = as.vector(fit$stan_data$y[match(areas, fit$areas), 1:n]),
      py = as.vector(diff[match(areas, fit$areas), ])
    )

    p <- p +
      ggplot2::geom_point(ggplot2::aes(x = px, y = py), data = point_data)
  }

  return(p)
}
