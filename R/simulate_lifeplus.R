#' Generate simulated data for testing lifeplus models.
#'
#' @param num_areas number of areas to simulate
#' @param num_times number of time points to simulate
#' @param seed optional random number seed
#' @return data frame
#' @importFrom stats runif rnorm
#' @export
simulate_lifeplus <- function(num_areas, num_times, seed = NULL) {
  checkmate::check_integer(num_areas, lower = 1)
  checkmate::check_integer(num_times, lower = 1)
  checkmate::check_integer(seed, null.ok = TRUE)

  if (!is.null(seed)) {
    set.seed(seed)
  }

  areas <- 1:num_areas
  times <- 1:num_times

  diffs <- matrix(
    stats::rnorm(num_areas * (num_times - 1), 1, 2),
    ncol = num_times - 1,
    nrow = num_areas
  )

  y <- t(apply(cbind(stats::runif(num_areas, 50, 70), diffs), 1, cumsum))

  data.frame(
    area = rep(areas, each = num_times),
    time = rep(times, times = num_areas),
    y = as.vector(t(y))
  )
}
