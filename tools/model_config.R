generate_config <- function(transition_model, data_model, shocks) {
  if (
    !(transition_model %in% c("double_logistic", "spline", "gaussian_process"))
  ) {
    stop(
      "transition_model must be one of: double_logistic, spline, gaussian_process"
    )
  }
  if (!(data_model %in% c("normal", "outlier", "mixture"))) {
    stop("data_model must be one of: normal, outlier, mixture")
  }
  if (!(shocks %in% c(FALSE, TRUE))) {
    stop("shocks must be one of: FALSE, TRUE")
  }

  inc <- function(x) here::here("src/modules", x)

  path <- glue::glue(
    "../src/stan/{transition_model}_{data_model}_{ifelse(shocks, 'regularized_horseshoe', 'none')}.stan"
  )

  deps <- list(
    output = path
  )
  deps[[inc("base.stan")]] <- "empty"

  if (
    transition_model == "gaussian_process" || data_model == "heteroskedastic"
  ) {
    deps[[inc("approximate_gp.stan")]] <- "empty"
  }

  if (data_model == "normal") {
    data_model_path <- inc("data_model_normal.stan")
  } else if (data_model == "outlier") {
    data_model_path <- inc("data_model_outlier.stan")
  } else if (data_model == "mixture") {
    data_model_path <- inc("data_model_mixture.stan")
  } else if (data_model == "heteroskedastic") {
    data_model_path <- inc("data_model_heteroskedastic.stan")
  }
  deps[[data_model_path]] <- "empty"

  if (shocks == TRUE) {
    deps[[inc("shock.stan")]] <- "empty"
  }

  if (transition_model == "double_logistic") {
    deps[[inc("Delta.stan")]] <- "empty"
    deps[[inc("hierarchical_matrix_cholesky.stan")]] <- list(
      "var" = "Delta",
      num = "D"
    )
    deps[[inc("transition_double_logistic.stan")]] <- "empty"
  } else if (transition_model == "gaussian_process") {
    deps[[inc("hierarchical_matrix.stan")]] <- list(
      "var" = "beta",
      "num" = "M"
    )
    deps[[inc("transition_gaussian_process.stan")]] <- "empty"
  } else if (transition_model == "spline") {
    deps[[inc("transition_spline_data.stan")]] <- "empty"
    deps[[inc("hierarchical_matrix.stan")]] <- list(
      "var" = "alpha",
      "num" = "num_basis"
    )
    deps[[inc("transition_spline.stan")]] <- "empty"
  }

  deps
}

configs <- tidyr::expand_grid(
  transition_model = c("double_logistic", "spline", "gaussian_process"),
  data_model = c("normal", "outlier"),
  shock = c(FALSE, TRUE)
) |>
  dplyr::filter(!(shock == TRUE & data_model == "outlier")) |>
  dplyr::mutate(
    config = purrr::pmap(
      list(transition_model, data_model, shock),
      generate_config
    )
  )

jsonlite::toJSON(configs$config, auto_unbox = TRUE) |>
  stringr::str_replace_all("\\.[0-9]", "") |>
  stringr::str_replace_all("\"empty\"", "{}") |>
  readr::write_file(here::here("tools/config.json"))
