load_model <- function(transition_model, data_model, shock_model) {
  model_name <- paste0(paste0(
    c(transition_model$name, data_model$name, shock_model$name),
    collapse = "_"
  ))

  return(instantiate::stan_package_model(
    model_name,
    package = "lifeplus"
  ))
}

check_model_compatibility <- function(
  transition_model,
  data_model,
  shock_model,
  error = TRUE
) {
  if (data_model$name == "outlier" & shock_model$name != "none") {
    if (error) {
      cli::cli_abort(
        "data_model_outlier can only be used with shock_model_none"
      )
    }
    return(FALSE)
  }
  return(TRUE)
}
