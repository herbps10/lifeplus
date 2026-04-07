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
