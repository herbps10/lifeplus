test_that("lifeplus successfully runs with all combinations of models", {
  data <- simulate_lifeplus(2, 15, 10016)

  transition_models <- list(
    transition_model_double_logistic(hierarchical = TRUE),
    transition_model_double_logistic(hierarchical = FALSE),
    transition_model_gaussian_process(hierarchical = TRUE),
    transition_model_gaussian_process(hierarchical = FALSE),
    transition_model_spline(hierarchical = TRUE),
    transition_model_spline(hierarchical = FALSE)
  )

  data_models <- list(
    data_model_normal(),
    data_model_outlier()
  )

  shock_models <- list(
    shock_model_none(),
    shock_model_regularized_horseshoe()
  )

  # Test it works with different hierarchical options
  for (transition_model_index in seq_along(transition_models)) {
    for (data_model_index in seq_along(data_models)) {
      for (shock_model_index in seq_along(shock_models)) {
        if (
          shock_models[[shock_model_index]]$name == "regularized_horseshoe" &
            data_models[[data_model_index]]$name == "outlier"
        ) {
          next
        }

        expect_no_error(
          fit <- lifeplus(
            data,
            "y",
            "time",
            area = "area",
            transition_model = transition_models[[transition_model_index]],
            data_model = data_models[[data_model_index]],
            shock_model = shock_models[[shock_model_index]],
            chains = 1,
            iter_warmup = 300,
            iter_sampling = 300,
            refresh = 0,
            show_messages = FALSE,
            show_exceptions = FALSE
          )
        )
      }

      expect_s3_class(fit, "lifeplus")
      expect_equal(sum(fit$posteriors$temporal$variable == "eta"), 30)
    }
  }
})
