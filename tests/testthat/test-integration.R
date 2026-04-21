test_that("lifeplus successfully runs with all model combinations", {
  skip_on_cran()
  skip_if_not_installed("cmdstanr")

  data <- lifeplus_simulate(num_areas = 2, num_times = 15, seed = 10016)

  transition_models <- list(
    dl_hier = transition_model_double_logistic(hierarchical = TRUE),
    dl_indep = transition_model_double_logistic(hierarchical = FALSE),
    gp_hier = transition_model_gaussian_process(hierarchical = TRUE),
    gp_indep = transition_model_gaussian_process(hierarchical = FALSE),
    sp_hier = transition_model_spline(hierarchical = TRUE),
    sp_indep = transition_model_spline(hierarchical = FALSE)
  )

  data_models <- list(
    normal = data_model_normal(),
    normal_fixed = data_model_normal(fixed_sd = 1),
    outlier = data_model_outlier(),
    outlier_fixed = data_model_outlier(fixed_sd = 1)
  )

  shock_models <- list(
    none = shock_model_none(),
    horseshoe = shock_model_regularized_horseshoe()
  )

  combos <- expand.grid(
    transition = names(transition_models),
    data = names(data_models),
    shock = names(shock_models),
    stringsAsFactors = FALSE
  )

  sampling_args <- list(
    data = data,
    y = "y",
    time = "time",
    area = "area",
    chains = 1,
    iter_warmup = 50,
    iter_sampling = 50,
    refresh = 0,
    show_messages = FALSE,
    show_exceptions = FALSE
  )

  for (index in seq_len(nrow(combos))) {
    combo <- combos[index, ]
    models <- list(
      transition_model = transition_models[[combo$transition]],
      data_model = data_models[[combo$data]],
      shock_model = shock_models[[combo$shock]]
    )

    if (
      !check_model_compatibility(
        models$transition_model,
        models$data_model,
        models$shock_model,
        error = FALSE
      )
    ) {
      next
    }

    label <- paste(combo$transition, combo$data, combo$shock, sep = " / ")

    args <- c(sampling_args, models)

    fit <- do.call(lifeplus, args)
    expect_equal(
      class(fit),
      "lifeplus",
      label = label,
      expected.label = "lifeplus"
    )
  }
})
