process_temporal <- function(fit, parallel_chains, probs) {
  temporal_variables <- c("eta")
  if (fit$shock_model$name != "none") {
    temporal_variables <- c(temporal_variables, "shock2", "eta_shockfree")
  }

  post <- fit$samples$summary(
    temporal_variables,
    ~ stats::quantile(
      .x,
      probs = probs
    ),
    .cores = parallel_chains
  )

  x <- stringr::str_match(
    post$variable,
    "([a-zA-Z0-9]+)\\[([0-9]+),([0-9]+)\\]"
  )
  post$variable <- x[, 2]
  c_index <- as.numeric(x[, 3])
  t_index <- as.numeric(x[, 4])

  post$area <- fit$areas[c_index]
  post$time <- fit$times[t_index]

  return(post)
}

process_transition_functions <- function(fit, parallel_chains, probs) {
  post <- fit$samples$summary(
    "transition_function_pred",
    ~ stats::quantile(.x, probs = probs),
    posterior::default_convergence_measures(),
    .cores = parallel_chains
  )

  x <- stringr::str_match(
    post$variable,
    "([a-zA-Z0-9]+)\\[([0-9]+),([0-9]+)\\]"
  )

  post$variable <- x[, 2]
  c_index <- as.numeric(x[, 3])
  grid_index <- as.numeric(x[, 4])

  post$area <- fit$areas[c_index]
  post$x <- fit$grid[grid_index]

  return(post)
}


process_transition_function_mean <- function(fit, parallel_chains, probs) {
  hierarchical <- !is.null(fit$stan_data$hierarchical) &&
    as.logical(fit$stan_data$hierarchical) == TRUE

  if (hierarchical == FALSE) {
    return(NULL)
  }

  post <- fit$samples$summary(
    "transition_function_pred_mean",
    ~ stats::quantile(.x, probs = probs),
    posterior::default_convergence_measures(),
    .cores = parallel_chains
  )

  x <- stringr::str_match(
    post$variable,
    "([a-zA-Z0-9]+)\\[([0-9]+)\\]"
  )

  post$variable <- x[, 2]
  grid_index <- as.numeric(x[, 3])

  post$x <- fit$grid[grid_index]

  return(post)
}

process_transition_params <- function(fit, parallel_chains, probs) {
  if (fit$transition_model$name == "logistic") {
    # Extract Delta parameters
    post <- fit$samples$summary(
      c("Delta"),
      ~ stats::quantile(
        .x,
        probs = probs
      ),
      posterior::default_convergence_measures(),
      .cores = parallel_chains
    )

    x <- stringr::str_match(
      post$variable,
      "([a-zA-Z0-9]+)\\[([0-9]+),([0-9]+)\\]"
    )

    vname <- function(x) c("Delta1", "Delta2", "Delta3", "Delta4", "k", "z")[x]

    c_index <- as.numeric(x[, 3])
    k_index <- as.numeric(x[, 4])
    post$variable <- vname(k_index)

    post$area <- fit$areas[c_index]

    hierarchical <- !is.null(fit$stan_data$hierarchical) &&
      as.logical(fit$stan_data$hierarchical) == TRUE

    correlation <- NULL
    if (hierarchical) {
      correlation <- fit$samples$summary("Omega_Delta")
      x <- stringr::str_match(
        correlation$variable,
        "([a-zA-Z0-9]+)\\[([0-9]+),([0-9]+),([0-9]+)\\]"
      )
      correlation$variable1 <- vname(as.numeric(x[, 3]))
      correlation$variable2 <- vname(as.numeric(x[, 4]))
    }

    return(list(
      params = post,
      correlation = correlation
    ))
  }
  return(NULL)
}

lifeplus_posteriors <- function(fit, parallel_chains = NULL) {
  if (is.null(parallel_chains)) {
    parallel_chains <- 1
  }

  probs <- c(
    0.001,
    0.01,
    0.025,
    0.1,
    0.25,
    0.5,
    0.75,
    0.9,
    0.975,
    0.99,
    0.999
  )

  temporal <- process_temporal(fit, parallel_chains, probs)
  transition_functions <- process_transition_functions(
    fit,
    parallel_chains,
    probs
  )

  transition_function_mean <- process_transition_function_mean(
    fit,
    parallel_chains,
    probs
  )

  transition_params <- process_transition_params(
    fit,
    parallel_chains,
    probs
  )

  ans <- list(
    temporal = temporal,
    transition_functions = transition_functions,
    transition_function_mean = transition_function_mean,
    transition_params = transition_params
  )

  ans
}
