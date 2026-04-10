functions {
  real error_rng(vector x, vector epsilon_params, int shock) {
    return normal_rng(0, epsilon_params[1]);
  }
}
data {
  real<lower=0> outlier_threshold;

  int<lower=0, upper=1> fix_epsilon_sigma;
  real<lower=0> epsilon_sigma_fixed;
  real<lower=0> epsilon_sigma_prior_mu;
  real<lower=0> epsilon_sigma_prior_sd;
}
transformed data {
  int n_below_threshold = 0;
  for(i in 1:(C * (T - 1))) {
    n_below_threshold += (abs(diff[i]) < outlier_threshold) ? 1 : 0;
  }
  array[n_below_threshold] int indices_below_threshold;
  
  {
    int index = 1;
    for(i in 1:(C * (T - 1))) {
      if(abs(diff[i]) < outlier_threshold) {
        indices_below_threshold[index] = i;
        index += 1;
      }
    }
  }
}
parameters {
  array[1 - fix_epsilon_sigma] real<lower=0> epsilon_sigma_raw;
}
transformed parameters {
  real epsilon_sigma;
  if(fix_epsilon_sigma == 1) {
    epsilon_sigma = epsilon_sigma_fixed;
  }
  else {
    epsilon_sigma = epsilon_sigma_raw[1];
  }
}
model {
  epsilon_sigma ~ normal(epsilon_sigma_prior_mu, epsilon_sigma_prior_sd);

  if(outlier_threshold < 1000) {
    diff[indices_below_threshold] ~ normal(to_vector(transition_function)[indices_below_threshold], epsilon_sigma);
  }
  else {
    diff ~ normal(to_vector(transition_function), epsilon_sigma);
  }
}
generated quantities {
  vector[1] epsilon_params;
  epsilon_params[1] = epsilon_sigma;

  if(outlier_threshold < 1000) {
    log_lik = rep_vector(0, size(log_lik));
    for(n in 1:size(indices_below_threshold)) {
      log_lik[indices_below_threshold[n]] = normal_lpdf(diff[indices_below_threshold[n]] | to_vector(transition_function)[indices_below_threshold[n]], epsilon_sigma);
    }
  }
  else {
    for(n in 1:size(log_lik)) {
      log_lik[n] = normal_lpdf(diff[n] | to_vector(transition_function)[n], epsilon_sigma);
    }
  }
}
