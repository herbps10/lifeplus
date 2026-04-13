functions {
  real error_rng(vector x, vector epsilon_params, int shock) {
    return normal_rng(0, epsilon_params[1]);
  }
}
data {
  int<lower=0, upper=1> fix_epsilon_sigma;
  real<lower=0> epsilon_sigma_fixed;
  real<lower=0> epsilon_sigma_prior_mu;
  real<lower=0> epsilon_sigma_prior_sd;
}
transformed data {
}
parameters {
  array[1 - fix_epsilon_sigma] real<lower=0> epsilon_sigma_raw;
}
transformed parameters{
  real epsilon_sigma;
  if(fix_epsilon_sigma == 1) {
    epsilon_sigma = epsilon_sigma_fixed;
  }
  else {
    epsilon_sigma = epsilon_sigma_raw[1];
  }

  if(tilted == 1) {
    ep_phi[1][, 1] = rep_vector(epsilon_sigma, C);
  }
}
model {
  epsilon_sigma ~ normal(epsilon_sigma_prior_mu, epsilon_sigma_prior_sd);
  if(shock_diff_mode == 1) {
    diff ~ normal(to_vector(transition_function + shock[, 2:T] - shock[, 1:(T - 1)]), epsilon_sigma);
  }
  else {
    diff ~ normal(to_vector(transition_function + shock), epsilon_sigma);
  }
}
generated quantities {
  vector[1] epsilon_params;
  epsilon_params[1] = epsilon_sigma;

  for(n in 1:size(log_lik)) {
    if(shock_diff_mode == 1) {
      log_lik[n] = normal_lpdf(diff[n] | to_vector(transition_function + shock[, 2:T] - shock[, 1:(T - 1)])[n], epsilon_sigma);
    }
    else {
      log_lik[n] = normal_lpdf(diff[n] | to_vector(transition_function + shock)[n], epsilon_sigma);
    }
  }
}
