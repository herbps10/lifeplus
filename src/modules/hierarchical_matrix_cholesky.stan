data {
  array[num] int<lower=0, upper=1> var_constrain;
  int<lower=0, upper=1> var_multi;
  int<lower=0, upper=1> var_log_scale;
  vector[num] var_lower;
  vector[num] var_upper;
  vector[num] var_prior_mean;
  vector[num] var_prior_sd;

  vector[num] var_sigma_prior_mean;
  vector[num] var_sigma_prior_sd;
}
parameters {
  matrix[C, num] raw_var;
  array[hierarchical] vector[num] mu_var;
  array[hierarchical] vector[num] sigma_var_raw;
  array[hierarchical * var_multi] cholesky_factor_corr[num] L_Omega_var;
}
transformed parameters {
  matrix[C, num] var;
  array[hierarchical] vector[num] sigma_var;
  if(hierarchical == 1) {
    if(var_log_scale == 1) {
      sigma_var[1] = exp(sigma_var_raw[1]);
    }
    else {
      sigma_var[1] = sigma_var_raw[1];
    }
  }

  if(hierarchical && var_multi) {
    var = rep_matrix(mu_var[1]', C) + (diag_pre_multiply(sigma_var[1], L_Omega_var[1]) * raw_var')';
  }
  else if(hierarchical && !var_multi) {
    var = rep_matrix(mu_var[1]', C) + (diag_matrix(sigma_var[1]) * raw_var')';
  }
  else {
    var = raw_var;
  }

  if(tilted == 1 && hierarchical == 1) {
    ep_phi[1][1:num] = mu_var[1];
    ep_phi[1][(num + 1):(2 * num)] = sigma_var_raw[1];
  }

  for(n in 1:num) {
    if(var_constrain[n] == 1) {
      var[, n] = inv_logit(var[, n]) * (var_upper[n] - var_lower[n]) + var_lower[n];
    }
  }
}
model {
  if(hierarchical == 1) {
    to_vector(raw_var) ~ std_normal();

    mu_var[1] ~ normal(var_prior_mean, var_prior_sd);
    sigma_var_raw[1] ~ normal(var_sigma_prior_mean, var_sigma_prior_sd);

    if(var_multi) {
      L_Omega_var[1] ~ lkj_corr_cholesky(1.0);
    }
  }
  else {
    for(d in 1:D) {
      to_vector(raw_var[, d]) ~ normal(var_prior_mean[d], var_prior_sd[d]);
    }
  }
}
generated quantities {
  corr_matrix[num * hierarchical * var_multi] Omega_var;
  if(hierarchical && var_multi == 1) Omega_var = multiply_lower_tri_self_transpose(L_Omega_var[1]);
}
