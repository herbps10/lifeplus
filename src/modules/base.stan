functions {
}

data {
  int C;     // Number of countries
  int T;     // Number of time points
  int Tpred; // Total number of timepoints

  matrix[C, T] y;
  
  int num_grid;
  vector[num_grid] grid;

  int<lower=0, upper=1> include_prior;
  int<lower=0, upper=1> hierarchical;  
  int<lower=0, upper=1> shock_diff_mode; // Whether the shock term is delta_{c,t} - delta_{c,t-1} or simply delta_{c,t}
}
transformed data {
  vector[C * (T - 1)] diff = to_vector(y[, 2:T] - y[, 1:(T - 1)]);
  int shock_term = 0;
  int generate_shock_free = 0;

  int intermediate_grid_index = 0;
  for(i in 1:num_grid) {
    if(grid[i] == 90) {
      intermediate_grid_index = i;
      break;
    }
  }

  int T_shocks = T;
  if(shock_diff_mode == 0) {
    T_shocks = T - 1;
  }
}
parameters {
}
transformed parameters {
  matrix[C, T_shocks] shock = rep_matrix(0, C, T_shocks);
  matrix[C, T - 1] transition_function = rep_matrix(0, C, T - 1);
  array[include_prior] vector[C] first_transition;
  array[include_prior] vector[C] intermediate_transition;
  array[include_prior] vector[C] final_transition;

  if(include_prior == 1) {
    first_transition[1]        = rep_vector(0, C);
    intermediate_transition[1] = rep_vector(0, C);
    final_transition[1]        = rep_vector(0, C);
  }
}

model {
  // Prior on transition function prediction
  if(include_prior == 1) {
    to_vector(first_transition[1]) ~ normal(0, 25);
    to_vector(intermediate_transition[1]) ~ normal(0, 4);
    to_vector(final_transition[1]) ~ normal(1.15 / 10, 0.5);
  }
}
generated quantities {
  matrix[C, Tpred] eta;

  matrix[generate_shock_free * C, generate_shock_free * Tpred] eta_shockfree;
  matrix[shock_term * C, shock_term * (T_shocks + Tpred - T)] shock2;
  if(shock_term == 1) shock2 = rep_matrix(0, C, T_shocks + Tpred - T);

  eta[1:C, 1:T] = y;
  if(shock_term == 1) {
    shock2[1:C, 1:T_shocks] = shock;
  }

  if(generate_shock_free == 1) {
    if(shock_diff_mode == 0) {
      eta_shockfree[, 1] = rep_vector(0, C);
      eta_shockfree[, 2:T] = y[, 2:T] - shock;
    }
    else {
      eta_shockfree[, 1:T] = y - shock;
    }
  }

  matrix[C, num_grid] transition_function_pred;
  vector[num_grid * hierarchical] transition_function_pred_mean;
}
