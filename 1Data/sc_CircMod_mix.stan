// generated with brms 2.22.0
functions {
  /* softplus link function inverse to 'log1p_exp'
   * Args:
   *   x: a positive scalar
   * Returns:
   *   a scalar in (-Inf, Inf)
   */
   real log_expm1(real x) {
     return log(expm1(x));
   }
  /* softplus link function inverse to 'log1p_exp' (vectorized)
   * Args:
   *   x: a positive vector
   * Returns:
   *   a vector in (-Inf, Inf)
   */
   vector log_expm1(vector x) {
     return log(expm1(x));
   }

 /* compute correlated group-level effects
  * Args:
  *   z: matrix of unscaled group-level effects
  *   SD: vector of standard deviation parameters
  *   L: cholesky factor correlation matrix
  * Returns:
  *   matrix of scaled group-level effects
  */
  matrix scale_r_cor(matrix z, vector SD, matrix L) {
    // r is stored in another dimension order than z
    return transpose(diag_pre_multiply(SD, L) * z);
  }




  
  real unwrap_von_mises_lpdf(real y, real mu, real kappa) {
    return von_mises_lpdf(y | mod_circular(mu), kappa);
  }
  real unwrap_von_mises_rng(real mu, real kappa) {
    return von_mises_rng( mod_circular(mu) , kappa);
  }

  
  real inv_mod_circular(real y) {
    return mod_circular(y);
  }

  
  real mod_circular(real y) {
    return fmod(y + pi(), 2*pi()) - pi();
  }

}
data {
  int<lower=1> N;  // total number of observations
  vector[N] Y;  // response variable
  int<lower=1> K_mu1;  // number of population-level effects
  matrix[N, K_mu1] X_mu1;  // population-level design matrix
  int<lower=1> Kc_mu1;  // number of population-level effects after centering
  int<lower=1> K_kappa1;  // number of population-level effects
  matrix[N, K_kappa1] X_kappa1;  // population-level design matrix
  int<lower=1> Kc_kappa1;  // number of population-level effects after centering
  int<lower=1> K_mu2;  // number of population-level effects
  matrix[N, K_mu2] X_mu2;  // population-level design matrix
  int<lower=1> Kc_mu2;  // number of population-level effects after centering
  int<lower=1> K_kappa2;  // number of population-level effects
  matrix[N, K_kappa2] X_kappa2;  // population-level design matrix
  int<lower=1> Kc_kappa2;  // number of population-level effects after centering
  int<lower=1> K_theta1;  // number of population-level effects
  matrix[N, K_theta1] X_theta1;  // population-level design matrix
  int<lower=1> Kc_theta1;  // number of population-level effects after centering
  // data for group-level effects of ID 1
  int<lower=1> N_1;  // number of grouping levels
  int<lower=1> M_1;  // number of coefficients per level
  array[N] int<lower=1> J_1;  // grouping indicator per observation
  // group-level predictor values
  vector[N] Z_1_mu1_1;
  vector[N] Z_1_mu1_2;
  vector[N] Z_1_mu1_3;
  vector[N] Z_1_mu1_4;
  int<lower=1> NC_1;  // number of group-level correlations
  // data for group-level effects of ID 2
  int<lower=1> N_2;  // number of grouping levels
  int<lower=1> M_2;  // number of coefficients per level
  array[N] int<lower=1> J_2;  // grouping indicator per observation
  // group-level predictor values
  vector[N] Z_2_mu2_1;
  vector[N] Z_2_mu2_2;
  vector[N] Z_2_mu2_3;
  vector[N] Z_2_mu2_4;
  int<lower=1> NC_2;  // number of group-level correlations
  // data for group-level effects of ID 3
  int<lower=1> N_3;  // number of grouping levels
  int<lower=1> M_3;  // number of coefficients per level
  array[N] int<lower=1> J_3;  // grouping indicator per observation
  // group-level predictor values
  vector[N] Z_3_kappa1_1;
  vector[N] Z_3_kappa1_2;
  vector[N] Z_3_kappa1_3;
  vector[N] Z_3_kappa1_4;
  int<lower=1> NC_3;  // number of group-level correlations
  // data for group-level effects of ID 4
  int<lower=1> N_4;  // number of grouping levels
  int<lower=1> M_4;  // number of coefficients per level
  array[N] int<lower=1> J_4;  // grouping indicator per observation
  // group-level predictor values
  vector[N] Z_4_kappa2_1;
  vector[N] Z_4_kappa2_2;
  vector[N] Z_4_kappa2_3;
  vector[N] Z_4_kappa2_4;
  int<lower=1> NC_4;  // number of group-level correlations
  // data for group-level effects of ID 5
  int<lower=1> N_5;  // number of grouping levels
  int<lower=1> M_5;  // number of coefficients per level
  array[N] int<lower=1> J_5;  // grouping indicator per observation
  // group-level predictor values
  vector[N] Z_5_theta1_1;
  vector[N] Z_5_theta1_2;
  vector[N] Z_5_theta1_3;
  vector[N] Z_5_theta1_4;
  int<lower=1> NC_5;  // number of group-level correlations
  int prior_only;  // should the likelihood be ignored?
}
transformed data {
  matrix[N, Kc_mu1] Xc_mu1;  // centered version of X_mu1 without an intercept
  vector[Kc_mu1] means_X_mu1;  // column means of X_mu1 before centering
  matrix[N, Kc_kappa1] Xc_kappa1;  // centered version of X_kappa1 without an intercept
  vector[Kc_kappa1] means_X_kappa1;  // column means of X_kappa1 before centering
  matrix[N, Kc_mu2] Xc_mu2;  // centered version of X_mu2 without an intercept
  vector[Kc_mu2] means_X_mu2;  // column means of X_mu2 before centering
  matrix[N, Kc_kappa2] Xc_kappa2;  // centered version of X_kappa2 without an intercept
  vector[Kc_kappa2] means_X_kappa2;  // column means of X_kappa2 before centering
  matrix[N, Kc_theta1] Xc_theta1;  // centered version of X_theta1 without an intercept
  vector[Kc_theta1] means_X_theta1;  // column means of X_theta1 before centering
  for (i in 2:K_mu1) {
    means_X_mu1[i - 1] = mean(X_mu1[, i]);
    Xc_mu1[, i - 1] = X_mu1[, i] - means_X_mu1[i - 1];
  }
  for (i in 2:K_kappa1) {
    means_X_kappa1[i - 1] = mean(X_kappa1[, i]);
    Xc_kappa1[, i - 1] = X_kappa1[, i] - means_X_kappa1[i - 1];
  }
  for (i in 2:K_mu2) {
    means_X_mu2[i - 1] = mean(X_mu2[, i]);
    Xc_mu2[, i - 1] = X_mu2[, i] - means_X_mu2[i - 1];
  }
  for (i in 2:K_kappa2) {
    means_X_kappa2[i - 1] = mean(X_kappa2[, i]);
    Xc_kappa2[, i - 1] = X_kappa2[, i] - means_X_kappa2[i - 1];
  }
  for (i in 2:K_theta1) {
    means_X_theta1[i - 1] = mean(X_theta1[, i]);
    Xc_theta1[, i - 1] = X_theta1[, i] - means_X_theta1[i - 1];
  }
}
parameters {
  vector[Kc_mu1] b_mu1;  // regression coefficients
  vector[Kc_kappa1] b_kappa1;  // regression coefficients
  real Intercept_kappa1;  // temporary intercept for centered predictors
  vector[Kc_mu2] b_mu2;  // regression coefficients
  vector[Kc_kappa2] b_kappa2;  // regression coefficients
  real Intercept_kappa2;  // temporary intercept for centered predictors
  vector[Kc_theta1] b_theta1;  // regression coefficients
  real Intercept_theta1;  // temporary intercept for centered predictors
  ordered[2] ordered_Intercept;  // to identify mixtures
  vector<lower=0>[M_1] sd_1;  // group-level standard deviations
  matrix[M_1, N_1] z_1;  // standardized group-level effects
  cholesky_factor_corr[M_1] L_1;  // cholesky factor of correlation matrix
  vector<lower=0>[M_2] sd_2;  // group-level standard deviations
  matrix[M_2, N_2] z_2;  // standardized group-level effects
  cholesky_factor_corr[M_2] L_2;  // cholesky factor of correlation matrix
  vector<lower=0>[M_3] sd_3;  // group-level standard deviations
  matrix[M_3, N_3] z_3;  // standardized group-level effects
  cholesky_factor_corr[M_3] L_3;  // cholesky factor of correlation matrix
  vector<lower=0>[M_4] sd_4;  // group-level standard deviations
  matrix[M_4, N_4] z_4;  // standardized group-level effects
  cholesky_factor_corr[M_4] L_4;  // cholesky factor of correlation matrix
  vector<lower=0>[M_5] sd_5;  // group-level standard deviations
  matrix[M_5, N_5] z_5;  // standardized group-level effects
  cholesky_factor_corr[M_5] L_5;  // cholesky factor of correlation matrix
}
transformed parameters {
  // identify mixtures via ordering of the intercepts
  real Intercept_mu1 = ordered_Intercept[1];
  // identify mixtures via ordering of the intercepts
  real Intercept_mu2 = ordered_Intercept[2];
  matrix[N_1, M_1] r_1;  // actual group-level effects
  // using vectors speeds up indexing in loops
  vector[N_1] r_1_mu1_1;
  vector[N_1] r_1_mu1_2;
  vector[N_1] r_1_mu1_3;
  vector[N_1] r_1_mu1_4;
  matrix[N_2, M_2] r_2;  // actual group-level effects
  // using vectors speeds up indexing in loops
  vector[N_2] r_2_mu2_1;
  vector[N_2] r_2_mu2_2;
  vector[N_2] r_2_mu2_3;
  vector[N_2] r_2_mu2_4;
  matrix[N_3, M_3] r_3;  // actual group-level effects
  // using vectors speeds up indexing in loops
  vector[N_3] r_3_kappa1_1;
  vector[N_3] r_3_kappa1_2;
  vector[N_3] r_3_kappa1_3;
  vector[N_3] r_3_kappa1_4;
  matrix[N_4, M_4] r_4;  // actual group-level effects
  // using vectors speeds up indexing in loops
  vector[N_4] r_4_kappa2_1;
  vector[N_4] r_4_kappa2_2;
  vector[N_4] r_4_kappa2_3;
  vector[N_4] r_4_kappa2_4;
  matrix[N_5, M_5] r_5;  // actual group-level effects
  // using vectors speeds up indexing in loops
  vector[N_5] r_5_theta1_1;
  vector[N_5] r_5_theta1_2;
  vector[N_5] r_5_theta1_3;
  vector[N_5] r_5_theta1_4;
  real lprior = 0;  // prior contributions to the log posterior
  // compute actual group-level effects
  r_1 = scale_r_cor(z_1, sd_1, L_1);
  r_1_mu1_1 = r_1[, 1];
  r_1_mu1_2 = r_1[, 2];
  r_1_mu1_3 = r_1[, 3];
  r_1_mu1_4 = r_1[, 4];
  // compute actual group-level effects
  r_2 = scale_r_cor(z_2, sd_2, L_2);
  r_2_mu2_1 = r_2[, 1];
  r_2_mu2_2 = r_2[, 2];
  r_2_mu2_3 = r_2[, 3];
  r_2_mu2_4 = r_2[, 4];
  // compute actual group-level effects
  r_3 = scale_r_cor(z_3, sd_3, L_3);
  r_3_kappa1_1 = r_3[, 1];
  r_3_kappa1_2 = r_3[, 2];
  r_3_kappa1_3 = r_3[, 3];
  r_3_kappa1_4 = r_3[, 4];
  // compute actual group-level effects
  r_4 = scale_r_cor(z_4, sd_4, L_4);
  r_4_kappa2_1 = r_4[, 1];
  r_4_kappa2_2 = r_4[, 2];
  r_4_kappa2_3 = r_4[, 3];
  r_4_kappa2_4 = r_4[, 4];
  // compute actual group-level effects
  r_5 = scale_r_cor(z_5, sd_5, L_5);
  r_5_theta1_1 = r_5[, 1];
  r_5_theta1_2 = r_5[, 2];
  r_5_theta1_3 = r_5[, 3];
  r_5_theta1_4 = r_5[, 4];
  lprior += normal_lpdf(b_mu1 | pi()/3, pi()/4);
  lprior += normal_lpdf(Intercept_mu1 | 0, pi()/12);
  lprior += normal_lpdf(b_kappa1 | 0, 1);
  lprior += normal_lpdf(Intercept_kappa1 | 5, 1);
  lprior += normal_lpdf(b_mu2 | -pi()/3, pi()/4);
  lprior += normal_lpdf(Intercept_mu2 | 0, pi()/12);
  lprior += normal_lpdf(b_kappa2 | 0, 1);
  lprior += normal_lpdf(Intercept_kappa2 | 5, 1);
  lprior += normal_lpdf(b_theta1 | 0, 10);
  lprior += normal_lpdf(Intercept_theta1 | 10, 1);
  lprior += lognormal_lpdf(sd_1[1] | log(pi()/15), 0.5);
  lprior += lognormal_lpdf(sd_1[2] | log(pi()/9), 0.3);
  lprior += lognormal_lpdf(sd_1[3] | log(pi()/9), 0.3);
  lprior += lognormal_lpdf(sd_1[4] | log(pi()/9), 0.3);
  lprior += lkj_corr_cholesky_lpdf(L_1 | 1);
  lprior += lognormal_lpdf(sd_2[1] | log(pi()/15), 0.5);
  lprior += lognormal_lpdf(sd_2[2] | log(pi()/9), 0.3);
  lprior += lognormal_lpdf(sd_2[3] | log(pi()/9), 0.3);
  lprior += lognormal_lpdf(sd_2[4] | log(pi()/9), 0.3);
  lprior += lkj_corr_cholesky_lpdf(L_2 | 1);
  lprior += student_t_lpdf(sd_3 | 3, 0, 1)
    - 4 * student_t_lccdf(0 | 3, 0, 1);
  lprior += lkj_corr_cholesky_lpdf(L_3 | 1);
  lprior += student_t_lpdf(sd_4 | 3, 0, 1)
    - 4 * student_t_lccdf(0 | 3, 0, 1);
  lprior += lkj_corr_cholesky_lpdf(L_4 | 1);
  lprior += student_t_lpdf(sd_5 | 3, 0, 0.5)
    - 4 * student_t_lccdf(0 | 3, 0, 0.5);
  lprior += lkj_corr_cholesky_lpdf(L_5 | 1);
}
model {
  // likelihood including constants
  if (!prior_only) {
    // initialize linear predictor term
    vector[N] mu1 = rep_vector(0.0, N);
    // initialize linear predictor term
    vector[N] kappa1 = rep_vector(0.0, N);
    // initialize linear predictor term
    vector[N] mu2 = rep_vector(0.0, N);
    // initialize linear predictor term
    vector[N] kappa2 = rep_vector(0.0, N);
    // initialize linear predictor term
    vector[N] theta1 = rep_vector(0.0, N);
    vector[N] theta2 = rep_vector(0.0, N);
    real log_sum_exp_theta;
    mu1 += Intercept_mu1 + Xc_mu1 * b_mu1;
    kappa1 += Intercept_kappa1 + Xc_kappa1 * b_kappa1;
    mu2 += Intercept_mu2 + Xc_mu2 * b_mu2;
    kappa2 += Intercept_kappa2 + Xc_kappa2 * b_kappa2;
    theta1 += Intercept_theta1 + Xc_theta1 * b_theta1;
    for (n in 1:N) {
      // add more terms to the linear predictor
      mu1[n] += r_1_mu1_1[J_1[n]] * Z_1_mu1_1[n] + r_1_mu1_2[J_1[n]] * Z_1_mu1_2[n] + r_1_mu1_3[J_1[n]] * Z_1_mu1_3[n] + r_1_mu1_4[J_1[n]] * Z_1_mu1_4[n];
    }
    for (n in 1:N) {
      // add more terms to the linear predictor
      kappa1[n] += r_3_kappa1_1[J_3[n]] * Z_3_kappa1_1[n] + r_3_kappa1_2[J_3[n]] * Z_3_kappa1_2[n] + r_3_kappa1_3[J_3[n]] * Z_3_kappa1_3[n] + r_3_kappa1_4[J_3[n]] * Z_3_kappa1_4[n];
    }
    for (n in 1:N) {
      // add more terms to the linear predictor
      mu2[n] += r_2_mu2_1[J_2[n]] * Z_2_mu2_1[n] + r_2_mu2_2[J_2[n]] * Z_2_mu2_2[n] + r_2_mu2_3[J_2[n]] * Z_2_mu2_3[n] + r_2_mu2_4[J_2[n]] * Z_2_mu2_4[n];
    }
    for (n in 1:N) {
      // add more terms to the linear predictor
      kappa2[n] += r_4_kappa2_1[J_4[n]] * Z_4_kappa2_1[n] + r_4_kappa2_2[J_4[n]] * Z_4_kappa2_2[n] + r_4_kappa2_3[J_4[n]] * Z_4_kappa2_3[n] + r_4_kappa2_4[J_4[n]] * Z_4_kappa2_4[n];
    }
    for (n in 1:N) {
      // add more terms to the linear predictor
      theta1[n] += r_5_theta1_1[J_5[n]] * Z_5_theta1_1[n] + r_5_theta1_2[J_5[n]] * Z_5_theta1_2[n] + r_5_theta1_3[J_5[n]] * Z_5_theta1_3[n] + r_5_theta1_4[J_5[n]] * Z_5_theta1_4[n];
    }
    kappa1 = log1p_exp(kappa1);
    kappa2 = log1p_exp(kappa2);
    for (n in 1:N) {
      // scale theta to become a probability vector
      log_sum_exp_theta = log(exp(theta1[n]) + exp(theta2[n]));
      theta1[n] = theta1[n] - log_sum_exp_theta;
      theta2[n] = theta2[n] - log_sum_exp_theta;
    }
    // likelihood of the mixture model
    for (n in 1:N) {
      array[2] real ps;
      ps[1] = theta1[n] + unwrap_von_mises_lpdf(Y[n] | mu1[n], kappa1[n]);
      ps[2] = theta2[n] + unwrap_von_mises_lpdf(Y[n] | mu2[n], kappa2[n]);
      target += log_sum_exp(ps);
    }
  }
  // priors including constants
  target += lprior;
  target += std_normal_lpdf(to_vector(z_1));
  target += std_normal_lpdf(to_vector(z_2));
  target += std_normal_lpdf(to_vector(z_3));
  target += std_normal_lpdf(to_vector(z_4));
  target += std_normal_lpdf(to_vector(z_5));
}
generated quantities {
  // actual population-level intercept
  real b_mu1_Intercept = Intercept_mu1 - dot_product(means_X_mu1, b_mu1);
  // actual population-level intercept
  real b_kappa1_Intercept = Intercept_kappa1 - dot_product(means_X_kappa1, b_kappa1);
  // actual population-level intercept
  real b_mu2_Intercept = Intercept_mu2 - dot_product(means_X_mu2, b_mu2);
  // actual population-level intercept
  real b_kappa2_Intercept = Intercept_kappa2 - dot_product(means_X_kappa2, b_kappa2);
  // actual population-level intercept
  real b_theta1_Intercept = Intercept_theta1 - dot_product(means_X_theta1, b_theta1);
  // compute group-level correlations
  corr_matrix[M_1] Cor_1 = multiply_lower_tri_self_transpose(L_1);
  vector<lower=-1,upper=1>[NC_1] cor_1;
  // compute group-level correlations
  corr_matrix[M_2] Cor_2 = multiply_lower_tri_self_transpose(L_2);
  vector<lower=-1,upper=1>[NC_2] cor_2;
  // compute group-level correlations
  corr_matrix[M_3] Cor_3 = multiply_lower_tri_self_transpose(L_3);
  vector<lower=-1,upper=1>[NC_3] cor_3;
  // compute group-level correlations
  corr_matrix[M_4] Cor_4 = multiply_lower_tri_self_transpose(L_4);
  vector<lower=-1,upper=1>[NC_4] cor_4;
  // compute group-level correlations
  corr_matrix[M_5] Cor_5 = multiply_lower_tri_self_transpose(L_5);
  vector<lower=-1,upper=1>[NC_5] cor_5;
  
  vector [Kc_mu1 +1] mu_circ1; //modulo circular estimate
  vector [Kc_mu2 +1] mu_circ2; //modulo circular estimate
  mu_circ1[1] = mod_circular(Intercept_mu1); //Intercept case
  mu_circ2[1] = mod_circular(Intercept_mu2); //Intercept case
  for (i in 1:Kc_mu1){
  mu_circ1[i+1] = mod_circular(b_mu1[i]);
  mu_circ2[i+1] = mod_circular(b_mu2[i]);
  }

  // extract upper diagonal of correlation matrix
  for (k in 1:M_1) {
    for (j in 1:(k - 1)) {
      cor_1[choose(k - 1, 2) + j] = Cor_1[j, k];
    }
  }
  // extract upper diagonal of correlation matrix
  for (k in 1:M_2) {
    for (j in 1:(k - 1)) {
      cor_2[choose(k - 1, 2) + j] = Cor_2[j, k];
    }
  }
  // extract upper diagonal of correlation matrix
  for (k in 1:M_3) {
    for (j in 1:(k - 1)) {
      cor_3[choose(k - 1, 2) + j] = Cor_3[j, k];
    }
  }
  // extract upper diagonal of correlation matrix
  for (k in 1:M_4) {
    for (j in 1:(k - 1)) {
      cor_4[choose(k - 1, 2) + j] = Cor_4[j, k];
    }
  }
  // extract upper diagonal of correlation matrix
  for (k in 1:M_5) {
    for (j in 1:(k - 1)) {
      cor_5[choose(k - 1, 2) + j] = Cor_5[j, k];
    }
  }
}

