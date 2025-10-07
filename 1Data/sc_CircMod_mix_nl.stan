// generated with brms 2.23.0
functions {
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
  real unwrap_von_mises_vect_lpdf(vector y, real mu, real kappa) {
    real tmp = 0;
    for(i in 1:size(y))
    {
    tmp = tmp + unwrap_von_mises_lpdf(y[i] | mu, kappa);
    }
      return tmp;
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
  int<lower=1> K_fmu1;  // number of population-level effects
  matrix[N, K_fmu1] X_fmu1;  // population-level design matrix
  int<lower=1> K_zmu1;  // number of population-level effects
  matrix[N, K_zmu1] X_zmu1;  // population-level design matrix
  int<lower=1> K_fmu2;  // number of population-level effects
  matrix[N, K_fmu2] X_fmu2;  // population-level design matrix
  int<lower=1> K_zmu2;  // number of population-level effects
  matrix[N, K_zmu2] X_zmu2;  // population-level design matrix
  int<lower=1> K_k1;  // number of population-level effects
  matrix[N, K_k1] X_k1;  // population-level design matrix
  int<lower=1> K_theta1;  // number of population-level effects
  matrix[N, K_theta1] X_theta1;  // population-level design matrix
  int<lower=1> Kc_theta1;  // number of population-level effects after centering
  // data for group-level effects of ID 1
  int<lower=1> N_1;  // number of grouping levels
  int<lower=1> M_1;  // number of coefficients per level
  array[N] int<lower=1> J_1;  // grouping indicator per observation
  // group-level predictor values
  vector[N] Z_1_theta1_1;
  vector[N] Z_1_theta1_2;
  vector[N] Z_1_theta1_3;
  vector[N] Z_1_theta1_4;
  vector[N] Z_1_theta1_5;
  int<lower=1> NC_1;  // number of group-level correlations
  // data for group-level effects of ID 2
  int<lower=1> N_2;  // number of grouping levels
  int<lower=1> M_2;  // number of coefficients per level
  array[N] int<lower=1> J_2;  // grouping indicator per observation
  // group-level predictor values
  vector[N] Z_2_k1_1;
  vector[N] Z_2_k1_2;
  vector[N] Z_2_k1_3;
  vector[N] Z_2_k1_4;
  int<lower=1> NC_2;  // number of group-level correlations
  int prior_only;  // should the likelihood be ignored?
}
transformed data {
  matrix[N, Kc_theta1] Xc_theta1;  // centered version of X_theta1 without an intercept
  vector[Kc_theta1] means_X_theta1;  // column means of X_theta1 before centering
  for (i in 2:K_theta1) {
    means_X_theta1[i - 1] = mean(X_theta1[, i]);
    Xc_theta1[, i - 1] = X_theta1[, i] - means_X_theta1[i - 1];
  }
}
parameters {
  vector[K_fmu1] b_fmu1;  // regression coefficients
  vector[K_zmu1] b_zmu1;  // regression coefficients
  vector[K_fmu2] b_fmu2;  // regression coefficients
  vector[K_zmu2] b_zmu2;  // regression coefficients
  vector[K_k1] b_k1;  // regression coefficients
  vector[Kc_theta1] b_theta1;  // regression coefficients
  real Intercept_theta1;  // temporary intercept for centered predictors
  vector<lower=0>[M_1] sd_1;  // group-level standard deviations
  matrix[M_1, N_1] z_1;  // standardized group-level effects
  cholesky_factor_corr[M_1] L_1;  // cholesky factor of correlation matrix
  vector<lower=0>[M_2] sd_2;  // group-level standard deviations
  matrix[M_2, N_2] z_2;  // standardized group-level effects
  cholesky_factor_corr[M_2] L_2;  // cholesky factor of correlation matrix
  
real kappamu1;
real kappamu1BR;
real kappamu1CL;
real kappamu1BR_CL;
real kappamu2;
real kappamu2BR;
real kappamu2CL;
real kappamu2BR_CL;
                           
}
transformed parameters {
  matrix[N_1, M_1] r_1;  // actual group-level effects
  // using vectors speeds up indexing in loops
  vector[N_1] r_1_theta1_1;
  vector[N_1] r_1_theta1_2;
  vector[N_1] r_1_theta1_3;
  vector[N_1] r_1_theta1_4;
  vector[N_1] r_1_theta1_5;
  matrix[N_2, M_2] r_2;  // actual group-level effects
  // using vectors speeds up indexing in loops
  vector[N_2] r_2_k1_1;
  vector[N_2] r_2_k1_2;
  vector[N_2] r_2_k1_3;
  vector[N_2] r_2_k1_4;
  // prior contributions to the log posterior
  real lprior = 0;
  // compute actual group-level effects
  r_1 = scale_r_cor(z_1, sd_1, L_1);
  r_1_theta1_1 = r_1[, 1];
  r_1_theta1_2 = r_1[, 2];
  r_1_theta1_3 = r_1[, 3];
  r_1_theta1_4 = r_1[, 4];
  r_1_theta1_5 = r_1[, 5];
  // compute actual group-level effects
  r_2 = scale_r_cor(z_2, sd_2, L_2);
  r_2_k1_1 = r_2[, 1];
  r_2_k1_2 = r_2[, 2];
  r_2_k1_3 = r_2[, 3];
  r_2_k1_4 = r_2[, 4];
  lprior += normal_lpdf(b_fmu1[1] | 0, pi()/3);
  lprior += normal_lpdf(b_fmu1[2] | 0, pi()/3);
  lprior += normal_lpdf(b_fmu1[3] | 0, pi()/3);
  lprior += normal_lpdf(b_fmu1[4] | pi()/3, pi()/2);
  lprior += normal_lpdf(b_fmu2[1] | 0, pi()/6);
  lprior += normal_lpdf(b_fmu2[2] | 0, pi()/6);
  lprior += normal_lpdf(b_fmu2[3] | 0, pi()/6);
  lprior += normal_lpdf(b_fmu2[4] | -pi()/3, pi()/2);
  lprior += normal_lpdf(b_k1[1] | 3, 3);
  lprior += normal_lpdf(b_k1[2] | 0, 3);
  lprior += normal_lpdf(b_k1[3] | 0, 3);
  lprior += normal_lpdf(b_k1[4] | 0, 3);
  lprior += normal_lpdf(b_theta1[1] | 0, 1);
  lprior += normal_lpdf(b_theta1[2] | 0, 1);
  lprior += normal_lpdf(b_theta1[3] | 0, 1);
  lprior += normal_lpdf(b_theta1[4] | 0, 5);
  lprior += normal_lpdf(Intercept_theta1 | 5, 2);
  lprior += student_t_lpdf(sd_1 | 3, 0, 0.5)
    - 5 * student_t_lccdf(0 | 3, 0, 0.5);
  lprior += lkj_corr_cholesky_lpdf(L_1 | 1);
  lprior += student_t_lpdf(sd_2 | 3, 0, 1)
    - 4 * student_t_lccdf(0 | 3, 0, 1);
  lprior += lkj_corr_cholesky_lpdf(L_2 | 1);
}
model {
  // likelihood including constants
  if (!prior_only) {
    // initialize linear predictor term
    vector[N] nlp_fmu1 = rep_vector(0.0, N);
    // initialize linear predictor term
    vector[N] nlp_zmu1 = rep_vector(0.0, N);
    // initialize linear predictor term
    vector[N] nlp_fmu2 = rep_vector(0.0, N);
    // initialize linear predictor term
    vector[N] nlp_zmu2 = rep_vector(0.0, N);
    // initialize linear predictor term
    vector[N] nlp_k1 = rep_vector(0.0, N);
    // initialize non-linear predictor term
    vector[N] mu1;
    // initialize non-linear predictor term
    vector[N] kappa1;
    // initialize non-linear predictor term
    vector[N] mu2;
    // initialize non-linear predictor term
    vector[N] kappa2;
    // initialize linear predictor term
    vector[N] theta1 = rep_vector(0.0, N);
    vector[N] theta2 = rep_vector(0.0, N);
    real log_sum_exp_theta;
    nlp_fmu1 += X_fmu1 * b_fmu1;
    nlp_zmu1 += X_zmu1 * b_zmu1;
    nlp_fmu2 += X_fmu2 * b_fmu2;
    nlp_zmu2 += X_zmu2 * b_zmu2;
    nlp_k1 += X_k1 * b_k1;
    theta1 += Intercept_theta1 + Xc_theta1 * b_theta1;
    for (n in 1:N) {
      // add more terms to the linear predictor
      nlp_k1[n] += r_2_k1_1[J_2[n]] * Z_2_k1_1[n] + r_2_k1_2[J_2[n]] * Z_2_k1_2[n] + r_2_k1_3[J_2[n]] * Z_2_k1_3[n] + r_2_k1_4[J_2[n]] * Z_2_k1_4[n];
    }
    for (n in 1:N) {
      // add more terms to the linear predictor
      theta1[n] += r_1_theta1_1[J_1[n]] * Z_1_theta1_1[n] + r_1_theta1_2[J_1[n]] * Z_1_theta1_2[n] + r_1_theta1_3[J_1[n]] * Z_1_theta1_3[n] + r_1_theta1_4[J_1[n]] * Z_1_theta1_4[n] + r_1_theta1_5[J_1[n]] * Z_1_theta1_5[n];
    }
    for (n in 1:N) {
      // compute non-linear predictor values
      mu1[n] = (nlp_fmu1[n] + nlp_zmu1[n]);
    }
    for (n in 1:N) {
      // compute non-linear predictor values
      kappa1[n] = log1p_exp(nlp_k1[n]);
    }
    for (n in 1:N) {
      // compute non-linear predictor values
      mu2[n] = (nlp_fmu1[n] + nlp_fmu2[n] + nlp_zmu1[n] + nlp_zmu2[n]);
    }
    for (n in 1:N) {
      // compute non-linear predictor values
      kappa2[n] = log1p_exp(nlp_k1[n]);
    }
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
  target += unwrap_von_mises_vect_lpdf(b_zmu1[1:10] | 0, log1p_exp(kappamu1)) + normal_lpdf(b_zmu1 | 0, 2*pi());
  target += normal_lpdf(kappamu1 | 3.0, 3.0);
  target += unwrap_von_mises_vect_lpdf(b_zmu2[1:10] | 0, log1p_exp(kappamu1 + kappamu2)) + normal_lpdf(b_zmu2 | 0, 2*pi());
  target += normal_lpdf(kappamu2 | 0.0, 1.0);
}
generated quantities {
  // actual population-level intercept
  real b_theta1_Intercept = Intercept_theta1 - dot_product(means_X_theta1, b_theta1);
  // compute group-level correlations
  corr_matrix[M_1] Cor_1 = multiply_lower_tri_self_transpose(L_1);
  vector<lower=-1,upper=1>[NC_1] cor_1;
  // compute group-level correlations
  corr_matrix[M_2] Cor_2 = multiply_lower_tri_self_transpose(L_2);
  vector<lower=-1,upper=1>[NC_2] cor_2;
  
real kappa_mu1 = log1p_exp(kappamu1);
real kappa_mu2 = log1p_exp(kappamu1+kappamu2);
          
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
}

