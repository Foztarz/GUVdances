// generated with brms 2.23.0
functions {
  
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
  // data for group-level effects of ID 1
  int<lower=1> N_1;  // number of grouping levels
  int<lower=1> M_1;  // number of coefficients per level
  array[N] int<lower=1> J_1;  // grouping indicator per observation
  // group-level predictor values
  vector[N] Z_1_theta1_1;
  // data for group-level effects of ID 2
  int<lower=1> N_2;  // number of grouping levels
  int<lower=1> M_2;  // number of coefficients per level
  array[N] int<lower=1> J_2;  // grouping indicator per observation
  // group-level predictor values
  vector[N] Z_2_k1_1;
  int prior_only;  // should the likelihood be ignored?
}
transformed data {
}
parameters {
  vector[K_fmu1] b_fmu1;  // regression coefficients
  vector[K_zmu1] b_zmu1;  // regression coefficients
  vector[K_fmu2] b_fmu2;  // regression coefficients
  vector[K_zmu2] b_zmu2;  // regression coefficients
  vector[K_k1] b_k1;  // regression coefficients
  real Intercept_theta1;  // temporary intercept for centered predictors
  vector<lower=0>[M_1] sd_1;  // group-level standard deviations
  array[M_1] vector[N_1] z_1;  // standardized group-level effects
  vector<lower=0>[M_2] sd_2;  // group-level standard deviations
  array[M_2] vector[N_2] z_2;  // standardized group-level effects
  
real kappamu1;
real kappamu2;
                           
}
transformed parameters {
  vector[N_1] r_1_theta1_1;  // actual group-level effects
  vector[N_2] r_2_k1_1;  // actual group-level effects
  // prior contributions to the log posterior
  real lprior = 0;
  r_1_theta1_1 = (sd_1[1] * (z_1[1]));
  r_2_k1_1 = (sd_2[1] * (z_2[1]));
  lprior += logistic_lpdf(Intercept_theta1 | 0, 1);
  lprior += student_t_lpdf(sd_1 | 3, 0, 2.5)
    - 1 * student_t_lccdf(0 | 3, 0, 2.5);
  lprior += student_t_lpdf(sd_2 | 3, 0, 2.5)
    - 1 * student_t_lccdf(0 | 3, 0, 2.5);
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
    theta1 += Intercept_theta1;
    for (n in 1:N) {
      // add more terms to the linear predictor
      nlp_k1[n] += r_2_k1_1[J_2[n]] * Z_2_k1_1[n];
    }
    for (n in 1:N) {
      // add more terms to the linear predictor
      theta1[n] += r_1_theta1_1[J_1[n]] * Z_1_theta1_1[n];
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
  target += std_normal_lpdf(z_1[1]);
  target += std_normal_lpdf(z_2[1]);
}
generated quantities {
  // actual population-level intercept
  real b_theta1_Intercept = Intercept_theta1;
  
real kappa_mu1 = log1p_exp(kappamu1);
real kappa_mu2 = log1p_exp(kappamu1+kappamu2);
          
}

