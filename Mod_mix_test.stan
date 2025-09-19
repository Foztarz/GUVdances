// generated with brms 2.22.0
functions {
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
  int prior_only;  // should the likelihood be ignored?
}
transformed data {
}
parameters {
  vector[K_fmu1] b_fmu1;  // regression coefficients
  vector[K_zmu1] b_zmu1;  // regression coefficients
  vector[K_fmu2] b_fmu2;  // regression coefficients
  vector[K_zmu2] b_zmu2;  // regression coefficients
  real<lower=0> sigma1;  // dispersion parameter
  real<lower=0> sigma2;  // dispersion parameter
  real Intercept_theta1;  // temporary intercept for centered predictors
  ordered[2] ordered_Intercept;  // to identify mixtures
}
transformed parameters {
  real lprior = 0;  // prior contributions to the log posterior
  lprior += student_t_lpdf(sigma1 | 3, 0, 2.5)
    - 1 * student_t_lccdf(0 | 3, 0, 2.5);
  lprior += student_t_lpdf(sigma2 | 3, 0, 2.5)
    - 1 * student_t_lccdf(0 | 3, 0, 2.5);
  lprior += logistic_lpdf(Intercept_theta1 | 0, 1);
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
    // initialize non-linear predictor term
    vector[N] mu1;
    // initialize non-linear predictor term
    vector[N] mu2;
    // initialize linear predictor term
    vector[N] theta1 = rep_vector(0.0, N);
    vector[N] theta2 = rep_vector(0.0, N);
    real log_sum_exp_theta;
    nlp_fmu1 += X_fmu1 * b_fmu1;
    nlp_zmu1 += X_zmu1 * b_zmu1;
    nlp_fmu2 += X_fmu2 * b_fmu2;
    nlp_zmu2 += X_zmu2 * b_zmu2;
    theta1 += Intercept_theta1;
    for (n in 1:N) {
      // compute non-linear predictor values
      mu1[n] = (mi(nlp_fmu1[n] + nlp_zmu1[n]) + mi(nlp_fmu2[n] + nlp_zmu2[n]));
    }
    for (n in 1:N) {
      // compute non-linear predictor values
      mu2[n] = (mi(nlp_fmu1[n] + nlp_zmu1[n]) + mi(nlp_fmu2[n] + nlp_zmu2[n]));
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
      ps[1] = theta1[n] + normal_lpdf(Y[n] | mu1[n], sigma1);
      ps[2] = theta2[n] + normal_lpdf(Y[n] | mu2[n], sigma2);
      target += log_sum_exp(ps);
    }
  }
  // priors including constants
  target += lprior;
}
generated quantities {
  // actual population-level intercept
  real b_theta1_Intercept = Intercept_theta1;
}

