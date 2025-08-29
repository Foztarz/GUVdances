// generated with brms 2.22.0
functions {
  /* compute the tan_half link
   * Args:
   *   x: a scalar in (-pi, pi)
   * Returns:
   *   a scalar in (-Inf, Inf)
   */
   real tan_half(real x) {
     return tan(x / 2);
   }
  /* compute the tan_half link (vectorized)
   * Args:
   *   x: a vector in (-pi, pi)
   * Returns:
   *   a vector in (-Inf, Inf)
   */
   vector tan_half(vector x) {
     return tan(x / 2);
   }
  /* compute the inverse of the tan_half link
   * Args:
   *   y: a scalar in (-Inf, Inf)
   * Returns:
   *   a scalar in (-pi, pi)
   */
   real inv_tan_half(real y) {
     return 2 * atan(y);
   }
  /* compute the inverse of the tan_half link (vectorized)
   * Args:
   *   y: a vector in (-Inf, Inf)
   * Returns:
   *   a vector in (-pi, pi)
   */
   vector inv_tan_half(vector y) {
     return 2 * atan(y);
   }

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
  
   // calculate the mean angle of a circular distribution (in radians)
   real mean_circular(vector y){
      int N = size(y);
      real sumsin = 0;
      real sumcos = 0;
     for (n in 1:N){
       sumsin += sin(y[n]);
       sumcos += cos(y[n]);
     }
     sumsin = sumsin/N;
     sumcos = sumcos/N;
     return(atan2(  sumsin, sumcos) );
   }

  
  real mod_circular(real y) {
    return fmod(y + pi(), 2*pi()) - pi();
  }

  
   // calculate the unwrapped version (no discontinuities)
vector unwrap_circular(vector y)
{
  int N = size(y);
  real mux = mean_circular(y);
  vector[N] centx = y - mux;
  for(n in 1:N)
  {
  centx[n] = mod_circular(centx[n]);
  }
  vector[N] unwrx = centx + mux;
  return(unwrx);
}

  
real von_mises3_lpdf(real y, real mu, real kappa) {
     if (kappa < 100) {
       return von_mises_lpdf(mod_circular(y) | mu, kappa);
     } else {
       return normal_lpdf(mod_circular(y) | mu, sqrt(1 / kappa));
     }
   }

  
real von_misesmix_lpdf(real y, real mu1, real kappa1, real mu2, real kappa2, real lambda) {
     if (kappa1 < 100 && kappa2 < 100) {
       return log( exp(von_mises_lpdf(mod_circular(y) | mu1, kappa1)+log(lambda)) + exp(von_mises_lpdf(mod_circular(y) | mu2, kappa2)+log(1-lambda)) );
     } else {
       return log( exp(normal_lpdf(mod_circular(y) | mu1, sqrt(1/kappa1))+log(lambda)) + exp(normal_lpdf(mod_circular(y) | mu2, 1/sqrt(kappa2))+log(1-lambda)) );
     }
   }

}
data {
  int<lower=1> N;  // total number of observations
  vector[N] Y;  // response variable
  int<lower=1> K_fmu;  // number of population-level effects
  matrix[N, K_fmu] X_fmu;  // population-level design matrix
  int<lower=1> K_zmu;  // number of population-level effects
  matrix[N, K_zmu] X_zmu;  // population-level design matrix
  int<lower=1> K_kappa;  // number of population-level effects
  matrix[N, K_kappa] X_kappa;  // population-level design matrix
  int<lower=1> Kc_kappa;  // number of population-level effects after centering
  // data for group-level effects of ID 1
  int<lower=1> N_1;  // number of grouping levels
  int<lower=1> M_1;  // number of coefficients per level
  array[N] int<lower=1> J_1;  // grouping indicator per observation
  // group-level predictor values
  vector[N] Z_1_kappa_1;
  vector[N] Z_1_kappa_2;
  vector[N] Z_1_kappa_3;
  vector[N] Z_1_kappa_4;
  int<lower=1> NC_1;  // number of group-level correlations
  int prior_only;  // should the likelihood be ignored?
}
transformed data {
  matrix[N, Kc_kappa] Xc_kappa;  // centered version of X_kappa without an intercept
  vector[Kc_kappa] means_X_kappa;  // column means of X_kappa before centering
  for (i in 2:K_kappa) {
    means_X_kappa[i - 1] = mean(X_kappa[, i]);
    Xc_kappa[, i - 1] = X_kappa[, i] - means_X_kappa[i - 1];
  }
}
parameters {
  vector[K_fmu] b_fmu;  // regression coefficients
  vector[K_zmu] b_zmu;  // regression coefficients
  vector[Kc_kappa] b_kappa;  // regression coefficients
  real Intercept_kappa;  // temporary intercept for centered predictors
  vector<lower=0>[M_1] sd_1;  // group-level standard deviations
  matrix[M_1, N_1] z_1;  // standardized group-level effects
  cholesky_factor_corr[M_1] L_1;  // cholesky factor of correlation matrix
  
real zkappa1;
real zkappa2;
real zkappa3;
real zkappa4;
                           
}
transformed parameters {
  matrix[N_1, M_1] r_1;  // actual group-level effects
  // using vectors speeds up indexing in loops
  vector[N_1] r_1_kappa_1;
  vector[N_1] r_1_kappa_2;
  vector[N_1] r_1_kappa_3;
  vector[N_1] r_1_kappa_4;
  real lprior = 0;  // prior contributions to the log posterior
  // compute actual group-level effects
  r_1 = scale_r_cor(z_1, sd_1, L_1);
  r_1_kappa_1 = r_1[, 1];
  r_1_kappa_2 = r_1[, 2];
  r_1_kappa_3 = r_1[, 3];
  r_1_kappa_4 = r_1[, 4];
  lprior += normal_lpdf(b_fmu[1] | 0, pi()/3);
  lprior += normal_lpdf(b_fmu[2] | 0, pi()/3);
  lprior += normal_lpdf(b_fmu[3] | 0, pi()/3);
  lprior += normal_lpdf(b_fmu[4] | 0, pi()/3);
  lprior += von_mises3_lpdf(b_zmu[1] | 0, log1p_exp(zkappa1));
  lprior += von_mises3_lpdf(b_zmu[2] | 0, log1p_exp(zkappa1));
  lprior += von_mises3_lpdf(b_zmu[3] | 0, log1p_exp(zkappa1));
  lprior += von_mises3_lpdf(b_zmu[4] | 0, log1p_exp(zkappa1));
  lprior += von_mises3_lpdf(b_zmu[5] | 0, log1p_exp(zkappa1));
  lprior += von_mises3_lpdf(b_zmu[6] | 0, log1p_exp(zkappa1));
  lprior += von_mises3_lpdf(b_zmu[7] | 0, log1p_exp(zkappa1));
  lprior += von_mises3_lpdf(b_zmu[8] | 0, log1p_exp(zkappa1));
  lprior += von_mises3_lpdf(b_zmu[9] | 0, log1p_exp(zkappa1));
  lprior += von_mises3_lpdf(b_zmu[10] | 0, log1p_exp(zkappa1));
  lprior += von_mises3_lpdf(b_zmu[11] | 0, log1p_exp(zkappa1));
  lprior += von_mises3_lpdf(b_zmu[12] | 0, log1p_exp(zkappa1));
  lprior += von_mises3_lpdf(b_zmu[13] | 0, log1p_exp(zkappa1));
  lprior += von_mises3_lpdf(b_zmu[14] | 0, log1p_exp(zkappa1));
  lprior += von_mises3_lpdf(b_zmu[15] | 0, log1p_exp(zkappa1));
  lprior += von_mises3_lpdf(b_zmu[16] | 0, log1p_exp(zkappa1));
  lprior += von_mises3_lpdf(b_zmu[17] | 0, log1p_exp(zkappa1));
  lprior += von_mises3_lpdf(b_zmu[18] | 0, log1p_exp(zkappa1+zkappa2));
  lprior += von_mises3_lpdf(b_zmu[19] | 0, log1p_exp(zkappa1));
  lprior += von_mises3_lpdf(b_zmu[20] | 0, log1p_exp(zkappa1+zkappa2));
  lprior += von_mises3_lpdf(b_zmu[21] | 0, log1p_exp(zkappa1));
  lprior += von_mises3_lpdf(b_zmu[22] | 0, log1p_exp(zkappa1+zkappa2));
  lprior += von_mises3_lpdf(b_zmu[23] | 0, log1p_exp(zkappa1));
  lprior += von_mises3_lpdf(b_zmu[24] | 0, log1p_exp(zkappa1+zkappa2));
  lprior += von_mises3_lpdf(b_zmu[25] | 0, log1p_exp(zkappa1+zkappa2));
  lprior += von_mises3_lpdf(b_zmu[26] | 0, log1p_exp(zkappa1+zkappa2));
  lprior += von_mises3_lpdf(b_zmu[27] | 0, log1p_exp(zkappa1));
  lprior += von_mises3_lpdf(b_zmu[28] | 0, log1p_exp(zkappa1+zkappa2));
  lprior += von_mises3_lpdf(b_zmu[29] | 0, log1p_exp(zkappa1));
  lprior += von_mises3_lpdf(b_zmu[30] | 0, log1p_exp(zkappa1));
  lprior += von_mises3_lpdf(b_zmu[31] | 0, log1p_exp(zkappa1+zkappa2));
  lprior += von_mises3_lpdf(b_zmu[32] | 0, log1p_exp(zkappa1+zkappa2));
  lprior += von_mises3_lpdf(b_zmu[33] | 0, log1p_exp(zkappa1));
  lprior += von_mises3_lpdf(b_zmu[34] | 0, log1p_exp(zkappa1+zkappa2));
  lprior += von_mises3_lpdf(b_zmu[35] | 0, log1p_exp(zkappa1+zkappa3));
  lprior += von_mises3_lpdf(b_zmu[36] | 0, log1p_exp(zkappa1));
  lprior += von_mises3_lpdf(b_zmu[37] | 0, log1p_exp(zkappa1+zkappa3));
  lprior += von_mises3_lpdf(b_zmu[38] | 0, log1p_exp(zkappa1));
  lprior += von_mises3_lpdf(b_zmu[39] | 0, log1p_exp(zkappa1+zkappa3));
  lprior += von_mises3_lpdf(b_zmu[40] | 0, log1p_exp(zkappa1));
  lprior += von_mises3_lpdf(b_zmu[41] | 0, log1p_exp(zkappa1+zkappa3));
  lprior += von_mises3_lpdf(b_zmu[42] | 0, log1p_exp(zkappa1+zkappa3));
  lprior += von_mises3_lpdf(b_zmu[43] | 0, log1p_exp(zkappa1+zkappa3));
  lprior += von_mises3_lpdf(b_zmu[44] | 0, log1p_exp(zkappa1));
  lprior += von_mises3_lpdf(b_zmu[45] | 0, log1p_exp(zkappa1+zkappa3));
  lprior += von_mises3_lpdf(b_zmu[46] | 0, log1p_exp(zkappa1));
  lprior += von_mises3_lpdf(b_zmu[47] | 0, log1p_exp(zkappa1));
  lprior += von_mises3_lpdf(b_zmu[48] | 0, log1p_exp(zkappa1+zkappa3));
  lprior += von_mises3_lpdf(b_zmu[49] | 0, log1p_exp(zkappa1+zkappa3));
  lprior += von_mises3_lpdf(b_zmu[50] | 0, log1p_exp(zkappa1));
  lprior += von_mises3_lpdf(b_zmu[51] | 0, log1p_exp(zkappa1+zkappa3));
  lprior += von_mises3_lpdf(b_zmu[52] | 0, log1p_exp(zkappa1+zkappa2+zkappa3+zkappa4));
  lprior += von_mises3_lpdf(b_zmu[53] | 0, log1p_exp(zkappa1));
  lprior += von_mises3_lpdf(b_zmu[54] | 0, log1p_exp(zkappa1+zkappa2+zkappa3+zkappa4));
  lprior += von_mises3_lpdf(b_zmu[55] | 0, log1p_exp(zkappa1));
  lprior += von_mises3_lpdf(b_zmu[56] | 0, log1p_exp(zkappa1+zkappa2+zkappa3+zkappa4));
  lprior += von_mises3_lpdf(b_zmu[57] | 0, log1p_exp(zkappa1));
  lprior += von_mises3_lpdf(b_zmu[58] | 0, log1p_exp(zkappa1+zkappa2+zkappa3+zkappa4));
  lprior += von_mises3_lpdf(b_zmu[59] | 0, log1p_exp(zkappa1+zkappa2+zkappa3+zkappa4));
  lprior += von_mises3_lpdf(b_zmu[60] | 0, log1p_exp(zkappa1+zkappa2+zkappa3+zkappa4));
  lprior += von_mises3_lpdf(b_zmu[61] | 0, log1p_exp(zkappa1));
  lprior += von_mises3_lpdf(b_zmu[62] | 0, log1p_exp(zkappa1+zkappa2+zkappa3+zkappa4));
  lprior += von_mises3_lpdf(b_zmu[63] | 0, log1p_exp(zkappa1));
  lprior += von_mises3_lpdf(b_zmu[64] | 0, log1p_exp(zkappa1));
  lprior += von_mises3_lpdf(b_zmu[65] | 0, log1p_exp(zkappa1+zkappa2+zkappa3+zkappa4));
  lprior += von_mises3_lpdf(b_zmu[66] | 0, log1p_exp(zkappa1+zkappa2+zkappa3+zkappa4));
  lprior += von_mises3_lpdf(b_zmu[67] | 0, log1p_exp(zkappa1));
  lprior += von_mises3_lpdf(b_zmu[68] | 0, log1p_exp(zkappa1+zkappa2+zkappa3+zkappa4));
  lprior += normal_lpdf(b_kappa[1] | 0.0, 2.0);
  lprior += normal_lpdf(b_kappa[2] | 0.0, 2.0);
  lprior += normal_lpdf(b_kappa[3] | 0.0, 2.0);
  lprior += normal_lpdf(Intercept_kappa | 3.0, 3.0);
  lprior += student_t_lpdf(sd_1[1] | 3, 0, 2.0)
    - 1 * student_t_lccdf(0 | 3, 0, 2.0);
  lprior += student_t_lpdf(sd_1[2] | 3, 0, 2.0)
    - 1 * student_t_lccdf(0 | 3, 0, 2.0);
  lprior += student_t_lpdf(sd_1[3] | 3, 0, 2.0)
    - 1 * student_t_lccdf(0 | 3, 0, 2.0);
  lprior += student_t_lpdf(sd_1[4] | 3, 0, 2.0)
    - 1 * student_t_lccdf(0 | 3, 0, 2.0);
  lprior += lkj_corr_cholesky_lpdf(L_1 | 1);
}
model {
  // likelihood including constants
  if (!prior_only) {
    // initialize linear predictor term
    vector[N] nlp_fmu = rep_vector(0.0, N);
    // initialize linear predictor term
    vector[N] nlp_zmu = rep_vector(0.0, N);
    // initialize non-linear predictor term
    vector[N] mu;
    // initialize linear predictor term
    vector[N] kappa = rep_vector(0.0, N);
    nlp_fmu += X_fmu * b_fmu;
    nlp_zmu += X_zmu * b_zmu;
    kappa += Intercept_kappa + Xc_kappa * b_kappa;
    for (n in 1:N) {
      // add more terms to the linear predictor
      kappa[n] += r_1_kappa_1[J_1[n]] * Z_1_kappa_1[n] + r_1_kappa_2[J_1[n]] * Z_1_kappa_2[n] + r_1_kappa_3[J_1[n]] * Z_1_kappa_3[n] + r_1_kappa_4[J_1[n]] * Z_1_kappa_4[n];
    }
    for (n in 1:N) {
      // compute non-linear predictor values
      mu[n] = (mod_circular(nlp_fmu[n] + nlp_zmu[n]));
    }
    kappa = log1p_exp(kappa);
    target += von_mises_lpdf(Y | mu, kappa);
  }
  // priors including constants
  target += lprior;
  target += std_normal_lpdf(to_vector(z_1));
  target += normal_lpdf(zkappa1 | 3, 3);
  target += normal_lpdf(zkappa2 | 0, 0.5);
  target += normal_lpdf(zkappa3 | 0, 0.5);
  target += normal_lpdf(zkappa4 | 0, 0.3);
}
generated quantities {
  // actual population-level intercept
  real b_kappa_Intercept = Intercept_kappa - dot_product(means_X_kappa, b_kappa);
  // compute group-level correlations
  corr_matrix[M_1] Cor_1 = multiply_lower_tri_self_transpose(L_1);
  vector<lower=-1,upper=1>[NC_1] cor_1;
  
vector [M_1] mu_circ; //modulo circular estimate
for (i in 1:size(b_fmu)){
mu_circ[i] = mod_circular(b_fmu[i]);
}

  
real kappa_id_condition1 = log1p_exp(zkappa1);
real kappa_id_condition2 = log1p_exp(zkappa1+zkappa2);
real kappa_id_condition3 = log1p_exp(zkappa1+zkappa3);
real kappa_id_condition4 = log1p_exp(zkappa1+zkappa2+zkappa3+zkappa4);
          
  
vector[N_1] zmu_id_condition1;  // condition on coefficients per indiv;
vector[N_1] zmu_id_condition2;  // condition on coefficients per indiv;
vector[N_1] zmu_id_condition3;  // condition on coefficients per indiv;
vector[N_1] zmu_id_condition4;  // condition on coefficients per indiv;
for (i in 1:N_1){
zmu_id_condition1[i] = mod_circular(b_zmu[i]); //each in modulus format
zmu_id_condition2[i] = mod_circular(b_zmu[i+N_1]); //each in modulus format
zmu_id_condition3[i] = mod_circular(b_zmu[i+N_1*2]); //each in modulus format
zmu_id_condition4[i] = mod_circular(b_zmu[i+N_1*3]); //each in modulus format
}
          
  // extract upper diagonal of correlation matrix
  for (k in 1:M_1) {
    for (j in 1:(k - 1)) {
      cor_1[choose(k - 1, 2) + j] = Cor_1[j, k];
    }
  }
}

