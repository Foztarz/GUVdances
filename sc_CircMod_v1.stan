// generated with brms 2.21.0
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

  /* von Mises log-PDF of a single response
   * for kappa > 100 the normal approximation is used
   * for reasons of numerial stability
   * Args:
   *   y: the response vector between -pi and pi
   *   mu: location parameter vector
   *   kappa: precision parameter
   * Returns:
   *   a scalar to be added to the log posterior
   */
   real von_mises2_lpdf(real y, real mu, real kappa) {
     if (kappa < 100) {
       return von_mises_lpdf(y | mu, kappa);
     } else {
       return normal_lpdf(y | mu, sqrt(1 / kappa));
     }
   }
  /* von Mises log-PDF of a response vector
   * for kappa > 100 the normal approximation is used
   * for reasons of numerial stability
   * Args:
   *   y: the response vector between -pi and pi
   *   mu: location parameter vector
   *   kappa: precision parameter
   * Returns:
   *   a scalar to be added to the log posterior
   */
   real von_mises2_lpdf(vector y, vector mu, real kappa) {
     if (kappa < 100) {
       return von_mises_lpdf(y | mu, kappa);
     } else {
       return normal_lpdf(y | mu, sqrt(1 / kappa));
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
  lprior += von_mises3_lpdf(b_fmu[1] | 0, 1.0);
  lprior += von_mises3_lpdf(b_fmu[2] | 0, 1.0);
  lprior += von_mises3_lpdf(b_fmu[3] | 0, 1.0);
  lprior += von_mises3_lpdf(b_fmu[4] | 0, 1.0);
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
  lprior += von_mises3_lpdf(b_zmu[18] | 0, log1p_exp(zkappa1));
  lprior += von_mises3_lpdf(b_zmu[19] | 0, log1p_exp(zkappa1));
  lprior += von_mises3_lpdf(b_zmu[20] | 0, log1p_exp(zkappa1));
  lprior += von_mises3_lpdf(b_zmu[21] | 0, log1p_exp(zkappa1));
  lprior += von_mises3_lpdf(b_zmu[22] | 0, log1p_exp(zkappa1));
  lprior += von_mises3_lpdf(b_zmu[23] | 0, log1p_exp(zkappa1));
  lprior += von_mises3_lpdf(b_zmu[24] | 0, log1p_exp(zkappa1));
  lprior += von_mises3_lpdf(b_zmu[25] | 0, log1p_exp(zkappa1));
  lprior += von_mises3_lpdf(b_zmu[26] | 0, log1p_exp(zkappa1));
  lprior += von_mises3_lpdf(b_zmu[27] | 0, log1p_exp(zkappa1));
  lprior += von_mises3_lpdf(b_zmu[28] | 0, log1p_exp(zkappa1));
  lprior += von_mises3_lpdf(b_zmu[29] | 0, log1p_exp(zkappa1));
  lprior += von_mises3_lpdf(b_zmu[30] | 0, log1p_exp(zkappa1));
  lprior += von_mises3_lpdf(b_zmu[31] | 0, log1p_exp(zkappa1));
  lprior += von_mises3_lpdf(b_zmu[32] | 0, log1p_exp(zkappa1));
  lprior += von_mises3_lpdf(b_zmu[33] | 0, log1p_exp(zkappa1));
  lprior += von_mises3_lpdf(b_zmu[34] | 0, log1p_exp(zkappa1));
  lprior += von_mises3_lpdf(b_zmu[35] | 0, log1p_exp(zkappa1));
  lprior += von_mises3_lpdf(b_zmu[36] | 0, log1p_exp(zkappa1));
  lprior += von_mises3_lpdf(b_zmu[37] | 0, log1p_exp(zkappa1));
  lprior += von_mises3_lpdf(b_zmu[38] | 0, log1p_exp(zkappa1));
  lprior += von_mises3_lpdf(b_zmu[39] | 0, log1p_exp(zkappa1));
  lprior += von_mises3_lpdf(b_zmu[40] | 0, log1p_exp(zkappa1));
  lprior += von_mises3_lpdf(b_zmu[41] | 0, log1p_exp(zkappa1));
  lprior += von_mises3_lpdf(b_zmu[42] | 0, log1p_exp(zkappa1));
  lprior += von_mises3_lpdf(b_zmu[43] | 0, log1p_exp(zkappa1));
  lprior += von_mises3_lpdf(b_zmu[44] | 0, log1p_exp(zkappa1));
  lprior += von_mises3_lpdf(b_zmu[45] | 0, log1p_exp(zkappa1));
  lprior += von_mises3_lpdf(b_zmu[46] | 0, log1p_exp(zkappa1));
  lprior += von_mises3_lpdf(b_zmu[47] | 0, log1p_exp(zkappa1));
  lprior += von_mises3_lpdf(b_zmu[48] | 0, log1p_exp(zkappa1));
  lprior += von_mises3_lpdf(b_zmu[49] | 0, log1p_exp(zkappa1));
  lprior += von_mises3_lpdf(b_zmu[50] | 0, log1p_exp(zkappa1));
  lprior += von_mises3_lpdf(b_zmu[51] | 0, log1p_exp(zkappa1));
  lprior += von_mises3_lpdf(b_zmu[52] | 0, log1p_exp(zkappa1));
  lprior += von_mises3_lpdf(b_zmu[53] | 0, log1p_exp(zkappa1));
  lprior += von_mises3_lpdf(b_zmu[54] | 0, log1p_exp(zkappa1));
  lprior += von_mises3_lpdf(b_zmu[55] | 0, log1p_exp(zkappa1));
  lprior += von_mises3_lpdf(b_zmu[56] | 0, log1p_exp(zkappa1));
  lprior += von_mises3_lpdf(b_zmu[57] | 0, log1p_exp(zkappa1));
  lprior += von_mises3_lpdf(b_zmu[58] | 0, log1p_exp(zkappa1));
  lprior += von_mises3_lpdf(b_zmu[59] | 0, log1p_exp(zkappa1));
  lprior += von_mises3_lpdf(b_zmu[60] | 0, log1p_exp(zkappa1));
  lprior += von_mises3_lpdf(b_zmu[61] | 0, log1p_exp(zkappa1));
  lprior += von_mises3_lpdf(b_zmu[62] | 0, log1p_exp(zkappa1));
  lprior += von_mises3_lpdf(b_zmu[63] | 0, log1p_exp(zkappa1));
  lprior += von_mises3_lpdf(b_zmu[64] | 0, log1p_exp(zkappa1));
  lprior += von_mises3_lpdf(b_zmu[65] | 0, log1p_exp(zkappa1));
  lprior += von_mises3_lpdf(b_zmu[66] | 0, log1p_exp(zkappa1));
  lprior += von_mises3_lpdf(b_zmu[67] | 0, log1p_exp(zkappa1));
  lprior += von_mises3_lpdf(b_zmu[68] | 0, log1p_exp(zkappa1));
  lprior += von_mises3_lpdf(b_zmu[69] | 0, log1p_exp(zkappa1));
  lprior += von_mises3_lpdf(b_zmu[70] | 0, log1p_exp(zkappa1));
  lprior += von_mises3_lpdf(b_zmu[71] | 0, log1p_exp(zkappa1));
  lprior += von_mises3_lpdf(b_zmu[72] | 0, log1p_exp(zkappa1));
  lprior += von_mises3_lpdf(b_zmu[73] | 0, log1p_exp(zkappa1));
  lprior += von_mises3_lpdf(b_zmu[74] | 0, log1p_exp(zkappa1));
  lprior += von_mises3_lpdf(b_zmu[75] | 0, log1p_exp(zkappa1));
  lprior += von_mises3_lpdf(b_zmu[76] | 0, log1p_exp(zkappa1));
  lprior += von_mises3_lpdf(b_zmu[77] | 0, log1p_exp(zkappa1));
  lprior += von_mises3_lpdf(b_zmu[78] | 0, log1p_exp(zkappa1));
  lprior += von_mises3_lpdf(b_zmu[79] | 0, log1p_exp(zkappa1));
  lprior += von_mises3_lpdf(b_zmu[80] | 0, log1p_exp(zkappa1));
  lprior += von_mises3_lpdf(b_zmu[81] | 0, log1p_exp(zkappa1));
  lprior += von_mises3_lpdf(b_zmu[82] | 0, log1p_exp(zkappa1));
  lprior += von_mises3_lpdf(b_zmu[83] | 0, log1p_exp(zkappa1));
  lprior += von_mises3_lpdf(b_zmu[84] | 0, log1p_exp(zkappa1));
  lprior += von_mises3_lpdf(b_zmu[85] | 0, log1p_exp(zkappa1));
  lprior += von_mises3_lpdf(b_zmu[86] | 0, log1p_exp(zkappa1));
  lprior += von_mises3_lpdf(b_zmu[87] | 0, log1p_exp(zkappa1));
  lprior += von_mises3_lpdf(b_zmu[88] | 0, log1p_exp(zkappa1));
  lprior += von_mises3_lpdf(b_zmu[89] | 0, log1p_exp(zkappa1));
  lprior += von_mises3_lpdf(b_zmu[90] | 0, log1p_exp(zkappa1));
  lprior += von_mises3_lpdf(b_zmu[91] | 0, log1p_exp(zkappa1));
  lprior += von_mises3_lpdf(b_zmu[92] | 0, log1p_exp(zkappa1));
  lprior += von_mises3_lpdf(b_zmu[93] | 0, log1p_exp(zkappa1));
  lprior += von_mises3_lpdf(b_zmu[94] | 0, log1p_exp(zkappa1));
  lprior += von_mises3_lpdf(b_zmu[95] | 0, log1p_exp(zkappa1));
  lprior += von_mises3_lpdf(b_zmu[96] | 0, log1p_exp(zkappa1));
  lprior += von_mises3_lpdf(b_zmu[97] | 0, log1p_exp(zkappa1));
  lprior += von_mises3_lpdf(b_zmu[98] | 0, log1p_exp(zkappa1));
  lprior += von_mises3_lpdf(b_zmu[99] | 0, log1p_exp(zkappa1));
  lprior += von_mises3_lpdf(b_zmu[100] | 0, log1p_exp(zkappa1));
  lprior += von_mises3_lpdf(b_zmu[101] | 0, log1p_exp(zkappa1));
  lprior += von_mises3_lpdf(b_zmu[102] | 0, log1p_exp(zkappa1));
  lprior += von_mises3_lpdf(b_zmu[103] | 0, log1p_exp(zkappa1));
  lprior += von_mises3_lpdf(b_zmu[104] | 0, log1p_exp(zkappa1));
  lprior += von_mises3_lpdf(b_zmu[105] | 0, log1p_exp(zkappa1));
  lprior += von_mises3_lpdf(b_zmu[106] | 0, log1p_exp(zkappa1));
  lprior += von_mises3_lpdf(b_zmu[107] | 0, log1p_exp(zkappa1));
  lprior += von_mises3_lpdf(b_zmu[108] | 0, log1p_exp(zkappa1));
  lprior += von_mises3_lpdf(b_zmu[109] | 0, log1p_exp(zkappa1));
  lprior += von_mises3_lpdf(b_zmu[110] | 0, log1p_exp(zkappa1));
  lprior += von_mises3_lpdf(b_zmu[111] | 0, log1p_exp(zkappa1));
  lprior += von_mises3_lpdf(b_zmu[112] | 0, log1p_exp(zkappa1));
  lprior += von_mises3_lpdf(b_zmu[113] | 0, log1p_exp(zkappa1));
  lprior += von_mises3_lpdf(b_zmu[114] | 0, log1p_exp(zkappa1));
  lprior += von_mises3_lpdf(b_zmu[115] | 0, log1p_exp(zkappa1));
  lprior += von_mises3_lpdf(b_zmu[116] | 0, log1p_exp(zkappa1));
  lprior += von_mises3_lpdf(b_zmu[117] | 0, log1p_exp(zkappa1));
  lprior += von_mises3_lpdf(b_zmu[118] | 0, log1p_exp(zkappa1));
  lprior += von_mises3_lpdf(b_zmu[119] | 0, log1p_exp(zkappa1));
  lprior += von_mises3_lpdf(b_zmu[120] | 0, log1p_exp(zkappa1));
  lprior += von_mises3_lpdf(b_zmu[121] | 0, log1p_exp(zkappa1));
  lprior += von_mises3_lpdf(b_zmu[122] | 0, log1p_exp(zkappa1));
  lprior += von_mises3_lpdf(b_zmu[123] | 0, log1p_exp(zkappa1));
  lprior += von_mises3_lpdf(b_zmu[124] | 0, log1p_exp(zkappa1));
  lprior += von_mises3_lpdf(b_zmu[125] | 0, log1p_exp(zkappa1));
  lprior += von_mises3_lpdf(b_zmu[126] | 0, log1p_exp(zkappa1));
  lprior += von_mises3_lpdf(b_zmu[127] | 0, log1p_exp(zkappa1));
  lprior += von_mises3_lpdf(b_zmu[128] | 0, log1p_exp(zkappa1));
  lprior += von_mises3_lpdf(b_zmu[129] | 0, log1p_exp(zkappa1));
  lprior += von_mises3_lpdf(b_zmu[130] | 0, log1p_exp(zkappa1));
  lprior += von_mises3_lpdf(b_zmu[131] | 0, log1p_exp(zkappa1));
  lprior += von_mises3_lpdf(b_zmu[132] | 0, log1p_exp(zkappa1));
  lprior += von_mises3_lpdf(b_zmu[133] | 0, log1p_exp(zkappa1));
  lprior += von_mises3_lpdf(b_zmu[134] | 0, log1p_exp(zkappa1));
  lprior += von_mises3_lpdf(b_zmu[135] | 0, log1p_exp(zkappa1));
  lprior += von_mises3_lpdf(b_zmu[136] | 0, log1p_exp(zkappa1));
  lprior += von_mises3_lpdf(b_zmu[137] | 0, log1p_exp(zkappa1));
  lprior += von_mises3_lpdf(b_zmu[138] | 0, log1p_exp(zkappa1));
  lprior += von_mises3_lpdf(b_zmu[139] | 0, log1p_exp(zkappa1));
  lprior += von_mises3_lpdf(b_zmu[140] | 0, log1p_exp(zkappa1));
  lprior += von_mises3_lpdf(b_zmu[141] | 0, log1p_exp(zkappa1));
  lprior += von_mises3_lpdf(b_zmu[142] | 0, log1p_exp(zkappa1));
  lprior += von_mises3_lpdf(b_zmu[143] | 0, log1p_exp(zkappa1));
  lprior += von_mises3_lpdf(b_zmu[144] | 0, log1p_exp(zkappa1));
  lprior += von_mises3_lpdf(b_zmu[145] | 0, log1p_exp(zkappa1));
  lprior += von_mises3_lpdf(b_zmu[146] | 0, log1p_exp(zkappa1));
  lprior += von_mises3_lpdf(b_zmu[147] | 0, log1p_exp(zkappa1));
  lprior += von_mises3_lpdf(b_zmu[148] | 0, log1p_exp(zkappa1));
  lprior += von_mises3_lpdf(b_zmu[149] | 0, log1p_exp(zkappa1));
  lprior += von_mises3_lpdf(b_zmu[150] | 0, log1p_exp(zkappa1));
  lprior += von_mises3_lpdf(b_zmu[151] | 0, log1p_exp(zkappa1));
  lprior += von_mises3_lpdf(b_zmu[152] | 0, log1p_exp(zkappa1));
  lprior += von_mises3_lpdf(b_zmu[153] | 0, log1p_exp(zkappa1));
  lprior += von_mises3_lpdf(b_zmu[154] | 0, log1p_exp(zkappa1));
  lprior += von_mises3_lpdf(b_zmu[155] | 0, log1p_exp(zkappa1));
  lprior += von_mises3_lpdf(b_zmu[156] | 0, log1p_exp(zkappa1));
  lprior += von_mises3_lpdf(b_zmu[157] | 0, log1p_exp(zkappa1));
  lprior += von_mises3_lpdf(b_zmu[158] | 0, log1p_exp(zkappa1));
  lprior += von_mises3_lpdf(b_zmu[159] | 0, log1p_exp(zkappa1));
  lprior += von_mises3_lpdf(b_zmu[160] | 0, log1p_exp(zkappa1));
  lprior += von_mises3_lpdf(b_zmu[161] | 0, log1p_exp(zkappa1));
  lprior += von_mises3_lpdf(b_zmu[162] | 0, log1p_exp(zkappa1));
  lprior += von_mises3_lpdf(b_zmu[163] | 0, log1p_exp(zkappa1));
  lprior += von_mises3_lpdf(b_zmu[164] | 0, log1p_exp(zkappa1));
  lprior += von_mises3_lpdf(b_zmu[165] | 0, log1p_exp(zkappa1));
  lprior += von_mises3_lpdf(b_zmu[166] | 0, log1p_exp(zkappa1));
  lprior += von_mises3_lpdf(b_zmu[167] | 0, log1p_exp(zkappa1));
  lprior += von_mises3_lpdf(b_zmu[168] | 0, log1p_exp(zkappa1));
  lprior += von_mises3_lpdf(b_zmu[169] | 0, log1p_exp(zkappa1));
  lprior += von_mises3_lpdf(b_zmu[170] | 0, log1p_exp(zkappa1));
  lprior += von_mises3_lpdf(b_zmu[171] | 0, log1p_exp(zkappa1));
  lprior += von_mises3_lpdf(b_zmu[172] | 0, log1p_exp(zkappa1));
  lprior += von_mises3_lpdf(b_zmu[173] | 0, log1p_exp(zkappa1));
  lprior += von_mises3_lpdf(b_zmu[174] | 0, log1p_exp(zkappa1));
  lprior += von_mises3_lpdf(b_zmu[175] | 0, log1p_exp(zkappa1));
  lprior += von_mises3_lpdf(b_zmu[176] | 0, log1p_exp(zkappa1));
  lprior += von_mises3_lpdf(b_zmu[177] | 0, log1p_exp(zkappa1));
  lprior += von_mises3_lpdf(b_zmu[178] | 0, log1p_exp(zkappa1));
  lprior += von_mises3_lpdf(b_zmu[179] | 0, log1p_exp(zkappa1));
  lprior += von_mises3_lpdf(b_zmu[180] | 0, log1p_exp(zkappa1));
  lprior += von_mises3_lpdf(b_zmu[181] | 0, log1p_exp(zkappa1));
  lprior += von_mises3_lpdf(b_zmu[182] | 0, log1p_exp(zkappa1));
  lprior += von_mises3_lpdf(b_zmu[183] | 0, log1p_exp(zkappa1));
  lprior += von_mises3_lpdf(b_zmu[184] | 0, log1p_exp(zkappa1));
  lprior += von_mises3_lpdf(b_zmu[185] | 0, log1p_exp(zkappa1));
  lprior += von_mises3_lpdf(b_zmu[186] | 0, log1p_exp(zkappa1));
  lprior += von_mises3_lpdf(b_zmu[187] | 0, log1p_exp(zkappa1));
  lprior += von_mises3_lpdf(b_zmu[188] | 0, log1p_exp(zkappa1+zkappa2));
  lprior += von_mises3_lpdf(b_zmu[189] | 0, log1p_exp(zkappa1));
  lprior += von_mises3_lpdf(b_zmu[190] | 0, log1p_exp(zkappa1));
  lprior += von_mises3_lpdf(b_zmu[191] | 0, log1p_exp(zkappa1+zkappa2));
  lprior += von_mises3_lpdf(b_zmu[192] | 0, log1p_exp(zkappa1));
  lprior += von_mises3_lpdf(b_zmu[193] | 0, log1p_exp(zkappa1+zkappa2));
  lprior += von_mises3_lpdf(b_zmu[194] | 0, log1p_exp(zkappa1));
  lprior += von_mises3_lpdf(b_zmu[195] | 0, log1p_exp(zkappa1));
  lprior += von_mises3_lpdf(b_zmu[196] | 0, log1p_exp(zkappa1));
  lprior += von_mises3_lpdf(b_zmu[197] | 0, log1p_exp(zkappa1));
  lprior += von_mises3_lpdf(b_zmu[198] | 0, log1p_exp(zkappa1));
  lprior += von_mises3_lpdf(b_zmu[199] | 0, log1p_exp(zkappa1));
  lprior += von_mises3_lpdf(b_zmu[200] | 0, log1p_exp(zkappa1));
  lprior += von_mises3_lpdf(b_zmu[201] | 0, log1p_exp(zkappa1));
  lprior += von_mises3_lpdf(b_zmu[202] | 0, log1p_exp(zkappa1));
  lprior += von_mises3_lpdf(b_zmu[203] | 0, log1p_exp(zkappa1));
  lprior += von_mises3_lpdf(b_zmu[204] | 0, log1p_exp(zkappa1));
  lprior += von_mises3_lpdf(b_zmu[205] | 0, log1p_exp(zkappa1));
  lprior += von_mises3_lpdf(b_zmu[206] | 0, log1p_exp(zkappa1));
  lprior += von_mises3_lpdf(b_zmu[207] | 0, log1p_exp(zkappa1));
  lprior += von_mises3_lpdf(b_zmu[208] | 0, log1p_exp(zkappa1));
  lprior += von_mises3_lpdf(b_zmu[209] | 0, log1p_exp(zkappa1));
  lprior += von_mises3_lpdf(b_zmu[210] | 0, log1p_exp(zkappa1));
  lprior += von_mises3_lpdf(b_zmu[211] | 0, log1p_exp(zkappa1));
  lprior += von_mises3_lpdf(b_zmu[212] | 0, log1p_exp(zkappa1));
  lprior += von_mises3_lpdf(b_zmu[213] | 0, log1p_exp(zkappa1));
  lprior += von_mises3_lpdf(b_zmu[214] | 0, log1p_exp(zkappa1));
  lprior += von_mises3_lpdf(b_zmu[215] | 0, log1p_exp(zkappa1));
  lprior += von_mises3_lpdf(b_zmu[216] | 0, log1p_exp(zkappa1));
  lprior += von_mises3_lpdf(b_zmu[217] | 0, log1p_exp(zkappa1));
  lprior += von_mises3_lpdf(b_zmu[218] | 0, log1p_exp(zkappa1));
  lprior += von_mises3_lpdf(b_zmu[219] | 0, log1p_exp(zkappa1));
  lprior += von_mises3_lpdf(b_zmu[220] | 0, log1p_exp(zkappa1));
  lprior += von_mises3_lpdf(b_zmu[221] | 0, log1p_exp(zkappa1));
  lprior += von_mises3_lpdf(b_zmu[222] | 0, log1p_exp(zkappa1));
  lprior += von_mises3_lpdf(b_zmu[223] | 0, log1p_exp(zkappa1));
  lprior += von_mises3_lpdf(b_zmu[224] | 0, log1p_exp(zkappa1+zkappa2));
  lprior += von_mises3_lpdf(b_zmu[225] | 0, log1p_exp(zkappa1));
  lprior += von_mises3_lpdf(b_zmu[226] | 0, log1p_exp(zkappa1));
  lprior += von_mises3_lpdf(b_zmu[227] | 0, log1p_exp(zkappa1));
  lprior += von_mises3_lpdf(b_zmu[228] | 0, log1p_exp(zkappa1));
  lprior += von_mises3_lpdf(b_zmu[229] | 0, log1p_exp(zkappa1));
  lprior += von_mises3_lpdf(b_zmu[230] | 0, log1p_exp(zkappa1));
  lprior += von_mises3_lpdf(b_zmu[231] | 0, log1p_exp(zkappa1));
  lprior += von_mises3_lpdf(b_zmu[232] | 0, log1p_exp(zkappa1));
  lprior += von_mises3_lpdf(b_zmu[233] | 0, log1p_exp(zkappa1));
  lprior += von_mises3_lpdf(b_zmu[234] | 0, log1p_exp(zkappa1));
  lprior += von_mises3_lpdf(b_zmu[235] | 0, log1p_exp(zkappa1));
  lprior += von_mises3_lpdf(b_zmu[236] | 0, log1p_exp(zkappa1));
  lprior += von_mises3_lpdf(b_zmu[237] | 0, log1p_exp(zkappa1));
  lprior += von_mises3_lpdf(b_zmu[238] | 0, log1p_exp(zkappa1));
  lprior += von_mises3_lpdf(b_zmu[239] | 0, log1p_exp(zkappa1));
  lprior += von_mises3_lpdf(b_zmu[240] | 0, log1p_exp(zkappa1+zkappa2));
  lprior += von_mises3_lpdf(b_zmu[241] | 0, log1p_exp(zkappa1));
  lprior += von_mises3_lpdf(b_zmu[242] | 0, log1p_exp(zkappa1+zkappa2));
  lprior += von_mises3_lpdf(b_zmu[243] | 0, log1p_exp(zkappa1));
  lprior += von_mises3_lpdf(b_zmu[244] | 0, log1p_exp(zkappa1));
  lprior += von_mises3_lpdf(b_zmu[245] | 0, log1p_exp(zkappa1));
  lprior += von_mises3_lpdf(b_zmu[246] | 0, log1p_exp(zkappa1));
  lprior += von_mises3_lpdf(b_zmu[247] | 0, log1p_exp(zkappa1));
  lprior += von_mises3_lpdf(b_zmu[248] | 0, log1p_exp(zkappa1));
  lprior += von_mises3_lpdf(b_zmu[249] | 0, log1p_exp(zkappa1));
  lprior += von_mises3_lpdf(b_zmu[250] | 0, log1p_exp(zkappa1));
  lprior += von_mises3_lpdf(b_zmu[251] | 0, log1p_exp(zkappa1));
  lprior += von_mises3_lpdf(b_zmu[252] | 0, log1p_exp(zkappa1));
  lprior += von_mises3_lpdf(b_zmu[253] | 0, log1p_exp(zkappa1+zkappa2));
  lprior += von_mises3_lpdf(b_zmu[254] | 0, log1p_exp(zkappa1));
  lprior += von_mises3_lpdf(b_zmu[255] | 0, log1p_exp(zkappa1));
  lprior += von_mises3_lpdf(b_zmu[256] | 0, log1p_exp(zkappa1));
  lprior += von_mises3_lpdf(b_zmu[257] | 0, log1p_exp(zkappa1));
  lprior += von_mises3_lpdf(b_zmu[258] | 0, log1p_exp(zkappa1));
  lprior += von_mises3_lpdf(b_zmu[259] | 0, log1p_exp(zkappa1));
  lprior += von_mises3_lpdf(b_zmu[260] | 0, log1p_exp(zkappa1));
  lprior += von_mises3_lpdf(b_zmu[261] | 0, log1p_exp(zkappa1));
  lprior += von_mises3_lpdf(b_zmu[262] | 0, log1p_exp(zkappa1));
  lprior += von_mises3_lpdf(b_zmu[263] | 0, log1p_exp(zkappa1));
  lprior += von_mises3_lpdf(b_zmu[264] | 0, log1p_exp(zkappa1));
  lprior += von_mises3_lpdf(b_zmu[265] | 0, log1p_exp(zkappa1));
  lprior += von_mises3_lpdf(b_zmu[266] | 0, log1p_exp(zkappa1));
  lprior += von_mises3_lpdf(b_zmu[267] | 0, log1p_exp(zkappa1));
  lprior += von_mises3_lpdf(b_zmu[268] | 0, log1p_exp(zkappa1));
  lprior += von_mises3_lpdf(b_zmu[269] | 0, log1p_exp(zkappa1));
  lprior += von_mises3_lpdf(b_zmu[270] | 0, log1p_exp(zkappa1));
  lprior += von_mises3_lpdf(b_zmu[271] | 0, log1p_exp(zkappa1+zkappa2));
  lprior += von_mises3_lpdf(b_zmu[272] | 0, log1p_exp(zkappa1));
  lprior += von_mises3_lpdf(b_zmu[273] | 0, log1p_exp(zkappa1));
  lprior += von_mises3_lpdf(b_zmu[274] | 0, log1p_exp(zkappa1));
  lprior += von_mises3_lpdf(b_zmu[275] | 0, log1p_exp(zkappa1));
  lprior += von_mises3_lpdf(b_zmu[276] | 0, log1p_exp(zkappa1));
  lprior += von_mises3_lpdf(b_zmu[277] | 0, log1p_exp(zkappa1));
  lprior += von_mises3_lpdf(b_zmu[278] | 0, log1p_exp(zkappa1));
  lprior += von_mises3_lpdf(b_zmu[279] | 0, log1p_exp(zkappa1));
  lprior += von_mises3_lpdf(b_zmu[280] | 0, log1p_exp(zkappa1));
  lprior += von_mises3_lpdf(b_zmu[281] | 0, log1p_exp(zkappa1));
  lprior += von_mises3_lpdf(b_zmu[282] | 0, log1p_exp(zkappa1));
  lprior += von_mises3_lpdf(b_zmu[283] | 0, log1p_exp(zkappa1));
  lprior += von_mises3_lpdf(b_zmu[284] | 0, log1p_exp(zkappa1));
  lprior += von_mises3_lpdf(b_zmu[285] | 0, log1p_exp(zkappa1));
  lprior += von_mises3_lpdf(b_zmu[286] | 0, log1p_exp(zkappa1));
  lprior += von_mises3_lpdf(b_zmu[287] | 0, log1p_exp(zkappa1));
  lprior += von_mises3_lpdf(b_zmu[288] | 0, log1p_exp(zkappa1));
  lprior += von_mises3_lpdf(b_zmu[289] | 0, log1p_exp(zkappa1));
  lprior += von_mises3_lpdf(b_zmu[290] | 0, log1p_exp(zkappa1));
  lprior += von_mises3_lpdf(b_zmu[291] | 0, log1p_exp(zkappa1));
  lprior += von_mises3_lpdf(b_zmu[292] | 0, log1p_exp(zkappa1));
  lprior += von_mises3_lpdf(b_zmu[293] | 0, log1p_exp(zkappa1));
  lprior += von_mises3_lpdf(b_zmu[294] | 0, log1p_exp(zkappa1));
  lprior += von_mises3_lpdf(b_zmu[295] | 0, log1p_exp(zkappa1));
  lprior += von_mises3_lpdf(b_zmu[296] | 0, log1p_exp(zkappa1));
  lprior += von_mises3_lpdf(b_zmu[297] | 0, log1p_exp(zkappa1));
  lprior += von_mises3_lpdf(b_zmu[298] | 0, log1p_exp(zkappa1));
  lprior += von_mises3_lpdf(b_zmu[299] | 0, log1p_exp(zkappa1));
  lprior += von_mises3_lpdf(b_zmu[300] | 0, log1p_exp(zkappa1));
  lprior += von_mises3_lpdf(b_zmu[301] | 0, log1p_exp(zkappa1));
  lprior += von_mises3_lpdf(b_zmu[302] | 0, log1p_exp(zkappa1));
  lprior += von_mises3_lpdf(b_zmu[303] | 0, log1p_exp(zkappa1));
  lprior += von_mises3_lpdf(b_zmu[304] | 0, log1p_exp(zkappa1));
  lprior += von_mises3_lpdf(b_zmu[305] | 0, log1p_exp(zkappa1));
  lprior += von_mises3_lpdf(b_zmu[306] | 0, log1p_exp(zkappa1));
  lprior += von_mises3_lpdf(b_zmu[307] | 0, log1p_exp(zkappa1));
  lprior += von_mises3_lpdf(b_zmu[308] | 0, log1p_exp(zkappa1));
  lprior += von_mises3_lpdf(b_zmu[309] | 0, log1p_exp(zkappa1));
  lprior += von_mises3_lpdf(b_zmu[310] | 0, log1p_exp(zkappa1));
  lprior += von_mises3_lpdf(b_zmu[311] | 0, log1p_exp(zkappa1));
  lprior += von_mises3_lpdf(b_zmu[312] | 0, log1p_exp(zkappa1));
  lprior += von_mises3_lpdf(b_zmu[313] | 0, log1p_exp(zkappa1));
  lprior += von_mises3_lpdf(b_zmu[314] | 0, log1p_exp(zkappa1));
  lprior += von_mises3_lpdf(b_zmu[315] | 0, log1p_exp(zkappa1));
  lprior += von_mises3_lpdf(b_zmu[316] | 0, log1p_exp(zkappa1));
  lprior += von_mises3_lpdf(b_zmu[317] | 0, log1p_exp(zkappa1+zkappa2));
  lprior += von_mises3_lpdf(b_zmu[318] | 0, log1p_exp(zkappa1));
  lprior += von_mises3_lpdf(b_zmu[319] | 0, log1p_exp(zkappa1));
  lprior += von_mises3_lpdf(b_zmu[320] | 0, log1p_exp(zkappa1));
  lprior += von_mises3_lpdf(b_zmu[321] | 0, log1p_exp(zkappa1));
  lprior += von_mises3_lpdf(b_zmu[322] | 0, log1p_exp(zkappa1));
  lprior += von_mises3_lpdf(b_zmu[323] | 0, log1p_exp(zkappa1));
  lprior += von_mises3_lpdf(b_zmu[324] | 0, log1p_exp(zkappa1));
  lprior += von_mises3_lpdf(b_zmu[325] | 0, log1p_exp(zkappa1));
  lprior += von_mises3_lpdf(b_zmu[326] | 0, log1p_exp(zkappa1+zkappa2));
  lprior += von_mises3_lpdf(b_zmu[327] | 0, log1p_exp(zkappa1));
  lprior += von_mises3_lpdf(b_zmu[328] | 0, log1p_exp(zkappa1));
  lprior += von_mises3_lpdf(b_zmu[329] | 0, log1p_exp(zkappa1));
  lprior += von_mises3_lpdf(b_zmu[330] | 0, log1p_exp(zkappa1));
  lprior += von_mises3_lpdf(b_zmu[331] | 0, log1p_exp(zkappa1));
  lprior += von_mises3_lpdf(b_zmu[332] | 0, log1p_exp(zkappa1));
  lprior += von_mises3_lpdf(b_zmu[333] | 0, log1p_exp(zkappa1));
  lprior += von_mises3_lpdf(b_zmu[334] | 0, log1p_exp(zkappa1));
  lprior += von_mises3_lpdf(b_zmu[335] | 0, log1p_exp(zkappa1));
  lprior += von_mises3_lpdf(b_zmu[336] | 0, log1p_exp(zkappa1));
  lprior += von_mises3_lpdf(b_zmu[337] | 0, log1p_exp(zkappa1));
  lprior += von_mises3_lpdf(b_zmu[338] | 0, log1p_exp(zkappa1));
  lprior += von_mises3_lpdf(b_zmu[339] | 0, log1p_exp(zkappa1));
  lprior += von_mises3_lpdf(b_zmu[340] | 0, log1p_exp(zkappa1));
  lprior += von_mises3_lpdf(b_zmu[341] | 0, log1p_exp(zkappa1));
  lprior += von_mises3_lpdf(b_zmu[342] | 0, log1p_exp(zkappa1));
  lprior += von_mises3_lpdf(b_zmu[343] | 0, log1p_exp(zkappa1));
  lprior += von_mises3_lpdf(b_zmu[344] | 0, log1p_exp(zkappa1));
  lprior += von_mises3_lpdf(b_zmu[345] | 0, log1p_exp(zkappa1));
  lprior += von_mises3_lpdf(b_zmu[346] | 0, log1p_exp(zkappa1));
  lprior += von_mises3_lpdf(b_zmu[347] | 0, log1p_exp(zkappa1));
  lprior += von_mises3_lpdf(b_zmu[348] | 0, log1p_exp(zkappa1));
  lprior += von_mises3_lpdf(b_zmu[349] | 0, log1p_exp(zkappa1));
  lprior += von_mises3_lpdf(b_zmu[350] | 0, log1p_exp(zkappa1));
  lprior += von_mises3_lpdf(b_zmu[351] | 0, log1p_exp(zkappa1));
  lprior += von_mises3_lpdf(b_zmu[352] | 0, log1p_exp(zkappa1));
  lprior += von_mises3_lpdf(b_zmu[353] | 0, log1p_exp(zkappa1));
  lprior += von_mises3_lpdf(b_zmu[354] | 0, log1p_exp(zkappa1));
  lprior += von_mises3_lpdf(b_zmu[355] | 0, log1p_exp(zkappa1));
  lprior += von_mises3_lpdf(b_zmu[356] | 0, log1p_exp(zkappa1));
  lprior += von_mises3_lpdf(b_zmu[357] | 0, log1p_exp(zkappa1+zkappa3));
  lprior += von_mises3_lpdf(b_zmu[358] | 0, log1p_exp(zkappa1));
  lprior += von_mises3_lpdf(b_zmu[359] | 0, log1p_exp(zkappa1));
  lprior += von_mises3_lpdf(b_zmu[360] | 0, log1p_exp(zkappa1+zkappa3));
  lprior += von_mises3_lpdf(b_zmu[361] | 0, log1p_exp(zkappa1));
  lprior += von_mises3_lpdf(b_zmu[362] | 0, log1p_exp(zkappa1+zkappa3));
  lprior += von_mises3_lpdf(b_zmu[363] | 0, log1p_exp(zkappa1));
  lprior += von_mises3_lpdf(b_zmu[364] | 0, log1p_exp(zkappa1));
  lprior += von_mises3_lpdf(b_zmu[365] | 0, log1p_exp(zkappa1));
  lprior += von_mises3_lpdf(b_zmu[366] | 0, log1p_exp(zkappa1));
  lprior += von_mises3_lpdf(b_zmu[367] | 0, log1p_exp(zkappa1));
  lprior += von_mises3_lpdf(b_zmu[368] | 0, log1p_exp(zkappa1));
  lprior += von_mises3_lpdf(b_zmu[369] | 0, log1p_exp(zkappa1));
  lprior += von_mises3_lpdf(b_zmu[370] | 0, log1p_exp(zkappa1));
  lprior += von_mises3_lpdf(b_zmu[371] | 0, log1p_exp(zkappa1));
  lprior += von_mises3_lpdf(b_zmu[372] | 0, log1p_exp(zkappa1));
  lprior += von_mises3_lpdf(b_zmu[373] | 0, log1p_exp(zkappa1));
  lprior += von_mises3_lpdf(b_zmu[374] | 0, log1p_exp(zkappa1));
  lprior += von_mises3_lpdf(b_zmu[375] | 0, log1p_exp(zkappa1));
  lprior += von_mises3_lpdf(b_zmu[376] | 0, log1p_exp(zkappa1));
  lprior += von_mises3_lpdf(b_zmu[377] | 0, log1p_exp(zkappa1));
  lprior += von_mises3_lpdf(b_zmu[378] | 0, log1p_exp(zkappa1));
  lprior += von_mises3_lpdf(b_zmu[379] | 0, log1p_exp(zkappa1));
  lprior += von_mises3_lpdf(b_zmu[380] | 0, log1p_exp(zkappa1));
  lprior += von_mises3_lpdf(b_zmu[381] | 0, log1p_exp(zkappa1));
  lprior += von_mises3_lpdf(b_zmu[382] | 0, log1p_exp(zkappa1));
  lprior += von_mises3_lpdf(b_zmu[383] | 0, log1p_exp(zkappa1));
  lprior += von_mises3_lpdf(b_zmu[384] | 0, log1p_exp(zkappa1));
  lprior += von_mises3_lpdf(b_zmu[385] | 0, log1p_exp(zkappa1));
  lprior += von_mises3_lpdf(b_zmu[386] | 0, log1p_exp(zkappa1));
  lprior += von_mises3_lpdf(b_zmu[387] | 0, log1p_exp(zkappa1));
  lprior += von_mises3_lpdf(b_zmu[388] | 0, log1p_exp(zkappa1));
  lprior += von_mises3_lpdf(b_zmu[389] | 0, log1p_exp(zkappa1));
  lprior += von_mises3_lpdf(b_zmu[390] | 0, log1p_exp(zkappa1));
  lprior += von_mises3_lpdf(b_zmu[391] | 0, log1p_exp(zkappa1));
  lprior += von_mises3_lpdf(b_zmu[392] | 0, log1p_exp(zkappa1));
  lprior += von_mises3_lpdf(b_zmu[393] | 0, log1p_exp(zkappa1+zkappa3));
  lprior += von_mises3_lpdf(b_zmu[394] | 0, log1p_exp(zkappa1));
  lprior += von_mises3_lpdf(b_zmu[395] | 0, log1p_exp(zkappa1));
  lprior += von_mises3_lpdf(b_zmu[396] | 0, log1p_exp(zkappa1));
  lprior += von_mises3_lpdf(b_zmu[397] | 0, log1p_exp(zkappa1));
  lprior += von_mises3_lpdf(b_zmu[398] | 0, log1p_exp(zkappa1));
  lprior += von_mises3_lpdf(b_zmu[399] | 0, log1p_exp(zkappa1));
  lprior += von_mises3_lpdf(b_zmu[400] | 0, log1p_exp(zkappa1));
  lprior += von_mises3_lpdf(b_zmu[401] | 0, log1p_exp(zkappa1));
  lprior += von_mises3_lpdf(b_zmu[402] | 0, log1p_exp(zkappa1));
  lprior += von_mises3_lpdf(b_zmu[403] | 0, log1p_exp(zkappa1));
  lprior += von_mises3_lpdf(b_zmu[404] | 0, log1p_exp(zkappa1));
  lprior += von_mises3_lpdf(b_zmu[405] | 0, log1p_exp(zkappa1));
  lprior += von_mises3_lpdf(b_zmu[406] | 0, log1p_exp(zkappa1));
  lprior += von_mises3_lpdf(b_zmu[407] | 0, log1p_exp(zkappa1));
  lprior += von_mises3_lpdf(b_zmu[408] | 0, log1p_exp(zkappa1));
  lprior += von_mises3_lpdf(b_zmu[409] | 0, log1p_exp(zkappa1+zkappa3));
  lprior += von_mises3_lpdf(b_zmu[410] | 0, log1p_exp(zkappa1));
  lprior += von_mises3_lpdf(b_zmu[411] | 0, log1p_exp(zkappa1+zkappa3));
  lprior += von_mises3_lpdf(b_zmu[412] | 0, log1p_exp(zkappa1));
  lprior += von_mises3_lpdf(b_zmu[413] | 0, log1p_exp(zkappa1));
  lprior += von_mises3_lpdf(b_zmu[414] | 0, log1p_exp(zkappa1));
  lprior += von_mises3_lpdf(b_zmu[415] | 0, log1p_exp(zkappa1));
  lprior += von_mises3_lpdf(b_zmu[416] | 0, log1p_exp(zkappa1));
  lprior += von_mises3_lpdf(b_zmu[417] | 0, log1p_exp(zkappa1));
  lprior += von_mises3_lpdf(b_zmu[418] | 0, log1p_exp(zkappa1));
  lprior += von_mises3_lpdf(b_zmu[419] | 0, log1p_exp(zkappa1));
  lprior += von_mises3_lpdf(b_zmu[420] | 0, log1p_exp(zkappa1));
  lprior += von_mises3_lpdf(b_zmu[421] | 0, log1p_exp(zkappa1));
  lprior += von_mises3_lpdf(b_zmu[422] | 0, log1p_exp(zkappa1+zkappa3));
  lprior += von_mises3_lpdf(b_zmu[423] | 0, log1p_exp(zkappa1));
  lprior += von_mises3_lpdf(b_zmu[424] | 0, log1p_exp(zkappa1));
  lprior += von_mises3_lpdf(b_zmu[425] | 0, log1p_exp(zkappa1));
  lprior += von_mises3_lpdf(b_zmu[426] | 0, log1p_exp(zkappa1));
  lprior += von_mises3_lpdf(b_zmu[427] | 0, log1p_exp(zkappa1));
  lprior += von_mises3_lpdf(b_zmu[428] | 0, log1p_exp(zkappa1));
  lprior += von_mises3_lpdf(b_zmu[429] | 0, log1p_exp(zkappa1));
  lprior += von_mises3_lpdf(b_zmu[430] | 0, log1p_exp(zkappa1));
  lprior += von_mises3_lpdf(b_zmu[431] | 0, log1p_exp(zkappa1));
  lprior += von_mises3_lpdf(b_zmu[432] | 0, log1p_exp(zkappa1));
  lprior += von_mises3_lpdf(b_zmu[433] | 0, log1p_exp(zkappa1));
  lprior += von_mises3_lpdf(b_zmu[434] | 0, log1p_exp(zkappa1));
  lprior += von_mises3_lpdf(b_zmu[435] | 0, log1p_exp(zkappa1));
  lprior += von_mises3_lpdf(b_zmu[436] | 0, log1p_exp(zkappa1));
  lprior += von_mises3_lpdf(b_zmu[437] | 0, log1p_exp(zkappa1));
  lprior += von_mises3_lpdf(b_zmu[438] | 0, log1p_exp(zkappa1));
  lprior += von_mises3_lpdf(b_zmu[439] | 0, log1p_exp(zkappa1));
  lprior += von_mises3_lpdf(b_zmu[440] | 0, log1p_exp(zkappa1+zkappa3));
  lprior += von_mises3_lpdf(b_zmu[441] | 0, log1p_exp(zkappa1));
  lprior += von_mises3_lpdf(b_zmu[442] | 0, log1p_exp(zkappa1));
  lprior += von_mises3_lpdf(b_zmu[443] | 0, log1p_exp(zkappa1));
  lprior += von_mises3_lpdf(b_zmu[444] | 0, log1p_exp(zkappa1));
  lprior += von_mises3_lpdf(b_zmu[445] | 0, log1p_exp(zkappa1));
  lprior += von_mises3_lpdf(b_zmu[446] | 0, log1p_exp(zkappa1));
  lprior += von_mises3_lpdf(b_zmu[447] | 0, log1p_exp(zkappa1));
  lprior += von_mises3_lpdf(b_zmu[448] | 0, log1p_exp(zkappa1));
  lprior += von_mises3_lpdf(b_zmu[449] | 0, log1p_exp(zkappa1));
  lprior += von_mises3_lpdf(b_zmu[450] | 0, log1p_exp(zkappa1));
  lprior += von_mises3_lpdf(b_zmu[451] | 0, log1p_exp(zkappa1));
  lprior += von_mises3_lpdf(b_zmu[452] | 0, log1p_exp(zkappa1));
  lprior += von_mises3_lpdf(b_zmu[453] | 0, log1p_exp(zkappa1));
  lprior += von_mises3_lpdf(b_zmu[454] | 0, log1p_exp(zkappa1));
  lprior += von_mises3_lpdf(b_zmu[455] | 0, log1p_exp(zkappa1));
  lprior += von_mises3_lpdf(b_zmu[456] | 0, log1p_exp(zkappa1));
  lprior += von_mises3_lpdf(b_zmu[457] | 0, log1p_exp(zkappa1));
  lprior += von_mises3_lpdf(b_zmu[458] | 0, log1p_exp(zkappa1));
  lprior += von_mises3_lpdf(b_zmu[459] | 0, log1p_exp(zkappa1));
  lprior += von_mises3_lpdf(b_zmu[460] | 0, log1p_exp(zkappa1));
  lprior += von_mises3_lpdf(b_zmu[461] | 0, log1p_exp(zkappa1));
  lprior += von_mises3_lpdf(b_zmu[462] | 0, log1p_exp(zkappa1));
  lprior += von_mises3_lpdf(b_zmu[463] | 0, log1p_exp(zkappa1));
  lprior += von_mises3_lpdf(b_zmu[464] | 0, log1p_exp(zkappa1));
  lprior += von_mises3_lpdf(b_zmu[465] | 0, log1p_exp(zkappa1));
  lprior += von_mises3_lpdf(b_zmu[466] | 0, log1p_exp(zkappa1));
  lprior += von_mises3_lpdf(b_zmu[467] | 0, log1p_exp(zkappa1));
  lprior += von_mises3_lpdf(b_zmu[468] | 0, log1p_exp(zkappa1));
  lprior += von_mises3_lpdf(b_zmu[469] | 0, log1p_exp(zkappa1));
  lprior += von_mises3_lpdf(b_zmu[470] | 0, log1p_exp(zkappa1));
  lprior += von_mises3_lpdf(b_zmu[471] | 0, log1p_exp(zkappa1));
  lprior += von_mises3_lpdf(b_zmu[472] | 0, log1p_exp(zkappa1));
  lprior += von_mises3_lpdf(b_zmu[473] | 0, log1p_exp(zkappa1));
  lprior += von_mises3_lpdf(b_zmu[474] | 0, log1p_exp(zkappa1));
  lprior += von_mises3_lpdf(b_zmu[475] | 0, log1p_exp(zkappa1));
  lprior += von_mises3_lpdf(b_zmu[476] | 0, log1p_exp(zkappa1));
  lprior += von_mises3_lpdf(b_zmu[477] | 0, log1p_exp(zkappa1));
  lprior += von_mises3_lpdf(b_zmu[478] | 0, log1p_exp(zkappa1));
  lprior += von_mises3_lpdf(b_zmu[479] | 0, log1p_exp(zkappa1));
  lprior += von_mises3_lpdf(b_zmu[480] | 0, log1p_exp(zkappa1));
  lprior += von_mises3_lpdf(b_zmu[481] | 0, log1p_exp(zkappa1));
  lprior += von_mises3_lpdf(b_zmu[482] | 0, log1p_exp(zkappa1));
  lprior += von_mises3_lpdf(b_zmu[483] | 0, log1p_exp(zkappa1));
  lprior += von_mises3_lpdf(b_zmu[484] | 0, log1p_exp(zkappa1));
  lprior += von_mises3_lpdf(b_zmu[485] | 0, log1p_exp(zkappa1));
  lprior += von_mises3_lpdf(b_zmu[486] | 0, log1p_exp(zkappa1+zkappa3));
  lprior += von_mises3_lpdf(b_zmu[487] | 0, log1p_exp(zkappa1));
  lprior += von_mises3_lpdf(b_zmu[488] | 0, log1p_exp(zkappa1));
  lprior += von_mises3_lpdf(b_zmu[489] | 0, log1p_exp(zkappa1));
  lprior += von_mises3_lpdf(b_zmu[490] | 0, log1p_exp(zkappa1));
  lprior += von_mises3_lpdf(b_zmu[491] | 0, log1p_exp(zkappa1));
  lprior += von_mises3_lpdf(b_zmu[492] | 0, log1p_exp(zkappa1));
  lprior += von_mises3_lpdf(b_zmu[493] | 0, log1p_exp(zkappa1));
  lprior += von_mises3_lpdf(b_zmu[494] | 0, log1p_exp(zkappa1));
  lprior += von_mises3_lpdf(b_zmu[495] | 0, log1p_exp(zkappa1+zkappa3));
  lprior += von_mises3_lpdf(b_zmu[496] | 0, log1p_exp(zkappa1));
  lprior += von_mises3_lpdf(b_zmu[497] | 0, log1p_exp(zkappa1));
  lprior += von_mises3_lpdf(b_zmu[498] | 0, log1p_exp(zkappa1));
  lprior += von_mises3_lpdf(b_zmu[499] | 0, log1p_exp(zkappa1));
  lprior += von_mises3_lpdf(b_zmu[500] | 0, log1p_exp(zkappa1));
  lprior += von_mises3_lpdf(b_zmu[501] | 0, log1p_exp(zkappa1));
  lprior += von_mises3_lpdf(b_zmu[502] | 0, log1p_exp(zkappa1));
  lprior += von_mises3_lpdf(b_zmu[503] | 0, log1p_exp(zkappa1));
  lprior += von_mises3_lpdf(b_zmu[504] | 0, log1p_exp(zkappa1));
  lprior += von_mises3_lpdf(b_zmu[505] | 0, log1p_exp(zkappa1));
  lprior += von_mises3_lpdf(b_zmu[506] | 0, log1p_exp(zkappa1));
  lprior += von_mises3_lpdf(b_zmu[507] | 0, log1p_exp(zkappa1));
  lprior += von_mises3_lpdf(b_zmu[508] | 0, log1p_exp(zkappa1));
  lprior += von_mises3_lpdf(b_zmu[509] | 0, log1p_exp(zkappa1));
  lprior += von_mises3_lpdf(b_zmu[510] | 0, log1p_exp(zkappa1));
  lprior += von_mises3_lpdf(b_zmu[511] | 0, log1p_exp(zkappa1));
  lprior += von_mises3_lpdf(b_zmu[512] | 0, log1p_exp(zkappa1));
  lprior += von_mises3_lpdf(b_zmu[513] | 0, log1p_exp(zkappa1));
  lprior += von_mises3_lpdf(b_zmu[514] | 0, log1p_exp(zkappa1));
  lprior += von_mises3_lpdf(b_zmu[515] | 0, log1p_exp(zkappa1));
  lprior += von_mises3_lpdf(b_zmu[516] | 0, log1p_exp(zkappa1));
  lprior += von_mises3_lpdf(b_zmu[517] | 0, log1p_exp(zkappa1));
  lprior += von_mises3_lpdf(b_zmu[518] | 0, log1p_exp(zkappa1));
  lprior += von_mises3_lpdf(b_zmu[519] | 0, log1p_exp(zkappa1));
  lprior += von_mises3_lpdf(b_zmu[520] | 0, log1p_exp(zkappa1));
  lprior += von_mises3_lpdf(b_zmu[521] | 0, log1p_exp(zkappa1));
  lprior += von_mises3_lpdf(b_zmu[522] | 0, log1p_exp(zkappa1));
  lprior += von_mises3_lpdf(b_zmu[523] | 0, log1p_exp(zkappa1));
  lprior += von_mises3_lpdf(b_zmu[524] | 0, log1p_exp(zkappa1));
  lprior += von_mises3_lpdf(b_zmu[525] | 0, log1p_exp(zkappa1));
  lprior += von_mises3_lpdf(b_zmu[526] | 0, log1p_exp(zkappa1+zkappa4));
  lprior += von_mises3_lpdf(b_zmu[527] | 0, log1p_exp(zkappa1));
  lprior += von_mises3_lpdf(b_zmu[528] | 0, log1p_exp(zkappa1));
  lprior += von_mises3_lpdf(b_zmu[529] | 0, log1p_exp(zkappa1+zkappa4));
  lprior += von_mises3_lpdf(b_zmu[530] | 0, log1p_exp(zkappa1));
  lprior += von_mises3_lpdf(b_zmu[531] | 0, log1p_exp(zkappa1+zkappa4));
  lprior += von_mises3_lpdf(b_zmu[532] | 0, log1p_exp(zkappa1));
  lprior += von_mises3_lpdf(b_zmu[533] | 0, log1p_exp(zkappa1));
  lprior += von_mises3_lpdf(b_zmu[534] | 0, log1p_exp(zkappa1));
  lprior += von_mises3_lpdf(b_zmu[535] | 0, log1p_exp(zkappa1));
  lprior += von_mises3_lpdf(b_zmu[536] | 0, log1p_exp(zkappa1));
  lprior += von_mises3_lpdf(b_zmu[537] | 0, log1p_exp(zkappa1));
  lprior += von_mises3_lpdf(b_zmu[538] | 0, log1p_exp(zkappa1));
  lprior += von_mises3_lpdf(b_zmu[539] | 0, log1p_exp(zkappa1));
  lprior += von_mises3_lpdf(b_zmu[540] | 0, log1p_exp(zkappa1));
  lprior += von_mises3_lpdf(b_zmu[541] | 0, log1p_exp(zkappa1));
  lprior += von_mises3_lpdf(b_zmu[542] | 0, log1p_exp(zkappa1));
  lprior += von_mises3_lpdf(b_zmu[543] | 0, log1p_exp(zkappa1));
  lprior += von_mises3_lpdf(b_zmu[544] | 0, log1p_exp(zkappa1));
  lprior += von_mises3_lpdf(b_zmu[545] | 0, log1p_exp(zkappa1));
  lprior += von_mises3_lpdf(b_zmu[546] | 0, log1p_exp(zkappa1));
  lprior += von_mises3_lpdf(b_zmu[547] | 0, log1p_exp(zkappa1));
  lprior += von_mises3_lpdf(b_zmu[548] | 0, log1p_exp(zkappa1));
  lprior += von_mises3_lpdf(b_zmu[549] | 0, log1p_exp(zkappa1));
  lprior += von_mises3_lpdf(b_zmu[550] | 0, log1p_exp(zkappa1));
  lprior += von_mises3_lpdf(b_zmu[551] | 0, log1p_exp(zkappa1));
  lprior += von_mises3_lpdf(b_zmu[552] | 0, log1p_exp(zkappa1));
  lprior += von_mises3_lpdf(b_zmu[553] | 0, log1p_exp(zkappa1));
  lprior += von_mises3_lpdf(b_zmu[554] | 0, log1p_exp(zkappa1));
  lprior += von_mises3_lpdf(b_zmu[555] | 0, log1p_exp(zkappa1));
  lprior += von_mises3_lpdf(b_zmu[556] | 0, log1p_exp(zkappa1));
  lprior += von_mises3_lpdf(b_zmu[557] | 0, log1p_exp(zkappa1));
  lprior += von_mises3_lpdf(b_zmu[558] | 0, log1p_exp(zkappa1));
  lprior += von_mises3_lpdf(b_zmu[559] | 0, log1p_exp(zkappa1));
  lprior += von_mises3_lpdf(b_zmu[560] | 0, log1p_exp(zkappa1));
  lprior += von_mises3_lpdf(b_zmu[561] | 0, log1p_exp(zkappa1));
  lprior += von_mises3_lpdf(b_zmu[562] | 0, log1p_exp(zkappa1+zkappa4));
  lprior += von_mises3_lpdf(b_zmu[563] | 0, log1p_exp(zkappa1));
  lprior += von_mises3_lpdf(b_zmu[564] | 0, log1p_exp(zkappa1));
  lprior += von_mises3_lpdf(b_zmu[565] | 0, log1p_exp(zkappa1));
  lprior += von_mises3_lpdf(b_zmu[566] | 0, log1p_exp(zkappa1));
  lprior += von_mises3_lpdf(b_zmu[567] | 0, log1p_exp(zkappa1));
  lprior += von_mises3_lpdf(b_zmu[568] | 0, log1p_exp(zkappa1));
  lprior += von_mises3_lpdf(b_zmu[569] | 0, log1p_exp(zkappa1));
  lprior += von_mises3_lpdf(b_zmu[570] | 0, log1p_exp(zkappa1));
  lprior += von_mises3_lpdf(b_zmu[571] | 0, log1p_exp(zkappa1));
  lprior += von_mises3_lpdf(b_zmu[572] | 0, log1p_exp(zkappa1));
  lprior += von_mises3_lpdf(b_zmu[573] | 0, log1p_exp(zkappa1));
  lprior += von_mises3_lpdf(b_zmu[574] | 0, log1p_exp(zkappa1));
  lprior += von_mises3_lpdf(b_zmu[575] | 0, log1p_exp(zkappa1));
  lprior += von_mises3_lpdf(b_zmu[576] | 0, log1p_exp(zkappa1));
  lprior += von_mises3_lpdf(b_zmu[577] | 0, log1p_exp(zkappa1));
  lprior += von_mises3_lpdf(b_zmu[578] | 0, log1p_exp(zkappa1+zkappa4));
  lprior += von_mises3_lpdf(b_zmu[579] | 0, log1p_exp(zkappa1));
  lprior += von_mises3_lpdf(b_zmu[580] | 0, log1p_exp(zkappa1+zkappa4));
  lprior += von_mises3_lpdf(b_zmu[581] | 0, log1p_exp(zkappa1));
  lprior += von_mises3_lpdf(b_zmu[582] | 0, log1p_exp(zkappa1));
  lprior += von_mises3_lpdf(b_zmu[583] | 0, log1p_exp(zkappa1));
  lprior += von_mises3_lpdf(b_zmu[584] | 0, log1p_exp(zkappa1));
  lprior += von_mises3_lpdf(b_zmu[585] | 0, log1p_exp(zkappa1));
  lprior += von_mises3_lpdf(b_zmu[586] | 0, log1p_exp(zkappa1));
  lprior += von_mises3_lpdf(b_zmu[587] | 0, log1p_exp(zkappa1));
  lprior += von_mises3_lpdf(b_zmu[588] | 0, log1p_exp(zkappa1));
  lprior += von_mises3_lpdf(b_zmu[589] | 0, log1p_exp(zkappa1));
  lprior += von_mises3_lpdf(b_zmu[590] | 0, log1p_exp(zkappa1));
  lprior += von_mises3_lpdf(b_zmu[591] | 0, log1p_exp(zkappa1+zkappa4));
  lprior += von_mises3_lpdf(b_zmu[592] | 0, log1p_exp(zkappa1));
  lprior += von_mises3_lpdf(b_zmu[593] | 0, log1p_exp(zkappa1));
  lprior += von_mises3_lpdf(b_zmu[594] | 0, log1p_exp(zkappa1));
  lprior += von_mises3_lpdf(b_zmu[595] | 0, log1p_exp(zkappa1));
  lprior += von_mises3_lpdf(b_zmu[596] | 0, log1p_exp(zkappa1));
  lprior += von_mises3_lpdf(b_zmu[597] | 0, log1p_exp(zkappa1));
  lprior += von_mises3_lpdf(b_zmu[598] | 0, log1p_exp(zkappa1));
  lprior += von_mises3_lpdf(b_zmu[599] | 0, log1p_exp(zkappa1));
  lprior += von_mises3_lpdf(b_zmu[600] | 0, log1p_exp(zkappa1));
  lprior += von_mises3_lpdf(b_zmu[601] | 0, log1p_exp(zkappa1));
  lprior += von_mises3_lpdf(b_zmu[602] | 0, log1p_exp(zkappa1));
  lprior += von_mises3_lpdf(b_zmu[603] | 0, log1p_exp(zkappa1));
  lprior += von_mises3_lpdf(b_zmu[604] | 0, log1p_exp(zkappa1));
  lprior += von_mises3_lpdf(b_zmu[605] | 0, log1p_exp(zkappa1));
  lprior += von_mises3_lpdf(b_zmu[606] | 0, log1p_exp(zkappa1));
  lprior += von_mises3_lpdf(b_zmu[607] | 0, log1p_exp(zkappa1));
  lprior += von_mises3_lpdf(b_zmu[608] | 0, log1p_exp(zkappa1));
  lprior += von_mises3_lpdf(b_zmu[609] | 0, log1p_exp(zkappa1+zkappa4));
  lprior += von_mises3_lpdf(b_zmu[610] | 0, log1p_exp(zkappa1));
  lprior += von_mises3_lpdf(b_zmu[611] | 0, log1p_exp(zkappa1));
  lprior += von_mises3_lpdf(b_zmu[612] | 0, log1p_exp(zkappa1));
  lprior += von_mises3_lpdf(b_zmu[613] | 0, log1p_exp(zkappa1));
  lprior += von_mises3_lpdf(b_zmu[614] | 0, log1p_exp(zkappa1));
  lprior += von_mises3_lpdf(b_zmu[615] | 0, log1p_exp(zkappa1));
  lprior += von_mises3_lpdf(b_zmu[616] | 0, log1p_exp(zkappa1));
  lprior += von_mises3_lpdf(b_zmu[617] | 0, log1p_exp(zkappa1));
  lprior += von_mises3_lpdf(b_zmu[618] | 0, log1p_exp(zkappa1));
  lprior += von_mises3_lpdf(b_zmu[619] | 0, log1p_exp(zkappa1));
  lprior += von_mises3_lpdf(b_zmu[620] | 0, log1p_exp(zkappa1));
  lprior += von_mises3_lpdf(b_zmu[621] | 0, log1p_exp(zkappa1));
  lprior += von_mises3_lpdf(b_zmu[622] | 0, log1p_exp(zkappa1));
  lprior += von_mises3_lpdf(b_zmu[623] | 0, log1p_exp(zkappa1));
  lprior += von_mises3_lpdf(b_zmu[624] | 0, log1p_exp(zkappa1));
  lprior += von_mises3_lpdf(b_zmu[625] | 0, log1p_exp(zkappa1));
  lprior += von_mises3_lpdf(b_zmu[626] | 0, log1p_exp(zkappa1));
  lprior += von_mises3_lpdf(b_zmu[627] | 0, log1p_exp(zkappa1));
  lprior += von_mises3_lpdf(b_zmu[628] | 0, log1p_exp(zkappa1));
  lprior += von_mises3_lpdf(b_zmu[629] | 0, log1p_exp(zkappa1));
  lprior += von_mises3_lpdf(b_zmu[630] | 0, log1p_exp(zkappa1));
  lprior += von_mises3_lpdf(b_zmu[631] | 0, log1p_exp(zkappa1));
  lprior += von_mises3_lpdf(b_zmu[632] | 0, log1p_exp(zkappa1));
  lprior += von_mises3_lpdf(b_zmu[633] | 0, log1p_exp(zkappa1));
  lprior += von_mises3_lpdf(b_zmu[634] | 0, log1p_exp(zkappa1));
  lprior += von_mises3_lpdf(b_zmu[635] | 0, log1p_exp(zkappa1));
  lprior += von_mises3_lpdf(b_zmu[636] | 0, log1p_exp(zkappa1));
  lprior += von_mises3_lpdf(b_zmu[637] | 0, log1p_exp(zkappa1));
  lprior += von_mises3_lpdf(b_zmu[638] | 0, log1p_exp(zkappa1));
  lprior += von_mises3_lpdf(b_zmu[639] | 0, log1p_exp(zkappa1));
  lprior += von_mises3_lpdf(b_zmu[640] | 0, log1p_exp(zkappa1));
  lprior += von_mises3_lpdf(b_zmu[641] | 0, log1p_exp(zkappa1));
  lprior += von_mises3_lpdf(b_zmu[642] | 0, log1p_exp(zkappa1));
  lprior += von_mises3_lpdf(b_zmu[643] | 0, log1p_exp(zkappa1));
  lprior += von_mises3_lpdf(b_zmu[644] | 0, log1p_exp(zkappa1));
  lprior += von_mises3_lpdf(b_zmu[645] | 0, log1p_exp(zkappa1));
  lprior += von_mises3_lpdf(b_zmu[646] | 0, log1p_exp(zkappa1));
  lprior += von_mises3_lpdf(b_zmu[647] | 0, log1p_exp(zkappa1));
  lprior += von_mises3_lpdf(b_zmu[648] | 0, log1p_exp(zkappa1));
  lprior += von_mises3_lpdf(b_zmu[649] | 0, log1p_exp(zkappa1));
  lprior += von_mises3_lpdf(b_zmu[650] | 0, log1p_exp(zkappa1));
  lprior += von_mises3_lpdf(b_zmu[651] | 0, log1p_exp(zkappa1));
  lprior += von_mises3_lpdf(b_zmu[652] | 0, log1p_exp(zkappa1));
  lprior += von_mises3_lpdf(b_zmu[653] | 0, log1p_exp(zkappa1));
  lprior += von_mises3_lpdf(b_zmu[654] | 0, log1p_exp(zkappa1));
  lprior += von_mises3_lpdf(b_zmu[655] | 0, log1p_exp(zkappa1+zkappa4));
  lprior += von_mises3_lpdf(b_zmu[656] | 0, log1p_exp(zkappa1));
  lprior += von_mises3_lpdf(b_zmu[657] | 0, log1p_exp(zkappa1));
  lprior += von_mises3_lpdf(b_zmu[658] | 0, log1p_exp(zkappa1));
  lprior += von_mises3_lpdf(b_zmu[659] | 0, log1p_exp(zkappa1));
  lprior += von_mises3_lpdf(b_zmu[660] | 0, log1p_exp(zkappa1));
  lprior += von_mises3_lpdf(b_zmu[661] | 0, log1p_exp(zkappa1));
  lprior += von_mises3_lpdf(b_zmu[662] | 0, log1p_exp(zkappa1));
  lprior += von_mises3_lpdf(b_zmu[663] | 0, log1p_exp(zkappa1));
  lprior += von_mises3_lpdf(b_zmu[664] | 0, log1p_exp(zkappa1+zkappa4));
  lprior += von_mises3_lpdf(b_zmu[665] | 0, log1p_exp(zkappa1));
  lprior += von_mises3_lpdf(b_zmu[666] | 0, log1p_exp(zkappa1));
  lprior += von_mises3_lpdf(b_zmu[667] | 0, log1p_exp(zkappa1));
  lprior += von_mises3_lpdf(b_zmu[668] | 0, log1p_exp(zkappa1));
  lprior += von_mises3_lpdf(b_zmu[669] | 0, log1p_exp(zkappa1));
  lprior += von_mises3_lpdf(b_zmu[670] | 0, log1p_exp(zkappa1));
  lprior += von_mises3_lpdf(b_zmu[671] | 0, log1p_exp(zkappa1));
  lprior += von_mises3_lpdf(b_zmu[672] | 0, log1p_exp(zkappa1));
  lprior += von_mises3_lpdf(b_zmu[673] | 0, log1p_exp(zkappa1));
  lprior += von_mises3_lpdf(b_zmu[674] | 0, log1p_exp(zkappa1));
  lprior += von_mises3_lpdf(b_zmu[675] | 0, log1p_exp(zkappa1));
  lprior += von_mises3_lpdf(b_zmu[676] | 0, log1p_exp(zkappa1));
  lprior += normal_lpdf(b_kappa[1] | 0.0, 5.0);
  lprior += normal_lpdf(b_kappa[2] | 0.0, 5.0);
  lprior += normal_lpdf(b_kappa[3] | 0.0, 5.0);
  lprior += normal_lpdf(Intercept_kappa | 2.0, 5.0);
  lprior += student_t_lpdf(sd_1[1] | 3, 0, 5.0)
    - 1 * student_t_lccdf(0 | 3, 0, 5.0);
  lprior += student_t_lpdf(sd_1[2] | 3, 0, 5.0)
    - 1 * student_t_lccdf(0 | 3, 0, 5.0);
  lprior += student_t_lpdf(sd_1[3] | 3, 0, 5.0)
    - 1 * student_t_lccdf(0 | 3, 0, 5.0);
  lprior += student_t_lpdf(sd_1[4] | 3, 0, 5.0)
    - 1 * student_t_lccdf(0 | 3, 0, 5.0);
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
    for (n in 1:N) {
      target += von_mises2_lpdf(Y[n] | mu[n], kappa[n]);
    }
  }
  // priors including constants
  target += lprior;
  target += std_normal_lpdf(to_vector(z_1));
  target += student_t_lpdf(zkappa1 | 3, 25, 5);
  target += student_t_lpdf(zkappa1+zkappa2 | 3, 25, 5);
  target += student_t_lpdf(zkappa1+zkappa3 | 3, 25, 5);
  target += student_t_lpdf(zkappa1+zkappa4 | 3, 25, 5);
}
generated quantities {
  // actual population-level intercept
  real b_kappa_Intercept = Intercept_kappa - dot_product(means_X_kappa, b_kappa);
  // compute group-level correlations
  corr_matrix[M_1] Cor_1 = multiply_lower_tri_self_transpose(L_1);
  vector<lower=-1,upper=1>[NC_1] cor_1;
  // extract upper diagonal of correlation matrix
  for (k in 1:M_1) {
    for (j in 1:(k - 1)) {
      cor_1[choose(k - 1, 2) + j] = Cor_1[j, k];
    }
  }
}

