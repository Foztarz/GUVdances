// generated with brms 2.21.0
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
   //circular functions
     // convert angle to modulo [-pi, pi]
   real mod_circular(real y){
     return(atan2(  sin(y), cos(y) ) ); // the output will always be in [-pi,pi]
   }   
   // convert array of angles to modulo [-pi, pi]
   vector mod_circular_array(vector y, int N){
     array[N] real yy ;
     for (n in 1:N) {
       yy[n] = atan2(  sin(y[n]), cos(y[n]) ) ; 
     }
     return(y ); // the output will always be in [-pi,pi]
   }
   // calculate the mean angle of a circular distribution (in radians)
   real mean_circular(vector y, int N){
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
    // calculate the mean array length for a circular distribution
   real rho_circular(vector y, int N){
      real sumsin = 0;
      real sumcos = 0;
     for (n in 1:N){
       sumsin += sin(y[n]);
       sumcos += cos(y[n]);
     }
     sumsin = sumsin/N;
     sumcos = sumcos/N;
     return(sqrt(  sumsin^2 + sumcos^2)/N );
   }
}
data {
  int<lower=1> N;  // total number of observations
  vector[N] Y;  // response variable
  int<lower=1> K;  // number of population-level effects
  matrix[N, K] X;  // population-level design matrix
  int<lower=1> Kc;  // number of population-level effects after centering
  int<lower=1> K_kappa;  // number of population-level effects
  matrix[N, K_kappa] X_kappa;  // population-level design matrix
  int<lower=1> Kc_kappa;  // number of population-level effects after centering
  // data for group-level effects of ID 1
  int<lower=1> N_1;  // number of grouping levels
  int<lower=1> M_1;  // number of coefficients per level
  array[N] int<lower=1> J_1;  // grouping indicator per observation
  // group-level predictor values
  vector[N] Z_1_1;
  // data for group-level effects of ID 2
  int<lower=1> N_2;  // number of grouping levels
  int<lower=1> M_2;  // number of coefficients per level
  array[N] int<lower=1> J_2;  // grouping indicator per observation
  // group-level predictor values
  vector[N] Z_2_kappa_1;
  int prior_only;  // should the likelihood be ignored?
}
transformed data {
  matrix[N, Kc] Xc;  // centered version of X without an intercept
  vector[Kc] means_X;  // column means of X before centering
  matrix[N, Kc_kappa] Xc_kappa;  // centered version of X_kappa without an intercept
  vector[Kc_kappa] means_X_kappa;  // column means of X_kappa before centering
  for (i in 2:K) {
    // means_X[i - 1] = mean(X[, i]);
     means_X[i - 1] = mean_circular(X[, i], num_elements(X[, i])); //this should be the circular mean
    Xc[, i - 1] = X[, i] - means_X[i - 1];
  }
  for (i in 2:K_kappa) {
    means_X_kappa[i - 1] = mean(X_kappa[, i]);
    Xc_kappa[, i - 1] = X_kappa[, i] - means_X_kappa[i - 1];
  }
}
parameters {
  vector[Kc] b;  // regression coefficients
  real Intercept;  // temporary intercept for centered predictors
  vector[Kc_kappa] b_kappa;  // regression coefficients
  real Intercept_kappa;  // temporary intercept for centered predictors
  //not needed for circular?
  // vector<lower=0>[M_1] sd_1;  // group-level standard deviations
  array[M_1] vector[N_1] z_1;  // standardized group-level effects
  vector<lower=0>[M_2] sd_2;  // group-level standard deviations
  array[M_2] vector[N_2] z_2;  // standardized group-level effects
}
transformed parameters {
  vector[N_1] r_1_1;  // actual group-level effects
  vector[N_2] r_2_kappa_1;  // actual group-level effects
  real lprior = 0;  // prior contributions to the log posterior
  // real cmean = mean_circular(z_1[1], N_1); //calculate the mean of all individuals
  real cmean = mean_circular(Y, N); //calculate the mean of all individuals
  real crho = rho_circular(Y, N); //calculate the mean of all individuals
  // r_1_1 = (sd_1[1] * (z_1[1]));
  r_1_1 =  mod_circular_array(z_1[1], N_1);
  r_2_kappa_1 = (sd_2[1] * (z_2[1]));
  // lprior += student_t_lpdf(Intercept | 3, 0, 2.5);//prior for intercept at 171°?!
  // lprior += student_t_lpdf(mod_circular(Intercept) | 3, cmean, 2.5);//prior for intercept at mean of raneff
  lprior += von_mises_lpdf(mod_circular(Intercept) | cmean, 2*crho+crho^3+(5*crho^5)/6);//prior for intercept at mean of raneff
  // lprior += normal_lpdf(Intercept_kappa | 5.0, 0.8);//this seems like a pathological prior!
  lprior += normal_lpdf(Intercept_kappa | 0, 0.35);//expect kappa in (0.5, 2) 
  // lprior += student_t_lpdf(sd_1 | 3, 0, 2.5)
  //   - 1 * student_t_lccdf(0 | 3, 0, 2.5);
  // lprior += student_t_lpdf(sd_2 | 3, 0, 2.5)
  //   - 1 * student_t_lccdf(0 | 3, 0, 2.5);
  lprior += student_t_lpdf(sd_2 | 3, 0, 0.35)
    - 1 * student_t_lccdf(0 | 3, 0, 0.35); //sigma same scale as for other kappa prior
}
model {
  // likelihood including constants
  if (!prior_only) {
    // initialize linear predictor term
    vector[N] mu = rep_vector(0.0, N);
    // initialize linear predictor term
    vector[N] kappa = rep_vector(0.0, N);
    mu += Intercept + Xc * b;
    kappa += Intercept_kappa + Xc_kappa * b_kappa;
    for (n in 1:N) {
      // add more terms to the linear predictor
      mu[n] += r_1_1[J_1[n]] * Z_1_1[n];
    }
    for (n in 1:N) {
      // add more terms to the linear predictor
      kappa[n] += r_2_kappa_1[J_2[n]] * Z_2_kappa_1[n];
    }
    // mu = inv_tan_half(mu);
    mu = mod_circular_array(mu, N);
    kappa = exp(kappa);
    for (n in 1:N) {
      target += von_mises2_lpdf(Y[n] | mu[n], kappa[n]);
    }
  }
  // priors including constants
  target += lprior;
  // target += std_normal_lpdf(z_1[1]); //no prior
  target += std_normal_lpdf(z_2[1]);
}
generated quantities {
  // actual population-level intercept
  // real b_Intercept = Intercept - dot_product(means_X, b);
  real b_Intercept = mod_circular( Intercept - dot_product(means_X, b) );
  // actual population-level intercept
  real b_kappa_Intercept = Intercept_kappa - dot_product(means_X_kappa, b_kappa);
   // actual group-level mu //syntax error?
  vector[N_1] animal_mu_offset = mod_circular_array(r_1_1, N_1);// find offset in [-pi, pi]
  // actual group-level rho //syntax error?
  real animal_mu_rho = rho_circular(r_1_1, N_1);//how do I extract this?
}

