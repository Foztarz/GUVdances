//
// This Stan program defines a simple model, with a
// array of values 'y' modeled as normally distributed
// with mean 'mu' and standard deviation 'sigma'.
//
// Learn more about model development with Stan at:
//
//    http://mc-stan.org/users/interfaces/rstan.html
//    https://github.com/stan-dev/rstan/wiki/RStan-Getting-Started
//

functions {
    //prior for a unit_array estimate of mean angle
   real von_mises_robust_lpdf(real angle, real mu, real kappa) {
        //extreme kappa correction for numeric stability
     if (kappa < 100) {
       return von_mises_lpdf(angle | mu, kappa);
     } else {
       return normal_lpdf(angle | mu, sqrt(1 / kappa));
     }
   }
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

// The input data is a array 'y' of length 'N'.
data {
  int<lower=1> N;  // number of observations
  array[N] real Y;  // response variable, angles
  // data for group-level effects of ID 1
  // effects on mu (added to temp_Intercept)
  int<lower=1> N_1;//number of groups
  int<lower=1> M_1;//number of random effects SD
  array[N] int<lower=1> J_1;//group identity
  array[N] real Z_1_1;//offset from population mean in SDs
  // data for group-level effects of ID 2
  // effects on kappa (added to temp_kappa_Intercept)
  int<lower=1> N_2;
  int<lower=1> M_2;
  array[N] int<lower=1> J_2;
  array[N] real Z_2_kappa_1;
  int prior_only;  // should the likelihood be ignored?
}
// The parameters accepted by the model. Our model
// accepts two parameters 'mu' and 'sigma'.
parameters {
  //A unit_array is used in place of temp_Intercept
  // this stops estimates of mu from walking around the circle multiple times
  real temp_Intercept;  // temporary intercept
  real temp_kappa_Intercept;  // temporary intercept
  // there is only one random effect (animal),
  array[N_1, M_1] real <lower=-pi(), upper=pi()> z_1;//[M_1];  // group-level mean angles effects
  array[M_2] real <lower=0> sd_2;  // group-level standard deviations
  array[N_2] real z_2[M_2];  // unscaled group-level effects
}


transformed parameters {
	// some kind of angle
	real modangle = mod_circular(temp_Intercept);
  // group-level effects (mu)
  array[N_1] real r_1_1 = z_1[1];
  // group-level effects (kappa)
  array[N_2] real r_2_kappa_1 = (sd_2[1] * (z_2[1]));
}
// The model to be estimated. We model the output
model {
  //population level circular mean estimate for mu_vec prior
  real cmean = mean_circular(Y,N);
  array[N] mu = modangle + rep_array(0, N);
  array[N] kappa = temp_kappa_Intercept + rep_array(0, N);

  for (n in 1:N) {
    mu[n] += r_1_1[J_1[n]] * Z_1_1[n];
    mu[n] = mod_circular(mu[n]); // find mu in [-pi, pi]
    kappa[n] += r_2_kappa_1[J_2[n]] * Z_2_kappa_1[n];
    kappa[n] = exp(kappa[n]); // kappa is estimated via the log link, avoiding kappa < 0
  }
  // priors including all constants
  target += uniform_lpdf(modangle | -pi(), pi()); //anything within one circle
  target += normal_lpdf(temp_kappa_Intercept | 0, 3); //most likely values of log(kappa)
  target += student_t_lpdf(sd_2 | 1 , 0, 3);
  target += normal_lpdf(z_2[1] | 0, 1);
  //mu_vec needs some kind of prior, but I'm not entirely sure what
  // My current theory, the population mean needs to be anchored to estimate
  // random effects means well. Because the von Mises distribution wraps,
  // it is possible to get good fits for multiple combinations of
  // population mean, random effects SD and random effects offsets.
  // A strong prior that the population mean should fall at the
  // pop-mean should be close to the mean of the random effects means.
  target += von_mises_robust_lpdf(modangle| mean_circular(z_1[1],N_1), 5); // strong prior
  // this has been proposed (https://discourse.mc-stan.org/t/von-mises-documentation-suggestion/1133/31)
  // likelihood including all constants
  if (!prior_only) {
    for (n in 1:N) {
      target += von_mises_robust_lpdf(Y[n] | mu[n], kappa[n]);
    }
  }
}

generated quantities {
  // return population-level mean angle component
  real b_mu0 = modangle;
  // actual population-level intercept
  real b_Intercept = temp_Intercept;//normally this would be temp_Intercept
  // actual population-level intercept
  real b_kappa_Intercept = temp_kappa_Intercept;
  // actual group-level mu //syntax error?
  array[N_1] real animal_mu_offset = mod_circular_array(r_1_1, N_1);// find offset in [-pi, pi]
  // actual group-level rho //syntax error?
  real animal_mu_rho = rho_circular(r_1_1, N_1);//how do I extract this?
  // actual group-level kappa
  array[N_2] real animal_kappa_offset = r_2_kappa_1;
}
