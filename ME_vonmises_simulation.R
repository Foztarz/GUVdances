#FOR A 'CLEAN' RUN, PRESS ctrl+shift+F10 to RESTART Rstudio
graphics.off()
# Details ---------------------------------------------------------------
#       AUTHOR:	James Foster              DATE: 2024 09 26
#     MODIFIED:	James Foster              DATE: 2024 09 26
#
#  DESCRIPTION: Simulate data from a set of von Mises distributions and 
#               fit a heirarchical von Mises model.
#               Modified from: ME_vonmises_test.R
#               
#       INPUTS: 
#               
#      OUTPUTS: Plots and test statistics
#
#	   CHANGES: - 
#
#   REFERENCES: Sayin, S., Graving, J., et al. in revision
#               
#               Gabry J, Češnovar R, Johnson A (2022). 
#               cmdstanr: R Interface to 'CmdStan'.
#               https://mc-stan.org/cmdstanr/
# 
#               Bürkner, P.-C. (2018). 
#               Advanced Bayesian Multilevel Modeling with the R Package brms. 
#               The R Journal 10, 395–411.
# 
#               Carpenter, B., Gelman, A., Hoffman, M. D., Lee, D., Goodrich, B., 
#               Betancourt, M., Brubaker, M., Guo, J., Li, P. and Riddell, A. (2017). 
#               Stan: A Probabilistic Programming Language. 
#               Journal of Statistical Software 76 doi: 10.18637/jss.v076.i01
# 
#
#    EXAMPLES:  
#
# 
#TODO   ---------------------------------------------
#TODO   
#- Histograms in descriptive plots
#- Plot funnel for kappa & mu
#- Interactions model
#- Test range of parameter values
#- Simulate continuous fixed effects


# Set up workspace --------------------------------------------------------


## Load packages ----------------------------------------------------------
#needs installing before first use (in Rstudio, see automatic message)
suppressMessages(#these are disturbing users unnecessarily
  {
    require(circular)#package for handling circular data
    require(CircStats)#package for circular hypothesis tests
    require(brms)#package for preparing Stan models
  }
)


## General functions -----------------------------------------------------

#convert angles to signed angles in (-180, 180)
Mod360.180 = function(x)
{#use atan2 to convert any angle to the range (-180,180)
  deg(
    atan2(y = sin(rad(x)),
          x = cos(rad(x))
    )
  )
}

#the radian equivalent
mod_circular = function(x)
{
  atan2(y = sin(x),
        x = cos(x))
}

#the unwrapped (no discontinuities)
unwrap_circular = function(x)
{
  mux = mean.circular(x = circular(x = x, template = 'none'))
  centx = atan2(y = sin(x - mux),
                x = cos(x  - mux))
  unwrx = centx + mux
}

#invert the softplus link
#https://en.wikipedia.org/wiki/Softplus
inv_softplus = function(x)
{
  log(exp(x)+1) 
}

#convert softplus scaled kappa to mean vector estimate
Softpl_to_meanvec = function(x)
{
  circular::A1(
    inv_softplus(x)
  )
}

#convert circular to normalised
NormCirc = function(x,
                    plusmean = TRUE)
{
  mn = mean.circular(x) * as.numeric(plusmean)
  return(mod_circular(x - mn) + mn)
}
## von Mises model inspection functions ----------------------------------
#function to construct the transformation list
Make_vmtrans = function(mod,
                        angfun = NormCirc,
                        kapfun = Softpl_to_meanvec)
{
  fx_nm = names(fixef(mod)[,1])
  nm_ang = fx_nm[grepl(pattern = 'mu', x = fx_nm)]
  nm_kap = fx_nm[grepl(pattern = 'kappa', x = fx_nm)]
  tlst = c(replicate(n = length(nm_ang),
                     expr = angfun,
                     simplify = FALSE),
           replicate(n = length(nm_kap),
                     expr = kapfun,
                     simplify = FALSE) )
  names(tlst) = paste0('b_', c(nm_ang, nm_kap))
  return(tlst)
}

PlotVMfix = function(mod,
                     angfun = NormCirc,
                     kapfun = Softpl_to_meanvec)
{
  plot(mod,
       nvariables = 10,
       variable = "^b_", 
       transformations = Make_vmtrans(mod, angfun, kapfun),
       regex = TRUE)
}

## Stan variables ---------------------------------------------------------


### Functions ------------------------------------------------------------


#set up the modulo function
mod_circular_fun = stanvar(scode = "
  real mod_circular(real y) {
    return atan2(sin(y), cos(y));
  }
",
                           block = 'functions')
#set up a von Mises PDF that converts to modulo (adapted from BRMS default)
von_mises3_fun = stanvar(scode = "
real von_mises3_lpdf(real y, real mu, real kappa) {
     if (kappa < 100) {
       return von_mises_lpdf(mod_circular(y) | mu, kappa);
     } else {
       return normal_lpdf(mod_circular(y) | mu, sqrt(1 / kappa));
     }
   }
",
                         block = 'functions')
#a circular mean that might be useful
meancirc_fun = stanvar(scode = "
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
",
                       block = 'functions')
#could unwrap estimate; not especially useful as only one estimate per iteration
unwrapcirc_fun = stanvar(scode = "
   // calculate the unwrapped version (no discontinuities)
vector unwrap_circular(vector y, int N)
{

  mux = mean_circular(y, N);
  centx = y - mux;
  for(n in 1:N)
  {
  centx[n] = mod_circular(centx[n])
  }
  unwrx = centx + mux
  return(unwrx)
}
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
",
                         block = 'functions')

stan_mvm_fun = meancirc_fun + 
                mod_circular_fun + 
                unwrapcirc_fun +
                von_mises3_fun

### Parameters -----------------------------------------------------------


#generate modulo outputs of fixed effects
mu_gen = stanvar(scode = "
real mu_circ = mod_circular(b_muangle[1]);
real mu_offs = mod_circular(b_muangle[2]);
",
                 block = 'genquant')
# #or maybe?
# mu_gen = stanvar(scode = "
# vector[size(b_muangle)] mu_circ; //modulo circular estimate
# for (i in 1:size(b_muangle)){
# mu_circ[i] = mod_circular(b_muangle[i])
# }
# ",
#                  block = 'genquant')

#random effects on mean angle
zmu_var = stanvar(scode = "
vector[K_zmu] zmu_id;  // regression coefficients;
for (i in 1:K_zmu){
  zmu_id[i] = mod_circular(b_zmu[i]); //each in modulus format
}
", 
                  block = 'genquant')

#concentration of random effects on mean angle
zkappa_var = stanvar(scode = "
real zkappa;", 
                     block = "parameters") + #define in the parameters block
              stanvar(scode = "
real kappa_id = log1p_exp(zkappa);
          ", 
          block = 'genquant') #inverse softplus (estimated on softplus scale)

zmu_var_slope = stanvar(scode = "
vector[K_zmu] zmu_id;  // regression coefficients;
vector[N_1] zmu_id_condition;  // condition on coefficients per indiv;
for (i in 1:N_1){
zmu_id[i] = mod_circular(b_zmu[i]); //each in modulus format
zmu_id_condition[i] = mod_circular(b_zmu[i+N_1]); //each in modulus format
}
          ", 
                        block = 'genquant')

#concentration of random effects on fixed effect on mean angle
zkappa_var_slope = stanvar(scode = "
real zkappa_condition;",
                           block = "parameters") + 
                     stanvar(scode = "
real kappa_id_condition = log1p_exp(zkappa+zkappa_condition);
          ", 
                            block = 'genquant')



stanvars_intercepts = stan_mvm_fun + zkappa_var + zmu_var
stanvars_slopes = stan_mvm_fun + zkappa_var + zkappa_var_slope + 
                  zmu_var_slope #includes intercepts



# Kappa investigation -----------------------------------------------------
r_aligned = rvonmises(1e4,
               circular(0),
               kappa = 1e6)
r_uniform = rcircularuniform(1e4)

mle.vonmises(r_aligned,bias = TRUE)$kappa
# 1014180
mle.vonmises(r_uniform, bias = TRUE)$kappa
# 0.03146355
mle.vonmises(sample(x = c(r_aligned, r_uniform), 
                    size = 1e4,
                    replace = FALSE),
                    bias = TRUE)$kappa
# 1.150822
mle.vonmises(c(sample(r_aligned,
                      size =  round(1e4 * 2/(1+2)),
                      replace = FALSE),
               sample(r_uniform,
                      size = round(1e4 * 1/(1+2)),
                      replace =  FALSE) ),
                      bias = TRUE)$kappa
# 1.862195
mle.vonmises(c(sample(r_aligned,
                      size =  round(1e4 * 3/(1+3)),
                      replace = FALSE),
               sample(r_uniform,
                      size = round(1e4 * 1/(1+3)),
                      replace =  FALSE) ),
             bias = TRUE)$kappa
# 2.373249
Est_kappa_ratio = function(ratio,
                           n = 1e4,
                           maxkappa = 1e6,
                           bias = TRUE)
{
  mle.vonmises(x = c(rvonmises(n = round(n * ratio/(1+ratio)),
                        mu = circular(0, template = 'none'),
                        kappa = maxkappa),
                     rcircularuniform(n = round(n * 1/(1+ratio)) )
                               ),
               bias = bias)$kappa
}

rats = seq(from = 0,
           to  = 100,
           length.out = 1e3)
ekappas = sapply(X = rats,
                 FUN = Est_kappa_ratio)

#unimodal/uniform as a predictor of kappa estimates
plot(x =rats,
     y = ekappas,
     col = gray(0.5, 0.1),
     pch = 19)
abline(a = 0, b = 1,
       h = 0, v = 0,
       col = c(2,1,1))
abline( lm(ekappas ~rats),
        col = 5,
        lwd = 3)
#on a complementary log scale we see where the relationship becomes linear
plot(x =rats[-1],
     y = ekappas[-1],
     log = 'xy',
     col = gray(0.5, 0.1),
     pch = 19)
abline(a = 0, b = 1,
       col = 2)
lines(x = rats[-1],
      y = exp(predict(lm(log(ekappas[-1]) ~log(rats[-1]))) ),
      col = 5,
      lwd = 3)
#on a logarithmic scale we can interpret ratios as a predictor of SD
plot(x =rats,
     y = deg(sqrt(1/ekappas)),
     log = 'xy',
     col = gray(0.5, 0.1),
     pch = 19)
lines(x = rats[-1],
      y = exp( predict(
        lm(formula = log(deg(sqrt(1/ekappas[-1])) )~log(rats[-1]) )
      ) ),
      col = 5,
      lwd = 3)
#on converted to mean vectors, the log ratio becomes a good predictor
plot(x =rats[-1],
     y = A1(ekappas[-1]),
     log = 'x',
     ylim = c(0,1),
     col = gray(0.5, 0.1),
     pch = 19)
abline(h = c(0,1))
lines(x = rats[-1],
      y = plogis(log(rats[-1])),
        col = 2,
      lwd = 1)
lines(x = rats[-1],
      y = predict(
        glm(formula = A1(ekappas[-1])~log(rats[-1]), 
            family = quasibinomial(link = 'logit')),
        type = 'response'
        ),
        col = 5,
      lwd = 3,
      lty = 3)

plot(x =rats,
     y = A1(ekappas),
     ylim = c(0,1),
     col = gray(0.5, 0.1),
     pch = 19)
abline(h = c(0,1), v = 0)
lines(x = rats[-1],
      y = plogis(log(rats[-1])),
      col = 2,
      lwd = 1)
lines(x = rats[-1],
      y = predict(
        glm(formula = A1(ekappas[-1])~log(rats[-1]), 
            family = quasibinomial(link = 'logit')),
        type = 'response'
      ),
      col = 5,
      lwd = 3,
      lty = 3)

#in fact the relationship is near 1:1,
#the logit of the mean vector is equal to the ratio of unimodal to uniform
#mean vector = 1- 1/(1+ratio)
plot(x = 1-1/(1+rats),
     y = A1(ekappas),
     col = gray(0.5, 0.1),
     pch = 19)
abline(a = 0, b =1, h = 0, v = 0,
       col = c(2,1,1))

plot(x = rats+1,
     y = 1/(1-A1(ekappas)),
     col = gray(0.5, 0.1),
     pch = 19)
abline(a = 0, b =1, h = 0, v = 0,
       col = c(2,1,1))



# Input Variables ----------------------------------------------------------

all_plots = FALSE # to speed up
#  .  User input -----------------------------------------------------------
set.seed(20240905)#day the simulation was updated
n_iter = 1e4 #optimisation iterations

csv_sep = ','#Is the csv comma separated or semicolon separated? For tab sep, use "\t"
angle_name = "angles" #The title of the column with angles; NO SPACES PLEASE
angle_unit = "degrees" # "degrees" or "radians"
angle_rot = "counter" # "clock" or "counter"

#Check the operating system and assign a logical flag (T or F)
sys_win = Sys.info()[['sysname']] == 'Windows'
#User profile instead of home directory
if(sys_win)
{#get rid of all the backslashes
  ltp  =  gsub(pattern = '\\\\', replacement = '/', x = Sys.getenv('USERPROFILE'))#Why does windows have to make this so difficult
}else{#Root directory should be the "HOME" directory on a Mac (or Linux?)
  ltp  = Sys.getenv('HOME')#Life was easier on Mac
}

# Simulate data  ------------------------------------------------
n_indiv = 50#10 # Number of random-effects groups
n_trials = 10 # Number of observations per random-effects group
n_conditions = 2 # Number of different conditions
# minimum discriminable angle appears to be approx 35°
kappa_pop = A1inv(0.7) #concentration around each trial mean
logkappa_var = 0.5 #scale of random variation in concentration (log units)
logkappa_indiv_var = 0.1 #scale of random variation in concentration (log units)
mu_0 = circular(x = 175,
                units = angle_unit,
                rotation = angle_rot
                # mu_0 = rcircularuniform(n = 1,
                #                         control.circular = list(units = angle_unit,
                #                                                 rotation = angle_rot)
)#population intercept
mu_offset = c(0,
              circular(20, 
                       units = angle_unit,
                       rotation = angle_rot)
              # rcircularuniform(n = n_conditions - 1, 
              #                  control.circular = list(units = angle_unit,
              #                                          rotation = angle_rot)
)#condition change in heading
#condition change in precision
lk_offset = c(0, 
              rnorm(n = n_conditions-1,
                    mean = 0,
                    sd = 0)
)
#kappa across individuals
lk_indiv = log( A1inv(0.75) ) #log( A1inv(0.98) ) #concentration across individuals (pairs)
#kappa across individual offsets
lk_indiv_offset = log( A1inv(0.90) ) #concentration across individual offsets
#trials
trials = rep(x = 1:n_trials, times = n_indiv*n_conditions)
#individuals
indivs = rep( x = sort( rep(1:n_indiv, times = n_trials) ),
              times = n_conditions)
#conditions
conds = sort( rep(1:n_conditions - 1, times = n_indiv*n_trials) )
#mean angle in condition 2 for each individual 
mu0_sim = rvonmises(n = n_indiv,
                    mu = mu_0,#random angle
                    kappa = exp(lk_indiv)#the wider the distribution of individual biases, the greater the influence of pairing
)
#mean angle chagne to condition 2 for each individual 
mu1_sim = rvonmises(n = n_indiv,
                    mu = circular(0, 
                                  units = angle_unit,
                                  rotation = angle_rot),#population mean
                    kappa = exp(lk_indiv_offset)#the wider the distribution of individual biases, the greater the influence of pairing
)
#log precision in trail 1 for each individual 
lk0_sim = rnorm(n = n_indiv,
                mean = log(kappa_pop),
                sd = logkappa_var)
#log precision change for condition 2 for each individual 
lk1_sim = rnorm(n = n_indiv,
                mean = 0,
                sd = logkappa_indiv_var)
#TODO check this works!
#all conditional means
mu_all =  mu0_sim[indivs] + 
  mu_offset[conds+1] + 
  conds*mu1_sim
lk_all =  lk0_sim[indivs] + 
  lk_offset[conds+1] + 
  conds*lk1_sim
# mu_all =  c( 
#               t( 
#                 rep(1, times = n_conditions) %*% t(mu0_sim) + #repeat indiv mean across conditions
#                 rep(c(0,mu_offset), times = n_conditions) #add condition mean
#                 )#collapse lengthways
#             )
# #all conditional kappas
# lk_all =  c( 
#               t( 
#                 rep(1, times = n_conditions) %*% t(lk0_sim) + #repeat indiv mean across conditions
#                 rep(c(0,kappa_offset), times = n_conditions) #add condition mean
#                 )#collapse lengthways
#             )

#simulate the full dataset
sim = data.frame(
  indiv = indivs,
  trial = trials,
  condition = conds,
  angle = round( #round to the nearest angle (typical resolution)
    mapply(FUN = rvonmises,
           n = 1,
           mu = lapply(X = mu_all,
                       FUN = circular,
                       units = angle_unit,
                       rotation = angle_rot),
           kappa = exp(lk_all)
    )
  )
)
# #save somewhere the user likely keeps data
# write.table(x = sim,
#             file = file.path(ltp,'Documents', "simulated_angles.csv"),
#             sep = csv_sep,
#             row.names = FALSE
#             )
dt_dim = n_indiv*n_trials

# Plot simulated data -----------------------------------------------------
# if(all_plots)
# {
par(mar =rep(0,4))
with(subset(sim, condition %in% 0),
     {
       plot.circular(x = circular(x = angle, 
                                  type = 'angles',
                                  unit = angle_unit,
                                  template = 'geographics',
                                  modulo = '2pi',
                                  zero = pi/2,
                                  rotation = angle_rot
       ),
       stack = TRUE,
       bins = 360/5,
       sep = 0.5/dt_dim,
       col = 'cyan4'
       )
     }
)
par(new = T)
with(subset(sim, condition %in% 1),
     {
       plot.circular(x = circular(x = angle, 
                                  type = 'angles',
                                  unit = 'degrees',
                                  template = 'geographics',
                                  modulo = '2pi',
                                  zero = pi/2,
                                  rotation = angle_rot
       ),
       stack = TRUE,
       bins = 360/5,
       sep = -0.5/dt_dim[1],
       col = 'darkblue',
       shrink = 1.05,
       axes = F
       )
     }
)
# }
# Organise data ----------------------------------------------------------
IDs = paste0(sort(rep(LETTERS, times = 26)), letters)
sim = within(sim,
             {ID = IDs[indivs]}
)

sim = within(sim,
             {
               rad_angle = rad(Mod360.180(angle))
             }
)


# Nonlinear Stan version --------------------------------------------------
#load packages
require(cmdstanr)
require(brms)
