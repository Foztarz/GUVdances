#FOR A 'CLEAN' RUN, PRESS ctrl+shift+F10 to RESTART Rstudio
graphics.off()
# Details ---------------------------------------------------------------
#       AUTHOR:	James Foster              DATE: 2024 08 28
#     MODIFIED:	James Foster              DATE: 2024 09 12
#
#  DESCRIPTION: Fit hierarchical maximum-likelihood von Mises.
#               
#       INPUTS: 
#               
#      OUTPUTS: Plots and test statistics
#
#	   CHANGES: - More comprehensive simulation
#             - Added von Mises random effects      
#             - t-dist von Mises random effects      
#             - Weibull dist von Mises random effects kappa
#             - plotting predictions for random effects (working well)
#             - model hypothesis testing
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
#- Test for angles close to 180 +
#- Compare circular and Gaussian raneff +
#- Softplus raneff kappa  +
#- Plot predictions +
#- Check hypothesis tests +
#- Model comparison hypothesis tests + 
#- Test with more individuals
#- Unwrap circular estimates
#- Plot estimates against simulations
#- Extract circular fixed effects medians
#- Simulate random slopes
#- Model random slopes
#- Simulate continuous fixed effects

# . Load packages ----------------------------------------------------------
#needs installing before first use (in Rstudio, see automatic message)
suppressMessages(#these are disturbing users unnecessarily
  {
    require(circular)#package for handling circular data
    require(CircStats)#package for circular hypothesis tests
    require(brms)#package for preparing Stan models
  }
)


# . General functions -----------------------------------------------------

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
# . von Mises model inspection functions ----------------------------------
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
mu_0 = circular(x = 175,
                units = angle_unit,
                rotation = angle_rot
# mu_0 = rcircularuniform(n = 1,
#                         control.circular = list(units = angle_unit,
#                                                 rotation = angle_rot)
)#population intercept
mu_offset = c(0,
              circular(5, 
                       units = angle_unit,
                      rotation = angle_rot)
              # rcircularuniform(n = n_conditions - 1, 
              #                  control.circular = list(units = angle_unit,
              #                                          rotation = angle_rot)
              )#condition change in heading
lk_offset = c(0, 
              rnorm(n = n_conditions-1,
                    mean = 0,
                    sd = 0)
)#condition change in precision
lk_indiv = log( A1inv(0.98) ) #concentration across individuals (pairs)
#trials
trials = rep(x = 1:n_trials, times = n_indiv*n_conditions)
#individuals
indivs = rep( x = sort( rep(1:n_indiv, times = n_trials) ),
              times = n_conditions)
#conditions
conds = sort( rep(1:n_conditions - 1, times = n_indiv*n_trials) )
#mean angle in trail 1 for each individual 
mu0_sim = rvonmises(n = n_indiv,
                    mu = mu_0,#random angle
                    kappa = exp(lk_indiv)#the wider the distribution of individual biases, the greater the influence of pairing
)
#log precision in trail 1 for each individual 
lk0_sim = rnorm(n = n_indiv,
                mean = log(kappa_pop),
                sd = logkappa_var)
#all conditional means
mu_all =  mu0_sim[indivs] + mu_offset[conds+1]
lk_all =  lk0_sim[indivs] + lk_offset[conds+1]
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

#set up model fit
formula_nlvm = bf(
  #set up a formula for the curve as a whole,
  formula = rad_angle ~ mod_circular(muangle), #working on calling this mu
  muangle ~ condition, #N.B. this is similar to the slope, so all of its effects depend on stimulus level
  family = von_mises(link = "identity",
                     link_kappa = "softplus"),
  nl = TRUE)

prior_nlvm = get_prior(formula = formula_nlvm,
          data = sim)

prior_nlvm = within(prior_nlvm,
                    {
                    prior[nlpar %in% 'muangle'] = 'normal(0, pi())'
                    }
                    )
#set up the modulo function
mod_circular_fun = "
  real mod_circular(real y) {
    return atan2(sin(y), cos(y));
  }
"
#set up a von Mises PDF that converts to modulo (adapted from BRMS default)
von_mises3_fun = "
real von_mises3_lpdf(real y, real mu, real kappa) {
     if (kappa < 100) {
       return von_mises_lpdf(mod_circular(y) | mu, kappa);
     } else {
       return normal_lpdf(mod_circular(y) | mu, sqrt(1 / kappa));
     }
   }
"
#a circular mean that might be useful
meancirc_fun = "
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
"
#could unwrap estimate; not especially useful as only one estimate per iteration
unwrapcirc_fun = "
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
"
#generate modulo outputs
mu_gen <- "
real mu_circ = mod_circular(b_muangle[1]);
real mu_offs = mod_circular(b_muangle[2]);
"
  
stan_var = stanvar(scode = mod_circular_fun, block = "functions") + 
  stanvar(scode = von_mises3_fun, block = "functions") + 
              stanvar(scode = mu_gen, block = "genquant")

# . Short dummy run to check the influence of the priors ------------------


#double check that the prior distribution is viable by first setting up a short dummy run
# Dummy run
system.time(
  {
    dummy_fit = brm( formula = formula_nlvm, # using our nonlinear formula
                     data = sim, # our data
                     prior = prior_nlvm, # our priors 
                     stanvars = stan_var,
                     sample_prior = 'only', #ignore the data to check the influence of the priors
                     iter = 300, # short run for 300 iterations
                     chains = 4, # 4 chains in parallel
                     cores = 4, # on 4 CPUs
                     refresh = 0, # don't echo chain progress
                     backend = 'cmdstanr') # use cmdstanr (other compilers broken)
  }
)
# On my computer this takes <60s, each chain running for <1 seconds (mainly compile time)
if(all_plots)
{
plot(dummy_fit,
     variable = '^mu_',
     regex = TRUE)
plot(dummy_fit,
     variable = '^kappa',
     regex = TRUE)
}
#conditional_effects(dummy_fit)


# . To check that estimation works  ---------------------------------------
# Run the model for a small number of iterations to check that it is possible to 
# estimate all parameters. This may be the point where we encounter numerical errors
# if the formula or priors are misspecified.
# e.g. if the formula returns estimates of correct choice rate outside of [0,1].

# Short run
system.time(
  {
    short_fit = brm( formula = formula_nlvm, # using our nonlinear formula
                     data = sim, # our data
                     prior = prior_nlvm, # our priors 
                     stanvars = stan_var,
                     iter = 300, # short run for 300 iterations (less than 300 gives insufficient warmup time)
                     chains = 4, # 4 chains in parallel
                     cores = 4, # on 4 CPUs
                     refresh = 0, # don't echo chain progress
                     backend = 'cmdstanr') # use cmdstanr (other compilers broken)
  }
)
if(all_plots)
{
plot(short_fit)
plot(short_fit,
     variable = '^mu_',
     regex = TRUE)
summary(short_fit)
}
# plot(
#   #conditional_effects(x = short_fit, 
#                       spaghetti = TRUE, 
#                       ndraws = 2e2,
#                       effects = 'condition')
# )

# . Full runs -------------------------------------------------------------


# Full run
system.time(
  {
    full_fit = brm( formula = formula_nlvm, # using our nonlinear formula
                     data = sim, # our data
                     prior = prior_nlvm, # our priors 
                     stanvars = stan_var,
                     iter = 1000, # short run for 1000 iterations (less than 300 gives insufficient warmup time)
                     chains = 4, # 4 chains in parallel
                     cores = 4, # on 4 CPUs
                     refresh = 0, # don't echo chain progress
                     backend = 'cmdstanr') # use cmdstanr (other compilers broken)
  }
)

if(all_plots){
plot(full_fit,
     variable = '^mu_',
     regex = TRUE)
plot(full_fit)
summary(full_fit)
}
# plot(
#   #conditional_effects(x = short_fit, 
#                       spaghetti = TRUE, 
#                       ndraws = 2e2,
#                       effects = 'condition')
# )


# Linear model sanity check -----------------------------------------------
#set up model fit
formula_lm = bf(
  #set up a formula for the curve as a whole,
  formula = rad_angle ~ condition + (1|ID),
  nl = FALSE)#the joint distribution for these parameters is undefined, and therefore the parameters themselves are "nonlinear"

prior_lm = get_prior(formula = formula_lm,
                         data = sim)

prior_lm = within(prior_lm,
                      {
                        prior[prior %in% ''] = 'normal(0, 3)'
                      }
)


# . Short dummy run to check the influence of the priors ------------------


#double check that the prior distribution is viable by first setting up a short dummy run
# Dummy run
system.time(
  {
    dummy_fitlm = brm( formula = formula_lm, # using our nonlinear formula
                       data = sim, # our data
                       prior = prior_lm, # our priors 
                       sample_prior = 'only', #ignore the data to check the influence of the priors
                       iter = 300, # short run for 300 iterations
                       chains = 4, # 4 chains in parallel
                       cores = 4, # on 4 CPUs
                       refresh = 0, # don't echo chain progress
                       backend = 'cmdstanr') # use cmdstanr (other compilers broken)
  }
)


# . Short run to check the convergence ------------------


# short run
system.time(
  {
    short_fitlm = brm( formula = formula_lm, # using our nonlinear formula
                       data = sim, # our data
                       prior = prior_lm, # our priors 
                       iter = 300, # short run for 300 iterations
                       chains = 4, # 4 chains in parallel
                       cores = 4, # on 4 CPUs
                       refresh = 0, # don't echo chain progress
                       backend = 'cmdstanr') # use cmdstanr (other compilers broken)
  }
)# . Full run to check the convergence ------------------


# full run
system.time(
  {
    full_fitlm = brm( formula = formula_lm, # using our nonlinear formula
                       data = sim, # our data
                       prior = prior_lm, # our priors 
                       iter = 1000, # full run for 1000 iterations
                       chains = 4, # 4 chains in parallel
                       cores = 4, # on 4 CPUs
                       refresh = 0, # don't echo chain progress
                       backend = 'cmdstanr') # use cmdstanr (other compilers broken)
  }
)

# . Plot LM with simulated data -----------------------------------------------------
if(all_plots)
{
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


arrows.circular(x = circular(mu_0,  
                             type = 'angles',
                             unit = 'degrees',
                             modulo = '2pi',
                             zero = pi/2,
                             template = 'geographics',
                             rotation = angle_rot),
                y = A1(kappa_pop),
                lwd = 1,
                length = 0,
                col = adjustcolor('cyan4', alpha.f = 1.0)
)
arrows.circular(x = circular(mu_0 + mu_offset[2],  
                             type = 'angles',
                             unit = 'degrees',
                             modulo = '2pi',
                             zero = pi/2,
                             template = 'geographics',
                             rotation = angle_rot),
                y = A1(exp(log(kappa_pop) + lk_offset[2])),
                lwd = 1,
                length = 0,
                col = adjustcolor('blue2', alpha.f = 1.0)
)
}
sm_lm = summary(full_fitlm, robust = TRUE)
prms_lm = with(sm_lm, rbind(fixed, spec_pars))
est_lm = data.frame(t(t(prms_lm)['Estimate',]))
if(all_plots)
{
with(est_lm,
     {
       arrows.circular(x = circular(deg(Intercept ),  
                                    type = 'angles',
                                    unit = 'degrees',
                                    modulo = '2pi',
                                    zero = pi/2,
                                    template = 'geographics',
                                    rotation = angle_rot),
                       y = 1/(sigma^2),
                       lwd = 5,
                       length = 0,
                       col = adjustcolor('cyan2', alpha.f = 0.5)
       )
     }
)
with(est_lm,
     {
       arrows.circular(x = circular(deg(Intercept  + condition),  
                                    type = 'angles',
                                    unit = 'degrees',
                                    modulo = '2pi',
                                    zero = pi/2,
                                    template = 'geographics',
                                    rotation = angle_rot),
                       y = 1/(sigma^2),
                       lwd = 5,
                       length = 0,
                       col = adjustcolor('blue2', alpha.f = 0.5)
       )
     }
)
}


# Nonlinear ME Stan version --------------------------------------------------


#set up model fit
formula_nlvmme = bf(
  #set up a formula for the curve as a whole,
  formula = rad_angle ~ mod_circular(muangle),
  muangle ~  condition + (1|ID), #linear random effects that will be converted to (-pi,pi)
  kappa ~ condition + (1|ID), #linear random effects that will be converted to inv_softplus(kappa)
  family = von_mises(link = "identity",
                     link_kappa = 'softplus'),
  nl = TRUE)#the joint distribution for these parameters is undefined, and therefore the parameters themselves are "nonlinear"

prior_nlvmme = get_prior(formula = formula_nlvmme,
          data = sim)

prior_nlvmme = within(prior_nlvmme,
                    {
                    prior[nlpar %in% 'muangle' & coef %in% 'Intercept'] = 'normal(0, pi())'
                    prior[nlpar %in% 'muangle' & class %in% 'b'] = 'normal(0, pi())'
                    prior[dpar %in% 'kappa' & class %in% 'Intercept'] = 'normal(0.5, 1.0)'
                    prior[dpar %in% 'kappa' & class %in% 'b'] = 'normal(0.0, 1.0)'
                    # prior[nlpar %in% 'muangle'] = 'von_mises(0, 0.5)'
                    # lb[nlpar %in% 'muangle'] = -pi
                    # ub[nlpar %in% 'muangle'] = pi
                    }
                    )

# . Short dummy run to check the influence of the priors ------------------


#double check that the prior distribution is viable by first setting up a short dummy run
# Dummy run
system.time(
  {
    dummy_fitme = brm( formula = formula_nlvmme, # using our nonlinear formula
                     data = sim, # our data
                     prior = prior_nlvmme, # our priors 
                     stanvars = stan_var,
                     sample_prior = 'only', #ignore the data to check the influence of the priors
                     iter = 300, # short run for 300 iterations
                     chains = 4, # 4 chains in parallel
                     cores = 4, # on 4 CPUs
                     refresh = 0, # don't echo chain progress
                     backend = 'cmdstanr') # use cmdstanr (other compilers broken)
  }
)
# On my computer this takes <60s, each chain running for <1 seconds (mainly compile time)
if(all_plots){
plot(dummy_fitme,
     variable = '^mu_',
     regex = TRUE)
plot(dummy_fitme)
}
#conditional_effects(dummy_fit, spaghetti = TRUE)


# . To check that estimation works  ---------------------------------------
# Run the model for a small number of iterations to check that it is possible to 
# estimate all parameters. This may be the point where we encounter numerical errors
# if the formula or priors are misspecified.
# e.g. if the formula returns estimates of correct choice rate outside of [0,1].

# Short run
system.time(
  {
    short_fitme = brm( formula = formula_nlvmme, # using our nonlinear formula
                     data = sim, # our data
                     prior = prior_nlvmme, # our priors 
                     stanvars = stan_var,
                     iter = 300, # short run for 300 iterations (less than 300 gives insufficient warmup time)
                     chains = 4, # 4 chains in parallel
                     cores = 4, # on 4 CPUs
                     refresh = 0, # don't echo chain progress
                     backend = 'cmdstanr') # use cmdstanr (other compilers broken)
  }
)
if(all_plots)
{
plot(short_fitme,
     variable = '^mu_',
     regex = TRUE)
plot(short_fitme, 
     nvariables = 10,
     variable = "^b_", 
     regex = TRUE)
summary(short_fitme)
}
# plot(
  #conditional_effects(x = short_fit, 
#                       spaghetti = TRUE, 
#                       ndraws = 2e2,
#                       effects = 'condition')
# )

# . Full runs -------------------------------------------------------------


# Full run
system.time(
  {
    full_fitme = brm( formula = formula_nlvmme, # using our nonlinear formula
                     data = sim, # our data
                     prior = prior_nlvmme, # our priors 
                     stanvars = stan_var,
                     iter = 1000, # short run for 1000 iterations (less than 300 gives insufficient warmup time)
                     chains = 4, # 4 chains in parallel
                     cores = 4, # on 4 CPUs
                     init = 0,
                     refresh = 0, # don't echo chain progress
                     backend = 'cmdstanr') # use cmdstanr (other compilers broken)
  }
)

if(all_plots)
{
plot(full_fitme,
     variable = "^mu_", 
     regex = TRUE)
}
fx_names = names(fixef(full_fitme)[,1])
nm_angle_vars = fx_names[grepl(pattern = 'angle', x = fx_names)]
trans_lst = replicate(n = length(nm_angle_vars),
                      expr = mod_circular,
                      simplify = FALSE)
names(trans_lst) = paste0('b_', nm_angle_vars)



# . von Mises specific diagnostics ----------------------------------------
# #convert softplus-kappa to mean vector estimate
# A1sp = function(x){A1(inv_softplus(x))}
# #convert circular to normalised
# NZmod_circular = function(x)
# {
#   mn = mean.circular(x)
#   return(mod_circular(x - mn) + mn)
# }
# #function to construct the transformation list
# Make_vmtrans = function(mod,
#                         angfun = NZmod_circular,
#                         kapfun = A1sp)
# {
#   fx_nm = names(fixef(mod)[,1])
#   nm_ang = fx_nm[grepl(pattern = 'mu', x = fx_nm)]
#   nm_kap = fx_nm[grepl(pattern = 'kappa', x = fx_nm)]
#   tlst = c(replicate(n = length(nm_ang),
#                           expr = angfun,
#                           simplify = FALSE),
#                 replicate(n = length(nm_kap),
#                           expr = kapfun,
#                           simplify = FALSE) )
#   names(tlst) = paste0('b_', c(nm_ang, nm_kap))
#   return(tlst)
# }
# 
# PlotVMfix = function(mod,
#                      angfun = NZmod_circular,
#                      kapfun = A1sp)
# {
#   plot(mod,
#        nvariables = 10,
#        variable = "^b_", 
#        transformations = Make_vmtrans(mod, angfun, kapfun),
#        regex = TRUE)
# }
if(all_plots)
{
PlotVMfix(full_fit)
PlotVMfix(full_fitme)


plot(full_fitme,
     nvariables = 10,
     variable = "^b_", 
     transform = trans_lst,
     regex = TRUE)
plot(full_fitme,
     nvariables = 10,
     variable = "^sd_", 
     regex = TRUE)

summary(full_fitme)
}
rhat(full_fitme,
     pars = "^mu_", 
     regex = TRUE)

# plot(
  #conditional_effects(x = short_fit, 
#                       spaghetti = TRUE, 
#                       ndraws = 2e2,
#                       effects = 'condition')
# )

# Nonlinear vM-ME Stan version --------------------------------------------------


#set up model fit
formula_nlvmmevm = bf(
  #set up a formula for the curve as a whole,
  formula = rad_angle ~ mod_circular(muangle + zmu),
  muangle ~  condition, #N.B. this is similar to the slope, so all of its effects depend on stimulus level
  zmu ~  0+ID, #N.B. this is similar to the slope, so all of its effects depend on stimulus level
  kappa ~ condition + (1|ID), #N.B. this is similar to the slope, so all of its effects depend on stimulus level
  family = von_mises(link = "identity",
                     link_kappa = 'softplus'),
  nl = TRUE)#the joint distribution for these parameters is undefined, and therefore the parameters themselves are "nonlinear"

prior_nlvmmevm = get_prior(formula = formula_nlvmmevm,
                         data = sim)

prior_nlvmmevm = within(prior_nlvmmevm,
                      {
                        prior[nlpar %in% 'muangle' & coef %in% 'Intercept'] = 'von_mises3(0, 0.1)'#'von_mises3(0, 1.5)'#'normal(0, pi())'
                        prior[nlpar %in% 'muangle' & class %in% 'b'] = 'von_mises3(0, 1.5)'#'normal(0, pi())'
                        prior[nlpar %in% 'zmu' & coef %in% 'Intercept'] = 'von_mises3(0, log1p_exp(zkappa))'
                        prior[nlpar %in% 'zmu' & class %in% 'b'] = 'von_mises3(0, log1p_exp(zkappa))'
                        prior[dpar %in% 'kappa' & class %in% 'Intercept'] = 'normal(0.5, 1.0)'
                        prior[dpar %in% 'kappa' & class %in% 'b'] = 'normal(0.0, 1.0)'
                        # prior[nlpar %in% 'muangle'] = 'von_mises(0, 0.5)'
                        # lb[nlpar %in% 'muangle'] = -pi
                        # ub[nlpar %in% 'muangle'] = pi
                      }
) + 
  # by default sd prior student_t(3,0,2.5) would assume pooling, what is circular equivalent?
 # set_prior("target += normal_lpdf(log1p_exp(zkappa) | 1,3)", 
 #           check = FALSE) 
 set_prior("target += weibull_lpdf(log1p_exp(zkappa) | 1.6, 5)", 
           check = FALSE)
  # set_prior("target += student_t_lpdf(log1p_exp(zkappa) | 1, 0, 0.25)", check = FALSE) 
#circular SD (Mardia, 1972) = sqrt(-2*log(rho))
# if kappa = 1, SD = deg(sqrt(-2*log(A1(1.0)))) = 73°
#moderate pooling sd = 90°, kappa = A1inv(exp((rad(90)^2)/(-2))) = 0.61
#high pooling sd = 30°, kappa = A1inv(exp((rad(30)^2)/(-2))) = 4.21
#extreme pooling sd = 10°, kappa = A1inv(exp((rad(10)^2)/(-2))) = 33.33
xx = seq(from  = -15, to  = 15, length.out = 1e3)
if(all_plots)
{
plot(x = A1(log(1+exp(xx))),
     y = dnorm(log(1+exp(xx)), mean = 1, sd = 3),
     type = 'l',
     xlab = 'rho of kappa_id',
     ylab = 'probability density',
     main = 'zkappa prior, normal(0,3)',
     lwd = 2,
     col = 5,
     xlim = c(0,1),
     ylim = c(0,0.2)
     )
abline(h = 0, v = 0)
abline(v = A1(exp(lk_indiv)), col = 2)
}
#weibull median is shape*(log(2)^(1/scale))
#target median = 1.5
#target scale
wscale = 5
wshape = 1.5/(log(2)^(1/wscale))
if(all_plots)
{
plot(x = log(1+exp(xx)),
     y = dweibull(log(1+exp(xx)), shape = 1.6, scale = 5),
     type = 'l',
     xlab = 'kappa_id',
     ylab = 'probability density',
     main = 'zkappa prior, Weilbull(1,1)',
     lwd = 2,
     col = 6,
     xlim = c(0,10),
     ylim = c(0,0.2)
     )
abline(h = 0, v = 0)

plot(x = A1(log(1+exp(xx))),
     y = dweibull(log(1+exp(xx)), shape = 1.6, scale = 5),
     type = 'l',
     xlab = 'rho of kappa_id',
     ylab = 'probability density',
     main = 'zkappa prior, Weilbull(1,1)',
     lwd = 2,
     col = 6,
     xlim = c(0,1),
     ylim = c(0,0.2)
     )
abline(h = 0, v = 0)
abline(v = A1(exp(lk_indiv)), col = 2)
}
#t-distribution?
if(all_plots)
{
plot(x = A1(log(1+exp(xx))),
     y = brms::dstudent_t(log(1+exp(xx)), df = 3, mu = 1.5, sigma = 10),
     type = 'l',
     xlab = 'rho of kappa_id',
     ylab = 'probability density',
     main = 'zkappa prior, student_t(3,1.5,3)',
     lwd = 2,
     col = 5,
     xlim = c(0,1),
     ylim = c(0,0.2)
)
abline(h = 0, v = 0)
abline(v = A1(exp(lk_indiv)), col = 2)
}   

zkappa_var = stanvar(scode = "  real zkappa;", block = "parameters") + 
  stanvar(scode = "
          real kappa_id = log1p_exp(zkappa);
          ", 
          block = 'genquant')
zmu_var = stanvar(scode = "
           vector[K_zmu] zmu_id;  // regression coefficients;
           for (i in 1:K_zmu){
            zmu_id[i] = mod_circular(b_zmu[i]); //each in modulus format
           }
          ", 
          block = 'genquant')

stanvars = stan_var + zkappa_var + zmu_var

# . Short dummy run to check the influence of the priors ------------------


#double check that the prior distribution is viable by first setting up a short dummy run
# Dummy run
system.time(
  {
    dummy_fitmevm = brm( formula = formula_nlvmmevm, # using our nonlinear formula
                       data = sim, # our data
                       prior = prior_nlvmmevm, # our priors 
                       stanvars = stanvars,
                       sample_prior = 'only', #ignore the data to check the influence of the priors
                       iter = 1000, #300, # short run for 300 iterations
                       chains = 4, # 4 chains in parallel
                       cores = 4, # on 4 CPUs
                       refresh = 0, # don't echo chain progress
                       backend = 'cmdstanr') # use cmdstanr (other compilers broken)
  }
)
# On my computer this takes <60s, each chain running for <1 seconds (mainly compile time)
if(all_plots)
  {
plot(dummy_fitmevm,
     variable = '^mu_',
     regex = TRUE)
plot(dummy_fitmevm,
     variable = '^kappa_',
     regex = TRUE)
plot(dummy_fitmevm,
     variable = '^zmu_id',
     nvariables = 5,
     regex = TRUE)
plot(dummy_fitmevm)
}
#

#

#

#

#

#

#

#conditional_effects(dummy_fitmevm, spaghetti = TRUE)


# . To check that estimation works  ---------------------------------------
# Run the model for a small number of iterations to check that it is possible to 
# estimate all parameters. This may be the point where we encounter numerical errors
# if the formula or priors are misspecified.
# e.g. if the formula returns estimates of correct choice rate outside of [0,1].

# Short run
system.time(
  {
    short_fitmevm = brm( formula = formula_nlvmmevm, # using our nonlinear formula
                       data = sim, # our data
                       prior = prior_nlvmmevm, # our priors 
                       stanvars = stanvars,
                       iter = 300, # short run for 300 iterations (less than 300 gives insufficient warmup time)
                       chains = 4, # 4 chains in parallel
                       cores = 4, # on 4 CPUs
                       refresh = 0, # don't echo chain progress
                       backend = 'cmdstanr') # use cmdstanr (other compilers broken)
  }
)

if(all_plots)
{
plot(short_fitmevm,
     variable = '^mu_',
     regex = TRUE)
plot(short_fitmevm,
     variable = '^kappa_',
     regex = TRUE)
plot(short_fitmevm,
     variable = '^zmu_id',
     nvariables = 10,
     regex = TRUE)


plot(short_fitmevm, 
     nvariables = 10,
     variable = "^b_", 
     regex = TRUE)

summary(short_fitmevm)
}

# plot(
#   #conditional_effects(x = short_fitmevm, 
#                       spaghetti = TRUE, 
#                       ndraws = 2e2,
#                       effects = 'condition')
# )

# . Full runs -------------------------------------------------------------


# Full run
system.time(
  {
    full_fitmevm = brm( formula = formula_nlvmmevm, # using our nonlinear formula
                      data = sim, # our data
                      prior = prior_nlvmmevm, # our priors 
                      stanvars = stanvars,
                      iter = 1000, # short run for 1000 iterations (less than 300 gives insufficient warmup time)
                      chains = 4, # 4 chains in parallel
                      cores = 4, # on 4 CPUs
                      init = 0,
                      refresh = 0, # don't echo chain progress
                      backend = 'cmdstanr') # use cmdstanr (other compilers broken)
  }
)

if(all_plots)
{
plot(full_fitmevm,
     variable = "^mu_", 
     regex = TRUE)
plot(full_fitmevm,
     variable = "^kappa_", 
     regex = TRUE)
plot(full_fitmevm,
     variable = '^zmu_id',
     nvariables = 5,
     regex = TRUE)



plot(full_fitmevm,
     nvariables = 10,
     variable = "^sd_", 
     regex = TRUE)

summary(full_fitmevm)
}
rhat(full_fitmevm,
     pars = "^mu_", 
     regex = TRUE)

# plot(
#   #conditional_effects(x = full_fitmevm, 
#                       spaghetti = TRUE, 
#                       ndraws = 2e2,
#                       effects = 'condition')
# )



# Plot NLvM-ME predictions against simulated data --------------------------


# . Collect fixed effects predictions -------------------------------------
#Get fixef predictions
sm_vm = summary(full_fitmevm, robust = TRUE)
prms_vm = with(sm_vm, rbind(fixed, spec_pars))
est_vm = data.frame(t(t(prms_vm)['Estimate',]))
#all draws for circular variables
full_fitmevm_mu_circ_draws = brms::as_draws_df(full_fitmevm,
                                                    variable = 'mu_circ') 
uw_mu_circ = unwrap_circular(full_fitmevm_mu_circ_draws$mu_circ)
full_fitmevm_mu_offs_draws = brms::as_draws_df(full_fitmevm,
                                                    variable = 'mu_offs') 
uw_mu_offset = unwrap_circular(full_fitmevm_mu_offs_draws$mu_offs)
est_vm$muangle_Intercept = median.circular(uw_mu_circ)
est_vm$muangle_condition = median.circular(uw_mu_offset)
# sapply(X = full_fitmevm_zmu_draws[1:n_indiv,],
#        FUN = hist,
#        breaks = 1e2)


# . Collect random effects predictions ------------------------------------
Cpal = colorRampPalette(colors = c(2:6,
                                   'seagreen',
                                   # 'salmon',
                                   # 'slategray3',
                                   'orange',
                                   'navajowhite4'
                                   ))
id_cols = sample(x = Cpal(n = n_indiv),
                 size = n_indiv,
                 replace = FALSE)
# [1] "#E59B14" "#DF536B" "#3EB0A2" "#BA22C0" "#24B0E5" "#626078" "#8A9630"
# [8] "#8B795E" "#7DB455" "#3ACAE0"

#Get raneff predictions for kappa
cf_vm = coef(full_fitmevm)$ID[,,'kappa_Intercept']
#We can collect linear scaled predictions for zmu
#these could potentially be wrong if the estimate is close to +-pi
rn = rownames(prms_vm)
rn_zmu = rn[grep(pattern = 'zmu',
                 x = rn)]
colnames(cf_vm) = colnames(prms_vm)[1:4]
#these are likely wrong for zmu
ran_prms_vm = rbind(cf_vm, prms_vm[rn_zmu,1:4])

#all draws
full_fitmevm_zmu_draws = brms::as_draws_df(full_fitmevm,
                                       variable = 'zmu_id') 
# sapply(X = full_fitmevm_zmu_draws[1:n_indiv,],
#        FUN = hist,
#        breaks = 1e2)
deg_pred = circular(x = deg(full_fitmevm_zmu_draws +
                              est_vm$muangle_Intercept),
                    type = 'angles',
                    unit = angle_unit,
                    template = 'geographics',
                    modulo = '2pi',
                    zero = pi/2,
                    rotation = angle_rot)
deg_pred = deg_pred[,1:n_indiv]#exclude non-estimated parameters

# with(subset(sim, condition %in% 0),
#      {
#        plot.circular(x = circular(x = angle, 
#                                   type = 'angles',
#                                   unit = angle_unit,
#                                   template = 'geographics',
#                                   modulo = '2pi',
#                                   zero = pi/2,
#                                   rotation = angle_rot
#        ),
#        stack = TRUE,
#        bins = 360/5,
#        sep = 0.5/dt_dim,
#        col = 'cyan4'
#        )
#      }
# )


# . plot circular random effects ------------------------------------------


par(pty = 's')#sometimes gets skipped? Needs to come first
par(mar =rep(0,4),
    mfrow = rep(x = ceiling(sqrt(n_indiv)), 
                times = 2) )
for(i in unique(indivs))
{
  with(subset(x = sim,
              subset = indiv %in% i &
                condition %in% 0      ),
       {
  plot.circular(x = circular(x = angle,
                               type = 'angles',
                               unit = angle_unit,
                               template = 'geographics',
                               modulo = '2pi',
                               zero = pi/2,
                               rotation = angle_rot
  ),
  sep = 2/dt_dim[1],
  stack = TRUE,
  bins = 360/5,
  col = adjustcolor(id_cols[i], alpha.f = 200/256)
  ) 
  }
  )
 par(new = T)
 plot.circular(x = circular(x = NULL, 
                           type = 'angles',
                           unit = 'degrees',
                           template = 'geographics',
                           modulo = '2pi',
                           zero = pi/2,
                           rotation = angle_rot
              ),
              shrink = 1.05,
              axes = F
              )
 points.circular(x = circular(x = deg_pred[,i],
                               type = 'angles',
                               unit = angle_unit,
                               template = 'geographics',
                               modulo = '2pi',
                               zero = pi/2,
                               rotation = angle_rot
                             ),
       sep = -1e-3,
       stack = TRUE,
       bins = 360,
       col = adjustcolor(id_cols[i], alpha.f = 1/256)
       ) 
  arrows.circular(x = median.circular(circular(x = mu0_sim[i],
                               type = 'angles',
                               unit = angle_unit,
                               template = 'geographics',
                               modulo = '2pi',
                               zero = pi/2,
                               rotation = angle_rot
                             )),
                  y = A1(exp(lk_all[i])),
                  length =0, 
                  lwd = 1,
                   col = adjustcolor(id_cols[i], alpha.f = 1)
       )  
  arrows.circular(x = median.circular(circular(x = deg_pred[,i],
                               type = 'angles',
                               unit = angle_unit,
                               template = 'geographics',
                               modulo = '2pi',
                               zero = pi/2,
                               rotation = angle_rot
                             )),
                  y = A1(inv_softplus(
                    est_vm$kappa_Intercept + cf_vm[i,'Estimate']
                    )
                    ),
                  lwd = 5,
                  length = 0,
                   col = adjustcolor(id_cols[i], alpha.f = 0.2)
       )
}

# . circular fixed effects ------------------------------------------------


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


arrows.circular(x = circular(mu_0,  
                             type = 'angles',
                             unit = 'degrees',
                             modulo = '2pi',
                             zero = pi/2,
                             template = 'geographics',
                             rotation = angle_rot),
                y = A1(kappa_pop),
                lwd = 1,
                length = 0,
                col = adjustcolor('cyan4', alpha.f = 1.0)
)
arrows.circular(x = circular(mu_0 + mu_offset[2],  
                             type = 'angles',
                             unit = 'degrees',
                             modulo = '2pi',
                             zero = pi/2,
                             template = 'geographics',
                             rotation = angle_rot),
                y = A1(exp(log(kappa_pop) + lk_offset[2])),
                lwd = 1,
                length = 0,
                col = adjustcolor('blue2', alpha.f = 1.0)
)

with(est_vm,
     {
       arrows.circular(x = circular(deg(mu_circ),  
                                    type = 'angles',
                                    unit = 'degrees',
                                    modulo = '2pi',
                                    zero = pi/2,
                                    template = 'geographics',
                                    rotation = angle_rot),
                       y = A1(inv_softplus(kappa_Intercept)),
                       lwd = 5,
                       length = 0,
                       col = adjustcolor('cyan2', alpha.f = 0.5)
       )
     }
)
with(est_vm,
     {
       arrows.circular(x = circular(deg(mu_circ + mu_offs),  
                                    type = 'angles',
                                    unit = 'degrees',
                                    modulo = '2pi',
                                    zero = pi/2,
                                    template = 'geographics',
                                    rotation = angle_rot),
                       y = A1(inv_softplus(kappa_Intercept + kappa_condition)),
                       lwd = 5,
                       length = 0,
                       col = adjustcolor('blue2', alpha.f = 0.5)
       )
     }
)

# . Model estimates -------------------------------------------------------

#all draws
full_fitmevm_fixed_draws = brms::as_draws_df(full_fitmevm,
                                           regex ='^b_') 
par(mar = c(0,3,4,0),
    bty = 'n')
with(full_fitmevm_fixed_draws,
     {
stripchart(x = round(inv_softplus(b_kappa_Intercept), digits = 2),
           pch = 20,
           col = gray(level = 0.5, alpha = 0.25),
           method = 'stack',
           vertical = TRUE,
           ylim = c(0, 5),
           xlim = c(1,2)+c(-1,1)*0.25,
           offset = 0.01,
           mgp = c(0,1,0.5),
           cex.main = 0.5, 
           bty = 'n')
 title(main = 'Population kappa',
       line = 0.1,
       cex.main = 0.5)
       abline(h = median(inv_softplus(b_kappa_Intercept)),
              col = adjustcolor(col = 'cyan4',alpha.f = 0.5),
              lwd = 5)
     }
)
abline(h = kappa_pop,
       col = 2,
       lwd = 2)

with(full_fitmevm_fixed_draws,
     {
stripchart(x = round(
  inv_softplus(b_kappa_Intercept + b_kappa_condition) - inv_softplus(b_kappa_Intercept), 
  digits = 2),
           pch = 20,
           col = gray(level = 0.5, alpha = 0.25),
           method = 'stack',
           vertical = TRUE,
           ylim = c(-1, 1)*2.5,
           xlim = c(1,2)+c(-1,1)*0.25,
  offset = 0.01,
  mgp = c(0,1,0.5),
  cex.main = 0.5, 
  bty = 'n')
       title(main = 'Condition effect on kappa',
             line = 0.1,
             cex.main = 0.5)
       abline(h = 0)
       abline(h = median( 
         inv_softplus(b_kappa_Intercept + b_kappa_condition) - inv_softplus(b_kappa_Intercept)
             ),
              col = adjustcolor(col = 'cyan4',alpha.f = 0.5),
              lwd = 5)
     }
)
abline(h = exp(log(kappa_pop) + lk_offset) - kappa_pop,
       col = 2,
       lwd = 2)

with(full_fitmevm_fixed_draws,
     {
stripchart(x = round(
  kappa_id, 
  digits = 2),
           pch = 20,
           col = gray(level = 0.5, alpha = 0.25),
           method = 'stack',
           vertical = TRUE,
            ylim = c(0, 25),
           xlim = c(1,2)+c(-1,1)*0.25,
  offset = 0.01,
  mgp = c(0,1,0.5),
  cex.main = 0.5, 
  bty = 'n')
       title(main = 'Random effects kappa',
             line = 0.1,
             cex.main = 0.5)
       abline(h = 0)
       abline(h = median( 
         kappa_id
             ),
              col = adjustcolor(col = 'cyan4',alpha.f = 0.5),
              lwd = 5)
     }
)
abline(h = exp(lk_indiv),
       col = 2,
       lwd = 2)

# . Add legend ------------------------------------------------------------
par(mar = rep(0,4),
    bty = 'o')
plot(x = NULL,
     xlim = c(-1,1),
     ylim = c(-1,1),
     axes = F,
     xlab = '',
     ylab = '')

legend(x = 'center',
       legend = c('true mean',
                  'simulated data',
                  'model estimates',
                  'model median'),
       pch = c(NA, 20, 19 , NA),
       lwd = c(2, 2, 2, 5),
       lty = c(1,NA,NA,1),
       col = c(gray(level = rep(0, 3),
                    alpha = c(1, 1, 0.5)),
               adjustcolor(col = 'cyan4',
                           alpha.f = 0.5)),
       cex = 1.2*sqrt(10/n_indiv)
)



# Compare model types -----------------------------------------------------


loo_lm = loo(full_fitlm)
loo_vm = loo(full_fit)
loo_vmme = loo(full_fitme)
loo_vmmevm = loo(full_fitmevm)
loo_compare(loo_lm, loo_vm, loo_vmme, loo_vmmevm)
# mu_0 = 275
# mu_offset = 5
#             elpd_diff se_diff
# full_fitmevm   0.0       0.0  
# full_fitme    -0.3       0.7  
# full_fit     -19.7       6.2  
# full_fitlm   -61.0      13.8  

# mu_0 = 175
# mu_offset = 20
#             elpd_diff se_diff
# full_fitmevm    0.0       0.0 
# full_fitme     -1.7       0.4 
# full_fit       -7.1       4.9 
# full_fitlm   -178.7      13.3 

# Hypothesis tests via model comparison -----------------------------------


# . Refit without effect of condition -------------------------------------


# . . Linear model version --------------------------------------------------
formula_lm_nocond = update(object = formula_lm,
                           .~. - condition)
prior_lm_nocond = subset(prior_lm,
                         subset = !(coef %in% 'condition') & 
                                   !(class %in% 'b')
                           )

# Full run
system.time(
  {
    full_fit_lm_nocond = brm( formula = formula_lm_nocond, # using our nonlinear formula
                      data = sim, # our data
                      prior = prior_lm_nocond, # our priors 
                      # stanvars = stan_var,#not needed
                      iter = 1000, # short run for 1000 iterations (less than 300 gives insufficient warmup time)
                      chains = 4, # 4 chains in parallel
                      cores = 4, # on 4 CPUs
                      init = 0,
                      refresh = 0, # don't echo chain progress
                      backend = 'cmdstanr') # use cmdstanr (other compilers broken)
  }
)

if(all_plots)
{
plot(full_fit_lm_nocond)
}

# . . Nonlinear vM with linear mixed effects --------------------------------



formula_nlvmme_nocond = formula_nlvmme
formula_nlvmme_nocond$pforms = within(formula_nlvmme_nocond$pforms,
                                      {
                                      muangle = update(muangle, .~. - condition)
                                      # kappa = '1 + (1|ID)'
                                      }) 

prior_nlvmme_nocond = subset(prior_nlvmme,
                         subset = !(coef %in% 'condition' & nlpar %in% 'muangle')
                           )

# Full run
system.time(
  {
    full_fit_nlvmme_nocond = brm( formula = formula_nlvmme_nocond, # using our nonlinear formula
                      data = sim, # our data
                      prior = prior_nlvmme_nocond, # our priors 
                      stanvars = stanvar(scode = mod_circular_fun, block = "functions") + 
                                  stanvar(scode = von_mises3_fun, block = "functions"),
                      iter = 1000, # short run for 1000 iterations (less than 300 gives insufficient warmup time)
                      chains = 4, # 4 chains in parallel
                      cores = 4, # on 4 CPUs
                      init = 0,
                      refresh = 0, # don't echo chain progress
                      backend = 'cmdstanr') # use cmdstanr (other compilers broken)
  }
)

if(all_plots)
{
plot(full_fit_nlvmme_nocond)
}
# . . Nonlinear vM with von Mises mixed effects --------------------------------



formula_nlvmmevm_nocond = formula_nlvmmevm
formula_nlvmmevm_nocond$pforms = within(formula_nlvmmevm_nocond$pforms,
                                      {
                                      muangle = update(muangle, .~. - condition)
                                      # kappa = '1 + (1|ID)'
                                      }) 

prior_nlvmmevm_nocond = subset(prior_nlvmmevm,
                         subset = !(coef %in% 'condition' & nlpar %in% 'muangle')
                           )

# Full run
system.time(
  {
    full_fit_nlvmmevm_nocond = brm( formula = formula_nlvmmevm_nocond, # using our nonlinear formula
                      data = sim, # our data
                      prior = prior_nlvmmevm_nocond, # our priors 
                      stanvars = stanvar(scode = mod_circular_fun, block = "functions") + 
                                  stanvar(scode = von_mises3_fun, block = "functions") + 
                                  zkappa_var + zmu_var,
                      iter = 1000, # short run for 1000 iterations (less than 300 gives insufficient warmup time)
                      chains = 4, # 4 chains in parallel
                      cores = 4, # on 4 CPUs
                      init = 0,
                      refresh = 0, # don't echo chain progress
                      backend = 'cmdstanr') # use cmdstanr (other compilers broken)
  }
)
if(all_plots)
{
plot(full_fit_nlvmmevm_nocond)
}

# . Model comparison ------------------------------------------------------


# . . Linear model version --------------------------------------------------
loo_compare(loo(full_fitlm), loo(full_fit_lm_nocond))
# > mu_offset[2]
# [1] 5°
#                   elpd_diff se_diff
# full_fit_lm_nocond  0.0       0.0   #no effect model chosen
# full_fitlm         -0.5       0.9   #no clear support for condition difference
# > mu_offset[2]
# [1] 10
#                   elpd_diff se_diff
# full_fit_lm_nocond  0.0       0.0   #no effect model chosen
# full_fitlm         -1.1       0.8   #clear support for no condition difference
# > mu_offset[2]
# [1] 15
#                   elpd_diff se_diff
# full_fit_lm_nocond  0.0       0.0   
# full_fitlm         -0.7       0.3   
# > mu_offset[2]
# [1] 20
#                   elpd_diff se_diff # pareto_k > 0.7 in model 'full_fitlm'
# full_fit_lm_nocond  0.0       0.0   #no effect model chosen
# full_fitlm         -1.3       0.5   #clear support for no condition difference

# . . Nonlinear vM with linear mixed effects --------------------------------
loo_compare(loo(full_fitme), loo(full_fit_nlvmme_nocond))
# > mu_offset[2]
# [1] 5°
#                       elpd_diff se_diff
# full_fit_nlvmme_nocond  0.0       0.0   #no effect model chosen
# full_fitme             -0.8       0.9   #no clear support for condition difference
# > mu_offset[2]
# [1] 10
#                         elpd_diff se_diff
# full_fitme              0.0       0.0   #effect model chosen
# full_fit_nlvmme_nocond -0.2       1.7   #no clear support for condition difference
# > mu_offset[2]
# [1] 15
#                       elpd_diff se_diff
# full_fitme              0.0       0.0  #effect model chosen 
# full_fit_nlvmme_nocond -2.1       2.4   #support is not completely clear
#             Estimate Est.Error l-95% CI u-95% CI Rhat
# mu_offs     0.28      0.11     0.05     0.51 1.00 # summary does support effect
# > mu_offset[2]
# [1] 20
#                       elpd_diff se_diff
# full_fitme              0.0       0.0   #effect model chosen
# full_fit_nlvmme_nocond -4.6       3.1   #clear support for condition difference
#             Estimate Est.Error l-95% CI u-95% CI Rhat
# mu_offs     0.36      0.12     0.13     0.59 1.00

# . . Nonlinear vM with von Mises mixed effects --------------------------------
loo_compare(loo(full_fitmevm), loo(full_fit_nlvmmevm_nocond))
# > mu_offset[2]
# [1] 5°
#                         elpd_diff se_diff
# full_fit_nlvmmevm_nocond  0.0       0.0   #no effect model chosen
# full_fitmevm             -0.5       0.8   #no clear support for condition difference
# > mu_offset[2]
# [1] 10
#                           elpd_diff se_diff
# full_fitmevm              0.0       0.0   #effect model chosen
# full_fit_nlvmmevm_nocond -0.6       1.6   #no clear support for condition difference                               
#             Estimate Est.Error l-95% CI u-95% CI Rhat
# mu_offs      0.17      0.11    -0.04     0.39 1.00
# > mu_offset[2]
# [1] 15
#                         elpd_diff se_diff
# full_fitmevm              0.0       0.0   #effect model chosen
# full_fit_nlvmmevm_nocond -1.9       2.2   #support not completely clear
#             Estimate Est.Error l-95% CI u-95% CI Rhat
# mu_offs      0.26      0.12     0.01     0.50 1.01 # summary does support effect
# > mu_offset[2]
# [1] 20
#                         elpd_diff se_diff
# full_fitmevm              0.0       0.0   #effect model chosen
# full_fit_nlvmmevm_nocond -4.4       2.9   #clear support for condition difference
#             Estimate Est.Error l-95% CI u-95% CI Rhat
# mu_offs      0.34      0.11     0.11     0.57 1.00

# Test w/ real data -------------------------------------------------------


# . Select file ---------------------------------------------------------

# set path to file
if(sys_win){#choose.files is only available on Windows
  message('\n\nPlease select the ".csv" file\n\n')
  Sys.sleep(0.5)#goes too fast for the user to see the message on some computers
  path_file  <- choose.files(
    default = file.path(ltp,'Documents', "*.csv"),#For some reason this is not possible in the "root" user
    caption = 'Please select the ".csv" file'
  )
}else{
  message('\n\nPlease select the ".csv" file\n\n')
  Sys.sleep(0.5)#goes too fast for the user to see the message on some computers
  path_file <- file.choose(new=F)
}
#show the user the path they have selected
if(is.null(path_file) | !length(path_file))
{stop('No file selected.')}else
{print(path_file)}


# Read in the data and format ---------------------------------------------

lp = read.table(file = path_file, 
                header = T, 
                sep  = '\t')
View(lp)

lp = within(lp,
            {
              ID = as.factor(Beetle) # beedance identifier as a factor
              trial = as.factor(Trial) # date as a factor
              condition = ifelse(test = Experiment %in% 'LPexper.20191117.SecurityLightBergsigON',
                                           yes = 1,
                                           no = 0) # light ON as a binomial
              signed_angle = Mod360.180(Heading)  # bearing between -180 and 180
              rad_angle = circular(rad(signed_angle),# bearing between -pi and pi
                               rotation = 'clock') # circular format suppresses later warnings
            }
)


# . Fit the normal raneff model -------------------------------------------
# Full run
system.time(
  {
    full_fit_LP = brm( formula = formula_nlvmme, # using our nonlinear formula
                      data = lp, # our data
                      prior = prior_nlvmme, # our priors 
                      stanvars = stan_var,
                      iter = 1000, # short run for 1000 iterations (less than 300 gives insufficient warmup time)
                      chains = 4, # 4 chains in parallel
                      cores = 4, # on 4 CPUs
                      init = 0,
                      refresh = 0, # don't echo chain progress
                      backend = 'cmdstanr') # use cmdstanr (other compilers broken)
  }
)
fx_names = names(fixef(full_fit_LP)[,1])
nm_angle_vars = fx_names[grepl(pattern = 'angle', x = fx_names)]
nm_kappa_vars = fx_names[grepl(pattern = 'kappa', x = fx_names)]
A1sp = function(x){A1(inv_softplus(x))}
trans_lst = c(replicate(n = length(nm_angle_vars),
                      expr = mod_circular,
                      simplify = FALSE),
              replicate(n = length(nm_kappa_vars),
                      expr = A1sp,
                      simplify = FALSE) )
              
names(trans_lst) = paste0('b_', c(nm_angle_vars, nm_kappa_vars))



plot(full_fit_LP,
     nvariables = 10,
     variable = "^b_", 
     transformations = trans_lst,
     regex = TRUE)

NZmod_circular = function(x)
{
  mn = mean.circular(x)
  return(mod_circular(x - mn))
}
nztrans_lst = replicate(n = length(trans_lst),
                        expr = NZmod_circular, 
                        simplify = FALSE)
names(nztrans_lst) = names(trans_lst)
plot(full_fit_LP,
     variable = "^b_mu", 
     transformations = nztrans_lst,
     regex = TRUE)


plot(full_fit_LP,
     nvariables = 10,
     variable = "^sd_", 
     regex = TRUE)

summary(full_fit_LP, robust = TRUE)

# plot(
#   #conditional_effects(x = full_fit_LP, 
#                       spaghetti = TRUE, 
#                       ndraws = 2e2,
#                       effects = 'condition')
# )





