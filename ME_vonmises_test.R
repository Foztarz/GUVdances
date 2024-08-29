#FOR A 'CLEAN' RUN, PRESS ctrl+shift+F10 to RESTART Rstudio
graphics.off()
# Details ---------------------------------------------------------------
#       AUTHOR:	James Foster              DATE: 2024 08 28
#     MODIFIED:	James Foster              DATE: 2024 08 28
#
#  DESCRIPTION: Fit hierarchical maximum-likelihood von Mises.
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
#- Test for angles close to 180
#- Compare circular and Gaussian raneff
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

#  .  User input -----------------------------------------------------------
set.seed(20240708)#day the script was fixed
n_iter = 1e4#4 #optimisation iterations

paired_data = TRUE # Are the data in the two columns paired (each from the same animal or group)?
csv_sep = ','#Is the csv comma separated or semicolon separated? For tab sep, use "\t"
angle_name = "angles" #The title of the column with angles; NO SPACES PLEASE
angle_unit = "degrees" # "degrees" or "radians"
angle_rot = "counter" # "clock" or "counter"

#Check the operating system and assign a logical flag (T or F)
sys_win <- Sys.info()[['sysname']] == 'Windows'
#User profile instead of home directory
if(sys_win){
  #get rid of all the backslashes
  ltp <- gsub('\\\\', '/', Sys.getenv('USERPROFILE'))#Why does windows have to make this so difficult
}else{#Root directory should be the "HOME" directory on a Mac (or Linux?)
  ltp <- Sys.getenv('HOME')#Life was easier on Mac
}

# Simulate data  ------------------------------------------------
n_angles = 44
# minimum discriminable angle appears to be approx 35°
kappa_both = A1inv(0.9);A1inv(0.7) #concentration around each trial mean
logkappa_var = 1.0 #scale of random variation in concentration (log units)
mu_1 = rad(185);rcircularuniform(n = 1)#
mu_offset = rad(30);rcircularuniform(n = 1)
if(paired_data)
{
kappa_indiv = A1inv(0.98) #concentration across individuals (pairs)
#mean angle in trail 1 for each individual (pair)
mu1_sim = rvonmises(n = n_angles,
                      mu = circular(mu_1, units = "radians", zero = 0, rotation = "counter"),#random angle
                      kappa = kappa_indiv#the wider the distribution of individual biases, the greater the influence of pairing
                      )
#simulate the full dataset
sim = data.frame(
                 angle_1 = round(c(suppressWarnings( #rvonmises converts to circular and warns
                   sapply(X = circular(mu1_sim, units = "radians", zero = 0, rotation = "counter"),
                          FUN = rvonmises,
                          n = 1,
                           kappa = exp(log(kappa_both) + rnorm(n = 1, sd = logkappa_var))
                   )
                 ))*180/pi),#convert to angles and round to nearest degree
                 angle_2 = round(c(suppressWarnings( #rvonmises converts to circular and warns
                             sapply(X =mu1_sim + mu_offset,# true difference,
                                    FUN = rvonmises,
                                    n = 1,
                                    kappa = exp(log(kappa_both) + rnorm(n = 1, sd = logkappa_var))
                             )
                 ))*180/pi) #convert to angles and round to nearest degree
                 )
}else
{
n_angles2 = ceiling(0.75*n_angles)
mu1_sim = rcircularuniform(n = 1,control.circular = list(units = angle_unit))
sim = data.frame(
                 angle_1 = round(c(suppressWarnings( #rvonmises converts to circular and warns
                   rvonmises(n = n_angles,
                             mu = circular(x = mu1_sim, units = angle_unit),
                           kappa = exp(log(kappa_both) + rnorm(n = 1, sd = logkappa_var))
                   )
                 ))),#convert to angles and round to nearest degree
                 angle_2 = round(c(suppressWarnings( #rvonmises converts to circular and warns
                   rvonmises(n = n_angles2,
                             mu = circular(x = mu1_sim+deg(mu_offset), units = angle_unit),
                             kappa = exp(log(kappa_both) + rnorm(n = 1, sd = logkappa_var))
                   )
                 ),
                 circular(x = rep(x = NA, times = n_angles - n_angles2),
                         units = angle_unit) #convert to angles and round to nearest degree
                 ) )
                )
}
# #save somewhere the user likely keeps data
# write.table(x = sim,
#             file = file.path(ltp,'Documents', "simulated_angles.csv"),
#             sep = csv_sep,
#             row.names = FALSE
#             )
adata = sim
dt_dim = dim(adata)
# Plot simulated data -----------------------------------------------------

par(mar =rep(0,4))
plot.circular(x = circular(x = adata$angle_1, 
                           type = 'angles',
                           unit = angle_unit,
                           template = 'geographics',
                           modulo = '2pi',
                           zero = pi/2,
                           rotation = angle_rot
),
stack = TRUE,
bins = 360/5,
sep = 0.5/dt_dim[1],
col = 'cyan4'
)
par(new = T)
plot.circular(x = circular(x = adata$angle_2, 
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

# Organise data ----------------------------------------------------------
IDs = paste0(sort(rep(LETTERS, times = 26)), letters)
longdata = with(adata,
                data.frame(angle = c(angle_1, angle_2),
                           ID = rep(IDs[1:dt_dim[1]], times = 2),
                           condition = sort(rep(0:1, times = dt_dim[1]))
                           )
                )

longdata = within(longdata,
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
          data = longdata)

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
#generate 
mu_gen <- "
real mu_circ = mod_circular(b_muangle[1]);
real mu_offs = mod_circular(b_muangle[2]);
"
  
stan_var = stanvar(scode = mod_circular_fun, block = "functions") + 
              stanvar(scode = mu_gen, block = "genquant")

# . Short dummy run to check the influence of the priors ------------------


#double check that the prior distribution is viable by first setting up a short dummy run
# Dummy run
system.time(
  {
    dummy_fit = brm( formula = formula_nlvm, # using our nonlinear formula
                     data = longdata, # our data
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
plot(dummy_fit,
     variable = '^mu_',
     regex = TRUE)
plot(dummy_fit,
     variable = '^kappa',
     regex = TRUE)

conditional_effects(dummy_fit)


# . To check that estimation works  ---------------------------------------
# Run the model for a small number of iterations to check that it is possible to 
# estimate all parameters. This may be the point where we encounter numerical errors
# if the formula or priors are misspecified.
# e.g. if the formula returns estimates of correct choice rate outside of [0,1].

# Short run
system.time(
  {
    short_fit = brm( formula = formula_nlvm, # using our nonlinear formula
                     data = longdata, # our data
                     prior = prior_nlvm, # our priors 
                     stanvars = stan_var,
                     iter = 300, # short run for 300 iterations (less than 300 gives insufficient warmup time)
                     chains = 4, # 4 chains in parallel
                     cores = 4, # on 4 CPUs
                     refresh = 0, # don't echo chain progress
                     backend = 'cmdstanr') # use cmdstanr (other compilers broken)
  }
)
plot(short_fit)
plot(short_fit,
     variable = '^mu_',
     regex = TRUE)
summary(short_fit)

plot(
  conditional_effects(x = short_fit, 
                      spaghetti = TRUE, 
                      ndraws = 2e2,
                      effects = 'condition')
)

# . Full runs -------------------------------------------------------------


# Full run
system.time(
  {
    full_fit = brm( formula = formula_nlvm, # using our nonlinear formula
                     data = longdata, # our data
                     prior = prior_nlvm, # our priors 
                     stanvars = stan_var,
                     iter = 1000, # short run for 1000 iterations (less than 300 gives insufficient warmup time)
                     chains = 4, # 4 chains in parallel
                     cores = 4, # on 4 CPUs
                     refresh = 0, # don't echo chain progress
                     backend = 'cmdstanr') # use cmdstanr (other compilers broken)
  }
)

plot(full_fit,
     variable = '^mu_',
     regex = TRUE)
plot(full_fit)
summary(full_fit)

plot(
  conditional_effects(x = short_fit, 
                      spaghetti = TRUE, 
                      ndraws = 2e2,
                      effects = 'condition')
)


# Linear model sanity check -----------------------------------------------
#set up model fit
formula_lm = bf(
  #set up a formula for the curve as a whole,
  formula = rad_angle ~ condition + (1|ID),
  nl = FALSE)#the joint distribution for these parameters is undefined, and therefore the parameters themselves are "nonlinear"

prior_lm = get_prior(formula = formula_lm,
                         data = longdata)

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
                       data = longdata, # our data
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
                       data = longdata, # our data
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
                       data = longdata, # our data
                       prior = prior_lm, # our priors 
                       iter = 1000, # full run for 1000 iterations
                       chains = 4, # 4 chains in parallel
                       cores = 4, # on 4 CPUs
                       refresh = 0, # don't echo chain progress
                       backend = 'cmdstanr') # use cmdstanr (other compilers broken)
  }
)

# . Plot LM with simulated data -----------------------------------------------------

par(mar =rep(0,4))
plot.circular(x = circular(x = adata$angle_1, 
                           type = 'angles',
                           unit = angle_unit,
                           template = 'geographics',
                           modulo = '2pi',
                           zero = pi/2,
                           rotation = angle_rot
),
stack = TRUE,
bins = 360/5,
sep = 0.5/dt_dim[1],
col = 'cyan4'
)
par(new = T)
plot.circular(x = circular(x = adata$angle_2, 
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


arrows.circular(x = circular(deg(mu_1),  
                             type = 'angles',
                             unit = 'degrees',
                             modulo = '2pi',
                             zero = pi/2,
                             template = 'geographics',
                             rotation = angle_rot),
                y = A1(kappa_both),
                lwd = 1,
                length = 0,
                col = adjustcolor('cyan4', alpha.f = 1.0)
)
arrows.circular(x = circular(deg(mu_1 + mu_offset),  
                             type = 'angles',
                             unit = 'degrees',
                             modulo = '2pi',
                             zero = pi/2,
                             template = 'geographics',
                             rotation = angle_rot),
                y = A1(kappa_both),
                lwd = 1,
                length = 0,
                col = adjustcolor('blue2', alpha.f = 1.0)
)

sm_lm = summary(full_fitlm, robust = TRUE)
prms_lm = with(sm_lm, rbind(fixed, spec_pars))
est_lm = data.frame(t(t(prms_lm)['Estimate',]))
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



# Nonlinear ME Stan version --------------------------------------------------


#set up model fit
formula_nlvmme = bf(
  #set up a formula for the curve as a whole,
  formula = rad_angle ~ mod_circular(muangle),
  muangle ~  condition + (1|ID), #N.B. this is similar to the slope, so all of its effects depend on stimulus level
  kappa ~ condition + (1|ID), #N.B. this is similar to the slope, so all of its effects depend on stimulus level
  family = von_mises(link = "identity",
                     link_kappa = 'softplus'),
  nl = TRUE)#the joint distribution for these parameters is undefined, and therefore the parameters themselves are "nonlinear"

prior_nlvmme = get_prior(formula = formula_nlvmme,
          data = longdata)

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
                     data = longdata, # our data
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

plot(dummy_fitme,
     variable = '^mu_',
     regex = TRUE)
plot(dummy_fitme)

conditional_effects(dummy_fit, spaghetti = TRUE)


# . To check that estimation works  ---------------------------------------
# Run the model for a small number of iterations to check that it is possible to 
# estimate all parameters. This may be the point where we encounter numerical errors
# if the formula or priors are misspecified.
# e.g. if the formula returns estimates of correct choice rate outside of [0,1].

# Short run
system.time(
  {
    short_fitme = brm( formula = formula_nlvmme, # using our nonlinear formula
                     data = longdata, # our data
                     prior = prior_nlvmme, # our priors 
                     stanvars = stan_var,
                     iter = 300, # short run for 300 iterations (less than 300 gives insufficient warmup time)
                     chains = 4, # 4 chains in parallel
                     cores = 4, # on 4 CPUs
                     refresh = 0, # don't echo chain progress
                     backend = 'cmdstanr') # use cmdstanr (other compilers broken)
  }
)

plot(short_fitme,
     variable = '^mu_',
     regex = TRUE)
plot(short_fitme, 
     nvariables = 10,
     variable = "^b_", 
     regex = TRUE)
summary(short_fitme)

plot(
  conditional_effects(x = short_fit, 
                      spaghetti = TRUE, 
                      ndraws = 2e2,
                      effects = 'condition')
)

# . Full runs -------------------------------------------------------------


# Full run
system.time(
  {
    full_fitme = brm( formula = formula_nlvmme, # using our nonlinear formula
                     data = longdata, # our data
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

plot(full_fitme,
     variable = "^mu_", 
     regex = TRUE)

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
rhat(full_fitme,
     pars = "^mu_", 
     regex = TRUE)

plot(
  conditional_effects(x = short_fit, 
                      spaghetti = TRUE, 
                      ndraws = 2e2,
                      effects = 'condition')
)


# Plot simulated data -----------------------------------------------------

par(mar =rep(0,4))
plot.circular(x = circular(x = adata$angle_1, 
                           type = 'angles',
                           unit = angle_unit,
                           template = 'geographics',
                           modulo = '2pi',
                           zero = pi/2,
                           rotation = angle_rot
),
stack = TRUE,
bins = 360/5,
sep = 0.5/dt_dim[1],
col = 'cyan4'
)
par(new = T)
plot.circular(x = circular(x = adata$angle_2, 
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


arrows.circular(x = circular(deg(mu_1),  
                            type = 'angles',
                            unit = 'degrees',
                            modulo = '2pi',
                            zero = pi/2,
                            template = 'geographics',
                            rotation = angle_rot),
               y = A1(kappa_both),
               lwd = 1,
               length = 0,
               col = adjustcolor('cyan4', alpha.f = 1.0)
)
arrows.circular(x = circular(deg(mu_1 + mu_offset),  
                            type = 'angles',
                            unit = 'degrees',
                            modulo = '2pi',
                            zero = pi/2,
                            template = 'geographics',
                            rotation = angle_rot),
               y = A1(kappa_both),
               lwd = 1,
               length = 0,
               col = adjustcolor('blue2', alpha.f = 1.0)
)

sm = summary(full_fitme, robust = TRUE)
prms = with(sm, rbind(fixed, spec_pars))
est = data.frame(t(t(prms)['Estimate',]))
with(est,
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
with(est,
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

loo_lm = loo(full_fitlm)
loo_vm = loo(full_fit)
loo_vmme = loo(full_fitme)
loo_compare(loo_lm, loo_vm, loo_vmme)

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

plot(
  conditional_effects(x = full_fit_LP, 
                      spaghetti = TRUE, 
                      ndraws = 2e2,
                      effects = 'condition')
)





