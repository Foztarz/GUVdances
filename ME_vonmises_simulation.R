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
#- Histograms in descriptive plots  +
#- Choose appropriate prior for zkappa
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

#degree version
unwrap_circular_deg = function(x)
{
  mux = mean.circular(x = circular(x = x, template = 'none'))
  centx = atan2(y = sin(x - mux),
                x = cos(x  - mux))
  unwrx = centx + mux
  return(deg(unwrx))
}

#invert the softplus link
#https://en.wikipedia.org/wiki/Softplus
#TODO confusing name, this is softplussing, inverse would be log(exp(x)-1)
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

#histograms on a vertical axis

VertHist = function(data,
                    breaks = 1e2,
                    ylab = 'data',
                    xlab = 'density',
                    ylim = NULL,
                    main = '',
                    col = 'gray',
                    border = NA,
                    ...)
{
  hst = hist(x = data,
             breaks = breaks,
             plot = FALSE)
  with(hst,
       {
         plot(x = NULL,
              xlim = c(0, max(density)),
              ylim = if(is.null(ylim)){range(mids)}else{ylim},
              xlab = xlab,
              ylab = ylab,
              main = main)
         for(i in 1:length(mids))
         {
           rect(xleft = 0,
                xright = density[i],
                ybottom = breaks[i], 
                ytop = breaks[i + 1],
                col = col,
                border = border,
                ...
           )
         }
       }
  )
}

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
  real mux = mean_circular(y, N);
  vector[N] centx = y - mux;
  for(n in 1:N)
  {
  centx[n] = mod_circular(centx[n]);
  }
  vector[N] unwrx = centx + mux;
  return(unwrx);
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
real mu_circ = mod_circular(b_fmu[1]);
real mu_offs = mod_circular(b_fmu[2]);
",
                 block = 'genquant')
# #or maybe?
# mu_gen = stanvar(scode = "
# vector[size(b_fmu)] mu_circ; //modulo circular estimate
# for (i in 1:size(b_fmu)){
# mu_circ[i] = mod_circular(b_fmu[i])
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



stanvars_intercepts = stan_mvm_fun + mu_gen + zkappa_var + zmu_var
stanvars_slopes = stan_mvm_fun + mu_gen + zkappa_var + zkappa_var_slope + 
                  zmu_var_slope #includes intercepts

# Input Variables ----------------------------------------------------------

all_plots = FALSE # to speed up
#  .  User input -----------------------------------------------------------
set.seed(20240905)#day the simulation was updated
n_iter = 1e4 #optimisation iterations

csv_sep = ','#Is the csv comma separated or semicolon separated? For tab sep, use "\t"
angle_name = "angles" #The title of the column with angles; NO SPACES PLEASE
angle_unit = "degrees" # "degrees" or "radians"
angle_rot = "clock" # "clock" or "counter"

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
mu_0 = circular(x = 178,
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


# Intercepts & slopes model --------------------------------------------------


## Formula ---------------------------------------------------------------

#set up model fit
formula_int_slope = bf(
  formula = rad_angle ~ #set up a formula for the mean angle, modulus to (-pi,pi)
    mod_circular(fmu + zmu), # mean angle combines fixed and random effects
  fmu ~  condition, # fixed effects only change as a function of condition
  zmu ~  0+ID+ID:condition, #random effects change as a function of individual and condition and their interaction
  kappa ~ condition + (1+condition|ID), #for kappa this occurs in linear space, and BRMS can set it up automatically
  family = von_mises(link = "identity", # the mean angle will be returned as-is
                     link_kappa = 'softplus'),#kappa will be returned via the softplus link https://en.wikipedia.org/wiki/Softplus
  nl = TRUE)#to accept user-defined extra parameters (zmu) we need to treat the formula as nonlinear



## Priors ----------------------------------------------------------------
prior_int_slope = get_prior(formula = formula_int_slope,
                               data = sim,
                            check = FALSE)
#suggests 113 possible priors
#Not cooperating w/ BRMS vectorisation
# prior_int_slope = prior(prior = 'von_mises3(0, 0.1)',
#                         nlpar = 'fmu',
#                         coef = ''
#                         ) + 
#                   prior(prior = 'von_mises3(0, 1.0)',
#                         nlpar = 'fmu',
#                         coef = 'ID'
#                         ) + 
#                   prior(prior = 'von_mises3(0, log1p_exp(zkappa))',
#                         nlpar = 'zmu',
#                         coef = 'ID:condition'
#                         ) + 
#                   prior(prior = 'von_mises3(0, log1p_exp(zkappa_condition))',
#                         nlpar = 'zmu',
#                         coef = 'condition'
#                         ) + 
#                   prior(prior = 'normal(1.5, 5.0)',
#                         dpar = 'kappa',
#                         class = 'Intercept') + 
#                   prior(prior = 'normal(0.0, 5.0)',
#                         dpar = 'kappa',
#                         coef = 'condition') + 
#                   prior(prior = 'student_t(3.0, 0.0, 5.0)',
#                         dpar = 'kappa',
#                         coef = 'Intercept',
#                         class = 'sd',
#                         group = 'ID') + 
#                   prior(prior = 'student_t(0.0, 0.0, 5.0)',
#                         dpar = 'kappa',
#                         coef = 'condition',
#                         class = 'sd',
#                         group  = 'ID')

### assign BRMS default priors -------------------------------------------
prior_int_slope = within(prior_int_slope,
                            {
#fixed effects on mean angle are von Mises distributed (von_mises3 converts estimates to modulo (-pi,pi))                              
      prior[nlpar %in% 'fmu' & coef %in% 'Intercept'] = 'von_mises3(0, 0.1)'# very weak bias to zero (could be no bias?)
      prior[nlpar %in% 'fmu' & class %in% 'b'] = 'von_mises3(0, 1.0)'#moderate bias to zero, no effect
#random effects on mean angle are von Mises distributed, with a kappa parameter estimated from the data
      prior[nlpar %in% 'zmu' & coef %in% 'Intercept'] = 'von_mises3(0, log1p_exp(zkappa))'
      prior[nlpar %in% 'zmu' & class %in% 'b'] = 'von_mises3(0, log1p_exp(zkappa))'
      prior[nlpar %in% 'zmu' & class %in% 'b' 
            & grepl(pattern = 'condition', #the random effect of condition has a different kappa
                    x = coef)] = 'von_mises3(0, log1p_exp(zkappa_condition))'
#fixed effects on kappa are normally distributed on the softplus scale
      prior[dpar %in% 'kappa' & class %in% 'Intercept'] = 'normal(1.5, 5.0)'#weak expectation of kappa around 1.7 (mean vector around 0.65)
      prior[dpar %in% 'kappa' & class %in% 'b'] = 'normal(0.0, 5.0)'#weak expectation of condition effect around 0
#random effects on kappa are t-distributed on the softplus scale
      prior[dpar %in% 'kappa' & class %in% 'sd'] = 'student_t(3, 0, 5.0)' #wide prior 
                            }
)


### add extra priors for the random effects mean angles -----------------



prior_int_slope = prior_int_slope + #random effects kappas are t-distributed on a softplus scale
  set_prior("target += student_t_lpdf(zkappa | 3, 25, 5)", #expect high concentration (low variation) 
            check = FALSE)+
  set_prior("target += student_t_lpdf(zkappa_condition | 3, 25, 5)", #expect high concentration (low variation) 
            check = FALSE)

#is this a good prior for zkappa?
xx = seq(from  = -50,
         to  = 50,
         length.out = 1e3)
if(all_plots)
{
  zk_prior = c(3,25,5)
  #the t-distributed prior appears to bias towards low values
  plot(x = xx,
       y = brms::dstudent_t(xx,
                            df = zk_prior[1],
                            mu = zk_prior[2],
                            sigma = zk_prior[3]),
       type = 'l',
       xlab = 'softplus kappa_id',
       ylab = 'probability density',
       main = 'zkappa prior, student_t(3,5,2.5)',
       lwd = 2,
       col = 6,
       xlim = c(-50,50),
       ylim = c(0,0.2)
  )
  #plot equidistant mean vectors for reference
  abline(v = log(exp(A1inv(1:50 / 50)) - 1),
         col = gray(0.5, alpha = 0.5))
  #the t-distributed prior appears to bias towards low values
  plot(x = inv_softplus(xx),
       y = brms::dstudent_t(xx,
                            df = zk_prior[1],
                            mu = zk_prior[2],
                            sigma = zk_prior[3]),
       type = 'l',
       xlab = 'kappa_id',
       ylab = 'probability density',
       main = 'zkappa prior, student_t(3,5,2.5)',
       lwd = 2,
       col = 5,
       xlim = c(0,50),
       ylim = c(0,0.2)
  )
  abline(v = A1inv(1:50 / 50),
         col = gray(0.5, alpha = 0.5))
}

## Short dummy run to check the influence of the priors ------------------


#double check that the prior distribution is viable by first setting up a short dummy run
# Dummy run
#TODO work out why this samples less efficiently than with data
system.time( #currently takes about 2 minutes for 2000 iterations
  {
    dummy_int_slope = brm( formula = formula_int_slope, # using our nonlinear formula
                               data = sim, # our data
                               prior = prior_int_slope, # our priors 
                               stanvars = stanvars_slopes,
                               sample_prior = 'only', #ignore the data to check the influence of the priors
                               iter = 10000, # can only estimate with enough iterations for params
                               chains = 4, # 4 chains in parallel
                               cores = 4, # on 4 CPUs
                               refresh = 0, # don't echo chain progress
                               backend = 'cmdstanr') # use cmdstanr (other compilers broken)
  }
)

if(all_plots)
{
  plot(dummy_int_slope,
       variable = 'fmu',
       regex = TRUE,
       transform = unwrap_circular_deg)
  plot(dummy_int_slope,
       variable = '^kappa_id',
       regex = TRUE)
  #samples inefficiently?
  plot(dummy_int_slope,
       variable = '^zkappa',
       regex = TRUE)
  plot(dummy_int_slope,
       variable = 'zmu_id',
       transform = unwrap_circular_deg,
       nvariables = 5,
       ask = FALSE)
  plot(dummy_int_slope,
       variable = 'zmu_id_condition',
       transform = unwrap_circular_deg,
       nvariables = 5,
       ask = FALSE)
  plot(dummy_int_slope)
}


## Full run --------------------------------------------------------------

# full run
system.time(#takes less than 1 minute
  {
    full_int_slope = brm( formula = formula_int_slope, # using our nonlinear formula
                             data = sim, # our data
                             prior = prior_int_slope, # our priors 
                             stanvars = stanvars_slopes,
                             iter = 1000, #doesn't take a lot of runs
                             chains = 4, # 4 chains in parallel
                             cores = 4, # on 4 CPUs
                             refresh = 0, # don't echo chain progress
                             backend = 'cmdstanr') # use cmdstanr (other compilers broken)
  }
)
if(all_plots)
{
  plot(full_int_slope,
       variable = 'fmu',
       regex = TRUE,
       transform = unwrap_circular_deg)
  plot(full_int_slope,
       variable = '^kappa_id',
       regex = TRUE)
  plot(full_int_slope,
       variable = '^zkappa',
       regex = TRUE)
  plot(full_int_slope,
       variable = 'zmu_id',
       transform = unwrap_circular_deg,
       nvariables = 5,
       ask = FALSE)
  plot(full_int_slope,
       variable = 'zmu_id_condition',
       transform = unwrap_circular_deg,
       nvariables = 5,
       ask = FALSE)
  plot(full_int_slope)
}


# Plot NLvM-ME predictions against simulated data --------------------------


## Collect fixed effects predictions -------------------------------------
#Get fixef predictions
sm_vm = summary(full_int_slope, robust = TRUE)
prms_vm = with(sm_vm,
               rbind(fixed, #the fixed effects
                     spec_pars) #generated parameters 
               )
est_vm = data.frame(t(t(prms_vm)['Estimate',])) # extract just the estimate
#all draws for circular variables
#circular intercept
full_int_slope_mu_circ_draws = brms::as_draws_df(full_int_slope,
                                               variable = 'mu_circ') 
uw_mu_circ = unwrap_circular(full_int_slope_mu_circ_draws$mu_circ) # unwrap
#circular condition offset
full_int_slope_mu_offs_draws = brms::as_draws_df(full_int_slope,
                                               variable = 'mu_offs') 
uw_mu_offset = unwrap_circular(full_int_slope_mu_offs_draws$mu_offs)
est_vm$fmu_Intercept = median.circular(uw_mu_circ)
est_vm$fmu_condition = median.circular(uw_mu_offset)


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

# [1] "#2297E6" "#69923E" "#EE9E0D" "#3DAFA5" "#3FC3DF" "#57A4D9" "#BC8B35" "#9F2F9F" "#72BE53" "#B529C1"
# [11] "#DF536B" "#CD6467" "#61D04F" "#5B6673" "#FFA500" "#884190" "#58C764" "#8B795E" "#DD981A" "#8667CD"
# [21] "#28E2E5" "#447865" "#CD0BBC" "#A9885F" "#6E85D3" "#E1A10C" "#85AC57" "#2E8B57" "#4B8E4A" "#BB7663"
# [31] "#AC8543" "#2B9FD0" "#22A1E5" "#9B7F50" "#24B7E5" "#B61DAD" "#979A5B" "#CD9228" "#4FBF7A" "#A59925"
# [41] "#25C1E5" "#27D7E5" "#C39D18" "#23ACE5" "#725482" "#879631" "#9D48C7" "#34A7BA" "#26CCE5" "#46B78F"

#Get raneff predictions for kappa
cf_vm = coef(full_int_slope)$ID[,,'kappa_Intercept']
#We can collect linear scaled predictions for zmu
#these could potentially be wrong if the estimate is close to +-pi
rn = rownames(prms_vm)
rn_zmu = rn[grep(pattern = 'zmu',
                 x = rn)]
colnames(cf_vm) = colnames(prms_vm)[1:4]
#these are likely wrong for zmu
ran_prms_vm = rbind(cf_vm, prms_vm[rn_zmu,1:4])

#all draws
full_int_slope_zmu_draws = brms::as_draws_df(full_int_slope,
                                           variable = 'zmu_id') 
full_int_slope_zmu_condition_draws = brms::as_draws_df(full_int_slope,
                                           variable = 'zmu_id_condition') 
deg_pred = circular(x = deg(full_int_slope_zmu_draws +
                              est_vm$fmu_Intercept),
                    type = 'angles',
                    unit = angle_unit,
                    template = 'geographics',
                    modulo = '2pi',
                    zero = pi/2,
                    rotation = angle_rot)
deg_pred_condition = circular(x = deg(full_int_slope_zmu_draws[,1:n_indiv]+
                                      full_int_slope_zmu_condition_draws[,1:n_indiv]  +
                                        with(est_vm, 
                              fmu_Intercept + fmu_condition) ),
                    type = 'angles',
                    unit = angle_unit,
                    template = 'geographics',
                    modulo = '2pi',
                    zero = pi/2,
                    rotation = angle_rot)
deg_pred = deg_pred[,1:n_indiv]#exclude non-estimated parameters
deg_pred_condition = deg_pred_condition[,1:n_indiv]#exclude non-estimated parameters

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
title(main = 'fixed effects',
      line = -3,
      cex.main = 0.5
)
# . Model estimates -------------------------------------------------------

#all draws
full_int_slope_fixed_draws = brms::as_draws_df(full_int_slope,
                                             regex ='^b_') 
par(mar = c(0,3,4,0),
    bty = 'n')


par(mfrow = c(1,2))
with(full_int_slope_fixed_draws,
     {
      VertHist(data = inv_softplus(b_kappa_Intercept),
               ylab = 'kappa intercept',
               ylim = c(1,4))
      abline(h = median(inv_softplus(b_kappa_Intercept)),
             col = adjustcolor(col = 'cyan4',alpha.f = 0.5),
             lwd = 5)
      segments(x0 = 0,x1 = 0,
               y0 = quantile(x = inv_softplus(b_kappa_Intercept),
                              probs = 0.05),
               y1 = quantile(x = inv_softplus(b_kappa_Intercept),
                              probs = 1-0.05),
               col = 'cyan4',
               lwd = 2
               )
     }
)
abline(h = kappa_pop,
       col = 2,
       lwd = 2)

with(full_int_slope_fixed_draws,
     {
      VertHist(data = inv_softplus(b_kappa_Intercept+b_kappa_condition),
               ylab = 'kappa intercept + condition',
               ylim = c(1,4))
      abline(h = median(inv_softplus(b_kappa_Intercept+b_kappa_condition)),
             col = adjustcolor(col = 'cyan4',alpha.f = 0.5),
             lwd = 5)
      segments(x0 = 0,x1 = 0,
               y0 = quantile(x = inv_softplus(b_kappa_Intercept+b_kappa_condition),
                              probs = 0.05),
               y1 = quantile(x = inv_softplus(b_kappa_Intercept+b_kappa_condition),
                              probs = 1-0.05),
               col = 'cyan4',
               lwd = 2
               )
     }
)
abline(h = exp(log(kappa_pop)+lk_offset[2]),
       col = 2,
       lwd = 2)

par(mfrow = c(1,1))
with(full_int_slope_fixed_draws,
     {
       VertHist(data = (b_kappa_condition),
                ylab = 'condition effect on kappa',
                ylim = NULL)
       abline(h = median(b_kappa_condition),
              col = adjustcolor(col = 'cyan4',alpha.f = 0.5),
              lwd = 5)
       segments(x0 = 0,x1 = 0,
                y0 = quantile(x = b_kappa_condition,
                              probs = 0.05),
                y1 = quantile(x = b_kappa_condition,
                              probs = 1-0.05),
                col = 'cyan4',
                lwd = 2
       )
     }
)
abline(h = 0,
       col = 2,
       lwd = 2)
# with(full_int_slope_fixed_draws,
#      {
#        stripchart(x = round(inv_softplus(b_kappa_Intercept), digits = 2),
#                   pch = 20,
#                   col = gray(level = 0.5, alpha = 0.25),
#                   method = 'stack',
#                   vertical = TRUE,
#                   ylim = c(0, max(inv_softplus(b_kappa_Intercept))),
#                   xlim = c(1,2)+c(-1,1)*0.25,
#                   offset = 0.01,
#                   mgp = c(0,1,0.5),
#                   cex.main = 0.5, 
#                   bty = 'n')
#        title(main = 'Population kappa',
#              line = 0.1,
#              cex.main = 0.5)
#        abline(h = median(inv_softplus(b_kappa_Intercept)),
#               col = adjustcolor(col = 'cyan4',alpha.f = 0.5),
#               lwd = 5)
#      }
# )
# abline(h = kappa_pop,
#        col = 2,
#        lwd = 2)

with(full_int_slope_fixed_draws,
     {
       stripchart(x = round(
         inv_softplus(b_kappa_Intercept + b_kappa_condition) - inv_softplus(b_kappa_Intercept), 
         digits = 2),
         pch = 20,
         col = gray(level = 0.5, alpha = 0.25),
         method = 'stack',
         vertical = TRUE,
         ylim = c(-1, 1)*max(abs(inv_softplus(b_kappa_Intercept + b_kappa_condition) - inv_softplus(b_kappa_Intercept))),
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

with(full_int_slope_fixed_draws,
     {
       stripchart(x = round(
         sd_ID__kappa_Intercept, 
         digits = 2),
         pch = 20,
         col = gray(level = 0.5, alpha = 0.25),
         method = 'stack',
         vertical = TRUE,
         ylim = c(0, max(sd_ID__kappa_Intercept)),
         xlim = c(1,2)+c(-1,1)*0.25,
         offset = 0.01,
         mgp = c(0,1,0.5),
         cex.main = 0.5, 
         bty = 'n')
       title(main = 'Kappa SD',
             line = 0.1,
             cex.main = 0.5)
       abline(h = 0)
       abline(h = median( 
         sd_ID__kappa_Intercept
       ),
       col = adjustcolor(col = 'cyan4',alpha.f = 0.5),
       lwd = 5)
     }
)
abline(h = inv_softplus(exp(logkappa_var)),
       col = 2,
       lwd = 2)
with(full_int_slope_fixed_draws,
     {
       stripchart(x = round(
         kappa_id, 
         digits = 2),
         pch = 20,
         col = gray(level = 0.5, alpha = 0.25),
         method = 'stack',
         vertical = TRUE,
         ylim = c(0, max(kappa_id)),
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

