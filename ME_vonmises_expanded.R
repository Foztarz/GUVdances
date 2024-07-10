#FOR A 'CLEAN' RUN, PRESS ctrl+shift+F10 to RESTART Rstudio
graphics.off()
# Details ---------------------------------------------------------------
#       AUTHOR:	James Foster              DATE: 2024 07 09
#     MODIFIED:	James Foster              DATE: 2024 07 10
#
#  DESCRIPTION: Load dance angles, fit maximum-likelihood von Mises.
#               Adapted from ME_vonmises_test.R
#               
#       INPUTS: 
#               
#      OUTPUTS: Plots and test statistics
#
#	   CHANGES: - 
#             - 
#
#   REFERENCES: Sayin, S., Graving, J., et al. under review
#
#    EXAMPLES:  
#
# 
#TODO   ---------------------------------------------
#TODO   
#- Simulate one effect 
#- Optim approach to fixef 
#- Neaten up
#- Generalised simulation of multiple indivs w/ multiple obs
#- Expand simulated dataset
#- Test w/ quap
#- Simulate null model 
#- Simulate two way
#- Simulate interaction

# . Load packages ----------------------------------------------------------
#needs installing before first use (in Rstudio, see automatic message)
suppressMessages(#these are disturbing users unnecessarily
  {
    require(circular)#package for handling circular data
    require(CircStats)#package for circular hypothesis tests
    # require(cmdstanr)#package for Bayesian modelling via Stan
    require(brms)#package for preparing Stan models
    # require(extraDistr)#package for unusual distributions
  }
)


# Input Variables ----------------------------------------------------------

#  .  User input -----------------------------------------------------------
set.seed(20240709)#day the script was started
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
n_indiv = 10 # Number of random-effects groups
n_trials = 10 # Number of observations per random-effects group
n_conditions = 2 # Number of different conditions
# minimum discriminable angle appears to be approx 35°
kappa_pop = A1inv(0.7) #concentration around each trial mean
logkappa_var = 0.5 #scale of random variation in concentration (log units)
mu_0 = rcircularuniform(n = 1,
                        control.circular = list(units = angle_unit,
                                                rotation = angle_rot)
                        )#population intercept
mu_offset = c(0,
              rcircularuniform(n = n_conditions - 1, 
                             control.circular = list(units = angle_unit,
                                                     rotation = angle_rot)
                             )
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
conds = sort( rep(1:n_conditions, times = n_indiv*n_trials) )
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
mu_all =  mu0_sim[indivs] + mu_offset[conds]
lk_all =  lk0_sim[indivs] + lk_offset[conds]
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

par(mar =rep(0,4))
with(subset(sim, condition %in% 1),
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
with(subset(sim, condition %in% 2),
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


## Optim version ---------------------------------------------------------
#calculate negative log likelihood
VM_LL = function(x, m, k, 
                 au = "degrees",
                 ar = "clock")
{
  if(k<100)
    {
  -dvonmises(x = circular(x = x,
                          units = au,
                          rotation = ar), # probability density for each observed angle
             mu = circular(m, # conversions to circular is very important here!
                           units = au,
                           rotation = ar), # ML estimated mean
             kappa = k, # ML estimated concentration
             log = TRUE) # on a log scale (i.e. add instead of multiplying)
  }else
  {
    -dnorm(x = x,
           mean = m,
           sd = if(au %in% 'degrees')
             {deg(1/sqrt(k))}else
             {1/sqrt(k)}
           )
  }
}
  
  # Mixed-effects optimisation ----------------------------------------------
  
#TODO try m_kappa as prior (again?)
  
  ME_VM = function(x, # angle
                   cond, #condition !must be a number from 1:n_conditions!
                   ID, #individual ID !must be an index from 1:n_individuals!
                   m0, #intercept mean
                   m1, #offset for condition means
                   k0, #intercept kappa
                   k1, #offset for condition log-kappas 
                   # m_kappa, # random effects kappa for mu
                   k_sd,# random effects sd for log kappa
                   mz = NULL,# random effects for mu
                   kz = NULL,# random effects for log kappa
                   au = 'degrees', #angle unit
                   ar = 'clock', #angle rotation direction
                   plotit = FALSE
  )
  {
    #set up population level parameter vectors
    ln = length(x) # data length
    nc = length(unique(cond[!is.na(cond)])) #number of unique conditions
    ni = length(unique(ID[!is.na(ID)])) #number of unique individuals
    mm = rep(m0, times = ln) # population mean for all data
    kk = rep(k0, times = ln) # population kappa for all data
    #adjust by conditions
    if(nc >1)
    {
      mm = mm + (cond-1)*m1 # add a condition offset mean
      kk = kk + (cond-1)*k1 # add a condition offset kappa
    }
    #adjust by ID
    if(ni >1)
    {
    mm = mm + mz[ID]         # add indiv offset mean
    kk = kk + kz[ID] * k_sd  # add indiv offset kappa
    }
    #convert to degrees if necessary
    if(au %in% 'degrees')
    {
      mm = deg(mm)
      mz = deg(mz)
      m0 = deg(m0)
      m1 = deg(m1)
    }
    #return neg log likelihood
    #for all estimates
    # nll = sum( future.apply::future_mapply(FUN = VM_LL, #faster without parallel
    nll = sum( mapply(FUN = VM_LL,
                                           x = x,
                                           m = mm,
                                           k = exp(kk), # convert to kappa
                                           au = au,
                                           ar = ar) )
    #plot for troubleshooting
    if(plotit)
    {
      arrows.circular(x = circular(m0,  
                                   type = 'angles',
                                   unit = au,
                                   modulo = '2pi',
                                   zero = pi/2,
                                   template = 'geographics',
                                   rotation = ar),
                      y = A1(exp(k0)),
                      lwd = 1,
                      length = 0,
                      col = adjustcolor('cyan2', alpha.f = 0.1)
      )
      
      
      arrows.circular(x = circular(m0+m1,  
                                   type = 'angles',
                                   unit = au,
                                   modulo = '2pi',
                                   zero = pi/2,
                                   template = 'geographics',
                                   rotation = ar),
                      y = A1(exp(k0+k1)),
                      lwd = 1,
                      length = 0,
                      col = adjustcolor('blue2', alpha.f = 0.1)
      )
    }
    
    
    #for random effects
    #fit a von Mises to all random effects means
    mlvm = mle.vonmises(circular(m0+mz, units = au, rotation  = ar))
    #on mu
    #priors
    #priors on mu
    #mu hyperprior
    #try to keep m0 at the centre of mz  
    nll = nll - dvonmises(x = circular(m0, units = au, rotation = ar),
                          # mu = circular(0, units = au, rotation = ar),
                          mu = circular(mlvm$mu, units = au, rotation = ar),
                          # kappa = 0.1,
                          # kappa = mlvm$kappa,
                          kappa = A1inv(0.9), #strong prior
                          log = TRUE)
    #priors on kappa
    nll = nll - dnorm(x = k0, # this is a log kappa!
                      mean = 0,
                      sd = 0.5,
                      log = TRUE)
    nll = nll - dnorm(x = k1,
                      mean = 0,
                      sd = 0.5,
                      log = TRUE)  
    #priors on random mu
    nll = nll - dstudent_t(
                           # x = 1 / sqrt( mlvm$kappa ) , #brms rough convert to circular SD
                           x = sqrt( -2 * log( A1( mlvm$kappa ) ) ), #Mardia convert to circular SD
                           df = 1,
                           mu = 0,
                           # sigma = pi/3, #expect narrow distribution, avoid overfitting
                           sigma = pi/6, #expect narrow distribution, avoid overfitting
                           log = TRUE)
    #priors on random kappa
    nll = nll - dstudent_t(x = exp(k_sd),
                           df = 3,
                           mu = 0,
                           sigma = 1.0,
                           log = TRUE)
    return(nll)
  }
  
ME_VM_opt = function(prm, 
                     x = angle,
                     cond = cond,
                     ID = ID,
                     au = 'degrees',
                     ar = 'clock',
                     plotit = FALSE,
                     warmup = FALSE)
{
    if(warmup){prm = prm + rnorm(length(prm), sd = 0.1)}
    nll = ME_VM(x = x,cond = cond,ID = ID,
                m0 = prm[1],m1 = prm[2],k0 = prm[3],k1 = prm[4],#m_kappa = prm[5],
                k_sd = prm[5],
                mz = prm[5+1:length(unique(ID))],
                kz = prm[5+length(unique(ID)) + 1:length(unique(ID))],
                au = au,
                ar = ar,
                plotit = plotit
    )
    return(if(is.finite(nll)){nll}else{1e9})
  }
  
first_est = c( 
  m0 = 0, #intercept mean
  m1 = 1, #condition 2 mean
  k0 = 0, #intercept kappa
  k1 = 0, #condition 2 kappa 
  # m_kappa = 0, # random effects kappa for mu
  k_sd = 0,# random effects sd for log kappa
  mz_ = rep(0, times = n_indiv ),# random effects for mu
  kz_ = rep(0, times = n_indiv )# random effects for log kappa
)

#reset seed to distinguish differences in data from differences in starting params
set.seed(20240703)#day the script was started
norm_est = rnorm(length(first_est))
norm_adj = rnorm(length(first_est),sd = 0.1) #slight adjustment after warmup
names(norm_est) = names(first_est)
names(norm_adj) = names(first_est)






system.time(
  {
    oo_warmup = with(sim,
                     # optim(par = first_est,
                     optim(par = norm_est,
                           fn = ME_VM_opt,
                           x = angle,
                           cond = condition,
                           ID = indiv,
                           au = angle_unit,
                           ar = angle_rot,
                           plotit = FALSE,#don't plot warmup
                           warmup = TRUE,#parameters are not precise
                           method = 'SANN',
                           control = list(trace = 6,
                                          REPORT = 10,
                                          maxit = ceiling(n_iter/2),
                                          temp = 1)
                     )
    )  
    oo_sample = with(sim,
                     # optim(par = first_est,
                     optim(par = oo_warmup$par+norm_adj,
                           fn = ME_VM_opt,
                           x = angle,
                           cond = condition,
                           ID = indiv,
                           au = angle_unit,
                           ar = angle_rot,
                           plotit = FALSE,
                           method = 'SANN',
                           control = list(trace = 6,
                                          REPORT = 10,
                                          maxit = ceiling(n_iter/2),
                                          temp = 0.5)
                     )
    )
  }
)

dpar = with(oo_sample, data.frame(t(par)))

# Plot predictions --------------------------------------------------------


par(mar =rep(0,4))
with(subset(sim, condition %in% 1),
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
with(subset(sim, condition %in% 2),
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

par(new = T)
with(data.frame(t(oo_sample$par)),
     {
       plot.circular(x = circular(x = deg(
         m0 + sapply(X = paste0('mz_', 1:n_indiv), 
                     FUN = get, 
                     envir = as.environment(dpar))
       ), 
       type = 'angles',
       unit = angle_unit,
       template = 'geographics',
       modulo = '2pi',
       zero = pi/2,
       rotation = angle_rot
       ),
       pch = 3,
       stack = TRUE,
       bins = 360/5,
       sep = 0.5/dt_dim[1],
       col = 'cyan4',
       shrink = 1.30,
       axes = F
       )
     }
)



par(new = T)
with(dpar,
     {
       plot.circular(x = circular(x = deg(
         m0 + m1 + 
           sapply(X = paste0('mz_', 1:n_indiv), 
                  FUN = get, 
                  envir = as.environment(dpar)) 
       ), 
       type = 'angles',
       unit = 'degrees',
       modulo = '2pi',
       zero = pi/2,
       template = 'geographics',
       rotation = angle_rot
       ),
       pch = 3,
       stack = TRUE,
       bins = 360/5,
       sep = -0.5/dt_dim[1],
       col = 'darkblue',
       shrink = 1.35,
       axes = F
       )
     })
with(dpar,
     {
       arrows.circular(x = circular(deg(m0),  
                                    type = 'angles',
                                    unit = 'degrees',
                                    modulo = '2pi',
                                    zero = pi/2,
                                    template = 'geographics',
                                    rotation = angle_rot),
                       y = A1(exp(k0)),
                       lwd = 5,
                       length = 0,
                       col = adjustcolor('cyan2', alpha.f = 0.5)
       )
     }
)
with(dpar,
     {
       arrows.circular(x = circular(deg(m0+m1),  
                                    type = 'angles',
                                    unit = 'degrees',
                                    modulo = '2pi',
                                    zero = pi/2,
                                    template = 'geographics',
                                    rotation = angle_rot),
                       y = A1(exp(k0+k1)),
                       lwd = 5,
                       length = 0,
                       col = adjustcolor('blue2', alpha.f = 0.5)
       )
     }
)



arrows.circular(x = circular(deg(mu_0),  
                             type = 'angles',
                             unit = 'degrees',
                             modulo = '2pi',
                             template = 'geographics',
                             zero = 0,
                             rotation = angle_rot),
                y = A1(kappa_pop),
                lwd = 1.5,
                length = 0.1,
                col = adjustcolor('cyan4', alpha.f = 1.0)
)

arrows.circular(x = circular(deg(mu_0+mu_offset),  
                             type = 'angles',
                             unit = 'degrees',
                             modulo = '2pi',
                             zero = pi/2,
                             template = 'geographics',
                             rotation = angle_rot),
                y = A1(exp(log(kappa_pop)+lk_offset)),
                lwd = 1.5,
                length = 0.1,
                col = adjustcolor('blue2', alpha.f = 1.0)
)


# Discarded ---------------------------------------------------------------
# benchmarking 
# #faster w/ future or not?
# summary(t(replicate(60,
#                     {
#                       system.time(
#                         {
#                           oo_warmup = with(sim,
#                                            # optim(par = first_est,
#                                            optim(par = norm_est,
#                                                  fn = ME_VM_opt,
#                                                  x = angle,
#                                                  cond = condition,
#                                                  ID = indiv,
#                                                  au = angle_unit,
#                                                  ar = angle_rot,
#                                                  plotit = FALSE,#don't plot warmup
#                                                  warmup = TRUE,#parameters are not precise
#                                                  method = 'SANN',
#                                                  control = list(trace = 0,
#                                                                 REPORT = 10,
#                                                                 maxit = 10,
#                                                                 temp = 1)
#                                            )
#                           )  
#                         }
#                       )
#                     })))
