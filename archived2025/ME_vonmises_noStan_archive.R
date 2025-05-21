#FOR A 'CLEAN' RUN, PRESS ctrl+shift+F10 to RESTART Rstudio
graphics.off()
# Details ---------------------------------------------------------------
#       AUTHOR:	James Foster              DATE: 2024 07 03
#     MODIFIED:	James Foster              DATE: 2024 07 10
#
#  DESCRIPTION: Archived attempt at mixed effects von Mises without NUTS.
#               Load dance angles, fit maximum-likelihood von Mises.
#               
#       INPUTS: 
#               
#      OUTPUTS: Plots and test statistics
#
#	   CHANGES: - attempt to speed up with future_mapply and skipping sums
#             - fixed formatting error
#
#   REFERENCES: Sayin, S., Graving, J., et al. under review
#
#    EXAMPLES:  
#
# 
#TODO   ---------------------------------------------
#TODO   
#- Test w/ optim  +
#- Simulate small dataset +
#- Calculate rather than estimate m_kappa +
#- Optim approach to fixef +
#- Simulate one effect +
#- Neaten up +
#- Test w/ hmc  +
#- Test w/ quap
#- Simulate null model 

# options(warn = 2) #to catch bad circular conversions
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

Mod360.180 = function(x)
{#use atan2 to convert any angle to the range (-180,180)
  deg(
    atan2(y = sin(rad(x)),
          x = cos(rad(x))
    )
  )
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
kappa_both = A1inv(0.7) #concentration around each trial mean
logkappa_var = 1.0 #scale of random variation in concentration (log units)
mu_1 = rad(175);rcircularuniform(n = 1)#
mu_offset = rcircularuniform(n = 1)
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







# Fit two means -----------------------------------------------------------

two_m = aggregate(x = circular(angle,
                                   units = angle_unit,
                                   rotation = angle_rot) ~
                            condition,
                  data = longdata,
                  FUN = mle.vonmises,
                  bias = TRUE)
two_m = within(two_m,
               {
                 ml_est = `circular(angle, units = angle_unit, rotation = angle_rot)`
               }
               )

with(two_m$ml_est[1, ],
     {
arrows.circular(x = circular(mu,
                             type = 'angles',
                             unit = 'degrees',
                             template = 'geographics',
                             modulo = '2pi',
                             zero = pi/2,
                             rotation = angle_rot),
                y = A1(kappa),
                lwd = 3,
                col = 'cyan4')
     }
)
with(two_m$ml_est[2, ],
     {
arrows.circular(x = circular(mu,
                             type = 'angles',
                             unit = 'degrees',
                             template = 'geographics',
                             modulo = '2pi',
                             zero = pi/2,
                             rotation = angle_rot
                             ),
                y = A1(kappa),
                lwd = 3,
                col = 'darkblue')
     }
)


## Optim version ---------------------------------------------------------

#TODO#is this working?! #now it is, unit error!
VM_LL = function(x, m, k,
                 au = "degrees",
                 ar = "clock")
{
  -dvonmises(x = circular(x = x,
                          units = au,
                          rotation = ar), # probability density for each observed angle
             mu = circular(m, # conversions to circular is very important here!
                           units = au,
                           rotation = ar), # ML estimated mean
             kappa = k, # ML estimated concentration
             log = TRUE) # on a log scale (i.e. add instead of multiplying)
}

### Just one condition ---------------------------------------------------

#visualise LPD
fqt_lpd = sapply( 0:359,
             function(m){
               sum( future.apply::future_mapply(
                 FUN = VM_LL,
                  x = subset(longdata, condition %in% 1)$angle,
                  m = m,
                  k = A1inv(0.75),
                  au = angle_unit,
                  ar = angle_rot) ) })

#uniform reference
unif_lpd =sum( future.apply::future_mapply(
  FUN = VM_LL,
  x = subset(longdata, condition %in% 1)$angle,
  m = 0,
  k = 0,
  au = angle_unit,
  ar = angle_rot) )

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

lines.circular(x = circular(0:359,
                            units = angle_unit,
                            rotation = angle_rot,
                            template = 'geographics',
                            modulo = '2pi',
                            zero = pi/2
                            ),
               y = (1-fqt_lpd/max(fqt_lpd))-1,
               col = 'purple',
               lwd = 5)



One_VM = function(x, # angle
                  cond, #condition
                  m0, #intercept mean
                  k0, #intercept kappa
                  au = 'degrees', #angle unit
                  ar = 'clock', #angle rotation direction
                  plotit = FALSE
)
{
  #set up population level parameter vectors
  ln = length(x) # data length
  mm = rep(m0, times = ln) # population mean for all data
  kk = rep(k0, times = ln) # population kappa for all data
  #convert to degrees if necessary
  if(au %in% 'degrees')
  {
    mm = deg(mm)
  }
  #return neg log likelihood
  #for all estimates
  nll = sum( future.apply::future_mapply(FUN = VM_LL,
                                         x = x,
                                         m = mm,
                                         k = exp(kk),
                                         au = au,
                                         ar = ar) )
  #priors
  #priors on mu
  #mu hyperprior
  # mlvm = mle.vonmises(circular(x, units = au, rotation  = ar))
  # nll = nll - dvonmises(x = m0,
  #                       # mu = circular(0, units = au, rotation = ar),
  #                       mu = circular(mlvm$mu, units = au, rotation = ar),
  #                       # kappa = 0.1,
  #                       kappa = mlvm$kappa,
  #                       log = TRUE)
  #priors on kappa
  nll = nll - dnorm(x = k0, # this is a log kappa!
                    # mean = log(mlvm$kappa),
                    mean = 0,
                    sd = 0.25,
                    log = TRUE)


  if(plotit)
  {
  arrows.circular(x = circular(deg(m0),
                               type = 'angles',
                               unit = au,
                               modulo = '2pi',
                               zero = pi/2,
                               template = 'geographics',
                               rotation = ar),
                  y = A1(exp(k0)),
                  lwd = 0.5,
                  length = 0,
                  col = adjustcolor('cyan2', alpha.f = 0.1)
  )
  }

  return(nll)
}

One_VM_opt = function(prm,
                     x = angle,
                     cond = cond,
                     ID = ID,
                     au = 'degrees',
                     ar = 'clock',
                     plotit = FALSE)
{
  nll = One_VM(x = x,cond = cond,
              m0 = prm[1],k0 = prm[2],
              au = au,
              ar = ar,
              plotit
  )
  return(if(is.finite(nll)){nll}else{1e9})
}


#reset seed to distinguish differences in data from differences in starting params
set.seed(20240703)#day the script was started
norm_one = rnorm(2)
names(norm_one) = c('m0', 'k0')


#plotting version

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



system.time(
  {
    oo_one = with(subset(longdata, condition %in% 1),
                 optim(par = norm_one,
                       fn = One_VM_opt,
                       x = angle,
                       cond = condition,
                       au = angle_unit,
                       ar = angle_rot,
                       plotit = TRUE,
                       method = 'SANN',
                       control = list(trace = 1,
                                      REPORT = 10,
                                      maxit = n_iter/4,
                                      abstol = 0,
                                      temp = 0.1))
    )
  }
)

dpar_one = with(oo_one, data.frame(t(par)))

with(dpar_one,
     {
       print(
         list(
           deg(m0) %% 360,
           A1(exp(k0))
         ),
         digits = 3 )
     }
)



### Two conditions -------------------------------------------------------




ML_VM = function(x, # angle
                 cond, #condition
                 m0, #intercept mean
                 m1, #condition 2 mean
                 k0, #intercept kappa
                 k1, #condition 2 kappa
                 au = 'degrees', #angle unit
                 ar = 'clock', #angle rotation direction
                 plotit = FALSE
)
{
  #set up population level parameter vectors
  ln = length(x) # data length
  mm = rep(m0, times = ln) # population mean for all data
  kk = rep(k0, times = ln) # population kappa for all data
  #adjust by conditions
  cond_2 = cond %in% unique(cond)[2] # just two conditions, find index of 2nd
  mm[cond_2] = mm[cond_2] + m1 # add a condition offset mean
  kk[cond_2] = kk[cond_2] + k1
  #convert to degrees if necessary
  if(au %in% 'degrees')
  {
    mm = deg(mm)
  }
  #return neg log likelihood
  #for all estimates
  nll = sum( future.apply::future_mapply(FUN = VM_LL,
                                         x = x,
                                         m = mm,
                                         k = exp(kk),
                                         au = au,
                                         ar = ar) )
  #priors
  #priors on mu
  #mu hyperprior
  # mlvm = mle.vonmises(circular(x[!cond_2], units = au, rotation  = ar))
    # mlvm = mle.vonmises(circular(x, units = au, rotation  = ar))
    # nll = nll - dvonmises(x = m0,
    #                       # mu = circular(0, units = au, rotation = ar),
    #                       mu = circular(mlvm$mu, units = au, rotation = ar),
    #                       # kappa = 0.1,
    #                       kappa = mlvm$kappa,
    #                       log = TRUE)
    # nll = nll - dvonmises(x = m1,
    #                       mu = circular(0, units = au, rotation = ar),
    #                       # kappa = A1inv( sqrt(-log(0.05)/(sum(cond_2))) ), #kappa is smallest sig. mean vector
    #                       # kappa = 0.25, #mid level kappa
    #                       kappa = mlvm$kappa,
    #                       log = TRUE)
  #priors on kappa
  nll = nll - dnorm(x = k0, # this is a log kappa!
                    # mean = log(mlvm$kappa),
                    mean = 0,
                    sd = 0.25,
                    log = TRUE)
  nll = nll - dnorm(x = k1,
                    mean = 0,
                    sd = 0.25,
                    log = TRUE)


  if(plotit)
  {
         arrows.circular(x = circular(deg(m0),
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


         arrows.circular(x = circular(deg(m0+m1),
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

  return(nll)
}

ML_VM_opt = function(prm,
                     x = angle,
                     cond = cond,
                     ID = ID,
                     au = 'degrees',
                     ar = 'clock',
                     plotit = FALSE)
{
  nll = ML_VM(x = x,cond = cond,
              m0 = prm[1],m1 = prm[2],k0 = prm[3],k1 = prm[4],
              au = au,
              ar = ar,
              plotit
  )
  return(if(is.finite(nll)){nll}else{1e9})
}


#reset seed to distinguish differences in data from differences in starting params
set.seed(20240703)#day the script was started
norm_ml = rnorm(4)
names(norm_ml) = c('m0', 'm1', 'k0', 'k1')


#plotting version

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


system.time(
  {
    oo_ml = with(longdata,
              optim(par = norm_ml,
                    fn = ML_VM_opt,
                    x = angle,
                    cond = condition,
                    au = angle_unit,
                    ar = angle_rot,
                    plotit = TRUE,
                    method = 'SANN',
                    control = list(trace = 1,
                                   REPORT = 10,
                                   maxit = ceiling(n_iter/4),
                                   abstol = 0,
                                   temp = 1)) #temperature of 1 is a Metropolis-Hastings algorithm
    )
  }
)

dpar_ml = with(oo_ml, data.frame(t(par)))

with(dpar_ml,
     {
       print(
         list(
           deg(m0) %% 360,
           deg(m0+m1) %% 360,
           A1(exp(k0)),
           A1(exp(k0 + k1))
         ),
         digits = 3 )
     }
)


## Plot predictions together ---------------------------------------------

#
# with(dpar_ml,
#      {
#        arrows.circular(x = circular(deg(m0),
#                                     type = 'angles',
#                                     unit = 'degrees',
#                                     modulo = '2pi',
#                                     zero = pi/2,
#                                     template = 'geographics',
#                                     rotation = angle_rot),
#                        y = A1(exp(k0)),
#                        lwd = 5,
#                        length = 0,
#                        col = adjustcolor('cyan2', alpha.f = 0.5)
#        )
#      }
# )
# with(dpar_ml,
#      {
#        arrows.circular(x = circular(deg(m0+m1),
#                                     type = 'angles',
#                                     unit = 'degrees',
#                                     modulo = '2pi',
#                                     zero = pi/2,
#                                     template = 'geographics',
#                                     rotation = angle_rot),
#                        y = A1(exp(k0+k1)),
#                        lwd = 5,
#                        length = 0,
#                        col = adjustcolor('blue2', alpha.f = 0.5)
#        )
#      }
# )


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

with(two_m$ml_est[1, ],
     {
       arrows.circular(x = circular(mu,
                                    type = 'angles',
                                    unit = 'degrees',
                                    template = 'geographics',
                                    modulo = '2pi',
                                    zero = pi/2,
                                    rotation = angle_rot),
                       y = A1(kappa),
                       lwd = 3,
                       col = 'cyan4')
     }
)
with(two_m$ml_est[2, ],
     {
       arrows.circular(x = circular(mu,
                                    type = 'angles',
                                    unit = 'degrees',
                                    template = 'geographics',
                                    modulo = '2pi',
                                    zero = pi/2,
                                    rotation = angle_rot),
                       y = A1(kappa),
                       lwd = 3,
                       col = 'darkblue')
     }
)

with(dpar_ml,
     {
       arrows.circular(x = circular(deg(m0),
                                    type = 'angles',
                                    unit = 'degrees',
                                    template = 'geographics',
                                    modulo = '2pi',
                                    zero = pi/2,
                                    rotation = angle_rot),
                       y = A1(exp(k0)),
                       lwd = 5,
                       length = 0,
                       col = adjustcolor('cyan2', alpha.f = 0.5)
       )
     }
)
with(dpar_ml,
     {
       arrows.circular(x = circular(deg(m0+m1),
                                    type = 'angles',
                                    unit = 'degrees',
                                    template = 'geographics',
                                    modulo = '2pi',
                                    zero = pi/2,
                                    rotation = angle_rot),
                       y = A1(exp(k0+k1)),
                       lwd = 5,
                       length = 0,
                       col = adjustcolor('blue2', alpha.f = 0.5)
       )
     }
)


# Mixed-effects optimisation ----------------------------------------------


ME_VM = function(x, # angle
                 cond, #condition
                 ID, #individual ID
                 m0, #intercept mean
                 m1, #condition 2 mean
                 k0, #intercept kappa
                 k1, #condition 2 kappa
                 m_kappa, # random effects kappa for mu
                 k_sd,# random effects sd for log kappa
                 mz,# random effects for mu
                 kz,# random effects for log kappa
                 au = 'degrees', #angle unit
                 ar = 'clock', #angle rotation direction
                 plotit = FALSE
                )
{
  #set up population level parameter vectors
  ln = length(x) # data length
  mm = rep(m0, times = ln) # population mean for all data
  kk = rep(k0, times = ln) # population kappa for all data
  #adjust by conditions
  cond_2 = cond %in% unique(cond)[2] # just two conditions, find index of 2nd
  mm[cond_2] = mm[cond_2] + m1 # add a condition offset mean
  kk[cond_2] = kk[cond_2] + k1
  #adjust by ID
  mm = mm + mz
  # mm = mm + mz * sqrt(-2*log(A1(exp(m_kappa)))) #convert kappa to sd
  # kk = kk + kz * exp(k_sd) #exponentiate log_kappa #I think this was wrong
  kk = kk + kz * k_sd #exponentiate log_kappa
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
  nll = sum( future.apply::future_mapply(FUN = VM_LL,
               x = x,
               m = mm,
               k = exp(kk),
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
    nn = mle.vonmises(circular(mm, units = au, rotation  = ar))
    arrows.circular(x = circular(nn$mu,
                                 type = 'angles',
                                 unit = au,
                                 modulo = '2pi',
                                 zero = pi/2,
                                 template = 'geographics',
                                 rotation = ar),
                    y = A1(nn$kappa),
                    lwd = 1,
                    length = 0,
                    col = adjustcolor('salmon', alpha.f = 0.1)
    )
  }


  #for random effects
  #may be unnecessary
  #on mu
  #just find the raneff SD
  mle_kappa = mle.vonmises(x =
                            circular(mz,
                                    units = au,
                                    rotation = ar)
                           )$kappa
  # nll = nll - dnorm(x = m_kappa, # check that m_kappa is representative
  #                   mean = log(mle_kappa),
  #                   sd = 0.1,
  #                   log = TRUE) # should track as closely as possible
  #priors
  #priors on mu
  #mu hyperprior
    # mlvm = mle.vonmises(circular(x[!cond_2], units = au, rotation  = ar))
    # mlvm = mle.vonmises(circular(x, units = au, rotation  = ar))
  #try to keep m0 at the centre of mz
  mlvm = mle.vonmises(circular(m0+mz, units = au, rotation  = ar))
    nll = nll - dvonmises(x = circular(m0, units = au, rotation  = ar),
                          # mu = circular(0, units = au, rotation = ar),
                          mu = circular(mlvm$mu, units = au, rotation = ar),
                          # kappa = 0.1,
                          # kappa = mlvm$kappa,
                          # kappa = A1inv(0.90), # strong prior
                          kappa = A1inv(0.99), # very strong prior
                          log = TRUE)
    # nll = nll - dvonmises(x = m1,
    #                       mu = circular(0, units = au, rotation = ar),
    #                       kappa = A1inv( sqrt(-log(0.05)/(sum(cond_2))) ), #kappa is smallest sig. mean vector
    #                       log = TRUE)
  #priors on kappa
  #could be improved?
  nll = nll - dnorm(x = k0, # this is a log kappa!
                    mean = 0,
                    # mean = log(mlvm$kappa),
                    sd = 0.3,
                    # sd = 1.0,
                    log = TRUE)
  nll = nll - dnorm(x = k1,
                    mean = 0,
                    sd = 0.3,
                    # sd = 1.0,
                    log = TRUE)
  #priors on random mu
  # nll = nll - dnorm(x = m_kappa,
  #                   mean = 1,
  #                   sd = 0.5,
  #                   log = TRUE)
  # nll = nll - dstudent_t(x = sqrt( -2 * log( A1( exp(m_kappa) ) ) ), #convert to circular SD
  nll = nll - dstudent_t(x = sqrt( -2 * log( A1( mle_kappa ) ) ), #convert to circular SD
                    df = 1,
                    mu = 0,
                    sigma = pi/3, #expect narrow distribution, avoid overfitting
                    # sigma = pi/6, #expect narrow distribution, avoid overfitting
                    log = TRUE)
  #priors on random kappa
  nll = nll - dstudent_t(x = exp(k_sd),
                    df = 3,
                    mu = 0,
                    sigma = 0.3,
                    # sigma = 1.0,
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
        m0 = prm[1],m1 = prm[2],k0 = prm[3],k1 = prm[4],m_kappa = prm[5],k_sd = prm[6],
        mz = prm[6+1:length(unique(ID))],
        kz = prm[6+length(unique(ID)) + 1:length(unique(ID))],
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
              m_kappa = 0, # random effects kappa for mu
              k_sd = 0,# random effects sd for log kappa
              mz_ = rep(0, times = length(unique(longdata$ID)) ),# random effects for mu
              kz_ = rep(0, times = length(unique(longdata$ID)) )# random effects for log kappa
            )

#reset seed to distinguish differences in data from differences in starting params
set.seed(20240703)#day the script was started
norm_est = rnorm(length(first_est))
norm_adj = rnorm(length(first_est),sd = 0.1) #slight adjustment after warmup
names(norm_est) = names(first_est)
names(norm_adj) = names(first_est)
#
# #print uniform neg LL for reference
# with(longdata,
#      print(-sum(log(dcircularuniform(x =
#                                         circular(angle, units = angle_unit, rotation = angle_rot)
#                             ) ))
#            )
# )
# #pring MLE neg LL for reference
# print(
#   with(data.frame(two_m$ml_est),
#        -sum( c(
#          sum(
#            dvonmises(
#              circular(subset(longdata, condition ==1)$angle,
#                       units = angle_unit,
#                       rotation = angle_rot),
#              mu = circular(deg(mu[[1]]),
#                            units = angle_unit,
#                            rotation = angle_rot),
#              kappa = kappa[[1]],
#              log = TRUE)
#          )
#        ),
#
#        sum(
#          dvonmises(
#            circular(subset(longdata, condition ==2)$angle,
#                     units = angle_unit,
#                     rotation = angle_rot),
#            mu = circular(deg(mu[[2]]),
#                          units = angle_unit,
#                          rotation = angle_rot),
#            kappa = kappa[[2]],
#            log = TRUE)
#        )
#        )
#   ),
#   digits = 5
# )

# Fit model ---------------------------------------------------------------



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

# working poorly
system.time(
  {
    oo_warmup = with(longdata,
              # optim(par = first_est,
              optim(par = norm_est,
               fn = ME_VM_opt,
               x = angle,
               cond = condition,
               ID = ID,
               au = angle_unit,
               ar = angle_rot,
               plotit = FALSE,#don't plot warmup
               warmup = FALSE,#parameters are not precise
               method = 'SANN',#Nelder-Mead',
               control = list(trace = 6,
                              REPORT = 10,
                              maxit = ceiling(n_iter/2),
                              abstol = 0,
                              temp = 1)
              )
    )
    oo_sample = with(longdata,
              # optim(par = first_est,
              optim(par = oo_warmup$par,#+norm_adj,
               fn = ME_VM_opt,
               x = angle,
               cond = condition,
               ID = ID,
               au = angle_unit,
               ar = angle_rot,
               plotit = TRUE,
               method = 'SANN',#Nelder-Mead',
               control = list(trace = 6,
                              REPORT = 10,
                              maxit = ceiling(n_iter/2),
                              abstol = 0,
                              temp = 0.5#0.1
                              )
              )
    )
  }
)

# Test w/ Hamiltonian Monte Carlo -----------------------------------------
#extremely slow and difficult to find appropriate parameters

require(rhmc)
ME_VM_HMC = function(prm)
{
  with(longdata,
       {
        ME_VM_opt(prm = prm,
                  x = angle,
                  cond = condition,
                  ID = ID,
                  au = angle_unit,
                  ar = angle_rot,
                  plotit = TRUE
                  )
       }
  )
}

#warmup
Warmup_HMC = function(hmcp = c(2, 0.05, 0.0), fun, pr, iter = 1e2)
{
  -hmc(f = fun,
      init = pr,
      numit = iter,
      L = abs(ceiling(hmcp[1])),
      eps = hmcp[2],
      mass = hmcp[3])$ar
}

#TODO find good starting values to make warmup more efficient
# system.time(
#   {
#     hmc_par = optim(#par = c(16, 0.3, 0.1),
#                     par = c(2, 0.01, 0.9),
#                     fn = Warmup_HMC,
#                     fun = ME_VM_HMC,
#                     pr = oo_warmup$par,
#                     iter = 5e0,
#                     method = "L-BFGS-B",
#                     lower = c(1, 0, 0),
#                     upper = c(30, 0.99, 0.99),
#                     control = list(trace = 6,
#                                    REPORT = 10,
#                                    maxit = 10,
#                                    abstol = 1-0.9,
#                                    temp = 0.5#0.1
#                     )
#                     )
#   }
# )


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

system.time(
  {
oo_hmc = hmc(f = ME_VM_HMC,
                  init = oo_warmup$par,
                  numit = 10,
                  L = 5,
                  eps = 0.05,
                  mass = 0.9)
}
)

# Extract parameters ------------------------------------------------------


dpar = with(oo_sample, data.frame(t(par)))

with(dpar,
     {
      print(
        list(
          deg(m0) %% 360,
          deg(m0+m1) %% 360,
          A1(exp(k0)),
          A1(exp(k0 + k1)),
          A1(exp(m_kappa)),
          exp(k_sd)
          ),
        digits = 3 )
      }
     )
#for reference, these are the mle params
      with(dpar,
           {
             print(
               cbind(
               c(
                 deg(m0) %% 360,
                 deg(m0+m1) %% 360),
 deg(unlist(two_m$ml_est[,'mu'])) %% 360
 ),
 digits = 5)
           }
               )
      with(dpar,
           {
             print(
               cbind(
                 c(exp(k0),
                   exp(k0 + k1) ),
                   unlist(
two_m$ml_est[,'kappa']
)
),
digits = 5
)
           }
)

# # Check fixef LL ----------------------------------------------------------
# print(-sum( c(
#       with(dpar,
#      sum(
#        dvonmises(
#          circular(subset(longdata, condition ==1)$angle,
#                   units = angle_unit,
#                   rotation = angle_rot),
#          mu = circular(deg(m0),
#                        units = angle_unit,
#                        rotation = angle_rot),
#          kappa = exp(k0),
#          log = TRUE)
#      )
# ),
#
# with(dpar,
#      sum(
#        dvonmises(
#          circular(subset(longdata, condition ==2)$angle,
#                   units = angle_unit,
#                   rotation = angle_rot),
#          mu = circular(deg(m0+m1),
#                        units = angle_unit,
#                        rotation = angle_rot),
#          kappa = exp(k0+k1),
#          log = TRUE)
#      )
# ) )
# ),
# digits = 5
# )
#
# print(
#             with(data.frame(two_m$ml_est),
# -sum( c(
#      sum(
#        dvonmises(
#          circular(subset(longdata, condition ==1)$angle,
#                   units = angle_unit,
#                   rotation = angle_rot),
#          mu = circular(deg(mu[[1]]),
#                        units = angle_unit,
#                        rotation = angle_rot),
#          kappa = kappa[[1]],
#          log = TRUE)
#      )
# ),
#
#      sum(
#        dvonmises(
#          circular(subset(longdata, condition ==2)$angle,
#                   units = angle_unit,
#                   rotation = angle_rot),
#          mu = circular(deg(mu[[2]]),
#                        units = angle_unit,
#                        rotation = angle_rot),
#          kappa = kappa[[2]],
#          log = TRUE)
#      )
# )
# ),
# digits = 5
# )

# Plot predictions --------------------------------------------------------


par(mar =rep(0,4), pty = 's')
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
col = 'cyan4',
pty = 's'
)
par(new = T, pty = 's')
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




par(new = T, pty = 's')
with(data.frame(t(oo_sample$par)),
           {
             plot.circular(x = circular(x = deg(
               m0 + sapply(X = paste0('mz_', 1:n_angles),
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



par(new = T, pty = 's')
with(dpar,
     {
       plot.circular(x = circular(x = deg(
                                 m0 + m1 +
                                   sapply(X = paste0('mz_', 1:n_angles),
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



arrows.circular(x = circular(deg(mu_1),
                            type = 'angles',
                            unit = 'degrees',
                            modulo = '2pi',
                            template = 'geographics',
                            zero = 0,
                            rotation = angle_rot),
               y = A1(kappa_both),
               lwd = 1.5,
               length = 0.1,
               col = adjustcolor('cyan4', alpha.f = 1.0)
)

arrows.circular(x = circular(deg(mu_1+mu_offset),
                              type = 'angles',
                              unit = 'degrees',
                              modulo = '2pi',
                              zero = pi/2,
                              template = 'geographics',
                              rotation = angle_rot),
                 y = A1(kappa_both),
                 lwd = 1.5,
                 length = 0.1,
                 col = adjustcolor('blue2', alpha.f = 1.0)
 )

