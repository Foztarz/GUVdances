#FOR A 'CLEAN' RUN, PRESS ctrl+shift+F10 to RESTART Rstudio
graphics.off()
# Details ---------------------------------------------------------------
#       AUTHOR:	James Foster              DATE: 2024 07 03
#     MODIFIED:	James Foster              DATE: 2024 07 03
#     MODIFIED:	James Foster              DATE: 2024 07 04
#
#  DESCRIPTION: Load dance angles, fit maximum-likelihood von Mises.
#               
#       INPUTS: 
#               
#      OUTPUTS: Plots and test statistics
#
#	   CHANGES: - Comparisons w/ emmeans
#             - Adapt delta selection for max random effects model
#             - Organise data by colour and brightness
#             - Set and inspect priors
#	   CHANGES: - attempt to speed up with future_mapply and skipping sums
#
#   REFERENCES: Sayin, S., Graving, J., et al. under review
#
#    EXAMPLES:  
#
# 
#TODO   ---------------------------------------------
#TODO   
#- Simulate null model
#- Simulate one effect
#- Test w/ optim  +
#- Test w/ quap
#- Simulate null model 
#- Simulate one effect +
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
# paired_data = TRUE # Are the data in the two columns paired (each from the same animal or group)?
set.seed(20240703)#day the script was started
n_iter = 150 # optimisation iterations

paired_data = TRUE # Are the data in the two columns paired (each from the same animal or group)?
csv_sep = ','#Is the csv comma separated or semicolon separated? For tab sep, use "\t"
angle_name = "angles" #The title of the column with angles; NO SPACES PLEASE
angle_unit = "degrees" # "degrees" or "radians"
angle_rot = "counter" # "clock" or "counter"
paired_data = TRUE

#Check the operating system and assign a logical flag (T or F)
sys_win <- Sys.info()[['sysname']] == 'Windows'
#User profile instead of home directory
if(sys_win){
  #get rid of all the backslashes
  ltp <- gsub('\\\\', '/', Sys.getenv('USERPROFILE'))#Why does windows have to make this so difficult
}else{#Root directory should be the "HOME" directory on a Mac (or Linux?)
  ltp <- Sys.getenv('HOME')#Life was easier on Mac
}

# Simulate data (not used) ------------------------------------------------
n_angles = 44
# minimum discriminable angle appears to be approx 35°
kappa_both = A1inv(0.7) #concentration around each trial mean
logkappa_var = 1.0 #scale of random variation in concentration (log units)
mu_offset = deg(rcircularuniform(n = 1))
if(paired_data)
{
kappa_indiv = A1inv(0.98) #concentration across individuals (pairs)
#mean angle in trail 1 for each individual (pair)
mu1_sim = rvonmises(n = n_angles,
                      mu = rcircularuniform(1),#random angle
                      kappa = kappa_indiv#the wider the distribution of individual biases, the greater the influence of pairing
                      )
#simulate the full dataset
sim = data.frame(
                 angle_1 = round(c(suppressWarnings( #rvonmises converts to circular and warns
                   sapply(X = mu1_sim,
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
                           condition = sort(rep(1:2, times = dt_dim[1]))
                           )
                )


# Fit two means -----------------------------------------------------------

two_m = aggregate(x = rad(angle) ~ condition,
                  data = longdata,
                  FUN = mle.vonmises,
                  bias = TRUE)

with(two_m$`rad(angle)`[1, ],
     {
arrows.circular(x = circular(deg(mu),  
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
with(two_m$`rad(angle)`[2, ],
     {
arrows.circular(x = circular(deg(mu),  
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


# Mixed-effects optimisation ----------------------------------------------

VM_LL = function(x, m, k, 
                 au = "degrees",
                 ar = "clock")
{
  # -sum( # add together # can only handle one point at a time
    -dvonmises(x = circular(x = x,
                           units = au,
                           rotation = ar), # probability density for each observed angle
              mu = m, # ML estimated mean
              kappa = k, # ML estimated concentration
              log = TRUE) # on a log scale (i.e. add instead of multiplying)
  # )
}

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
                 au = 'degrees',
                 ar = 'clock'
                )
{
  #set up population level parameter vectors
  ln = length(x)
  mm = rep(m0, times = ln)
  kk = rep(k0, times = ln)
  #adjust by conditions
  cond_2 = cond %in% unique(cond)[2]
  mm[cond_2] = mm[cond_2] + m1
  kk[cond_2] = kk[cond_2] + k1
  #adjust by ID
  mm = mm + mz
  kk = kk + kz * exp(k_sd)
  #convert to degrees if necessary
  if(au %in% 'degrees')
  {
    mm = deg(mm)
    mz = deg(mz)
  }
  #return neg log likelihood
  #for all estimates
  nll = sum( mapply(FUN = VM_LL,
               x = x,
               m = mm,
               k = exp(kk),
               au = au,
               ar = ar) )
  # #if this is lower than a circular uniform, return circular uniform -LL
  # nll = min(c(nll,
  #             -log(dcircularuniform(circular(x, units = au, rotation = ar)))
  #           ))
  #for random effects
  #on mu
  nll = nll -sum(
              dvonmises(x = circular(mz,
                                     units = au,
                                     rotation = ar),
                        mu = circular(0,
                                      units = au,
                                      rotation = ar),
                        kappa = exp(m_kappa),
                        log = TRUE
                        )
                  )
#on kappa
  nll = nll -sum(
              dnorm(x = kz, 
                        mean = 0,
                        sd = 1.0,#exp(k_sd),
                        log = TRUE
                        )
                  )
  #priors
  #priors on mu
  #mu hyperprior
  mlvm = mle.vonmises(circular(x[!cond_2], units = au, rotation  = ar))
  nll = nll - dvonmises(x = m0,
                        # mu = circular(0, units = au, rotation = ar),
                        mu = circular(mlvm$mu, units = au, rotation = ar),
                        # kappa = 0.1,
                        kappa = mlvm$kappa,
                        log = TRUE)
  nll = nll - dvonmises(x = m1,
                        mu = circular(0, units = au, rotation = ar),
                        kappa = 1,
                        log = TRUE)
  #priors on kappa
  nll = nll - dnorm(x = k0,
                    mean = mlvm$kappa,
                    sd = 0.1,
                    log = TRUE)
  nll = nll - dnorm(x = k1,
                    mean = 0,
                    sd = 0.1,
                    log = TRUE)  
  #priors on random mu
  nll = nll - dnorm(x = m_kappa,
                    mean = 1,
                    sd = 0.1,
                    log = TRUE)
  #priors on random kappa
  nll = nll - dstudent_t(x = exp(k_sd),
                    df = 3,
                    mu = 0,
                    sigma = 0.1,
                    log = TRUE)
  return(nll)
}

ME_VM_opt = function(prm, 
                     x = angle,
                     cond = cond,
                     ID = ID,
                     au = 'degrees',
                     ar = 'clock')
{
 nll = ME_VM(x = x,cond = cond,ID = ID,
        m0 = prm[1],m1 = prm[2],k0 = prm[3],k1 = prm[4],m_kappa = prm[5],k_sd = prm[6],
        mz = prm[6+1:length(unique(ID))],
        kz = prm[6+length(unique(ID)) + 1:length(unique(ID))],
        au = au,
        ar = ar
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
set.seed(20240703)
norm_est = rnorm(length(first_est))
names(norm_est) = names(first_est)

system.time(
  {
    oo = with(longdata,
              # optim(par = first_est,
              optim(par = norm_est,
               fn = ME_VM_opt,
               x = angle,
               cond = condition,
               ID = ID,
               au = angle_unit,
               ar = angle_rot,
               method = 'BFGS',
               control = list(trace = 6,
                              REPORT = 1,
                              maxit = 100))
    )
  }
)
dpar = with(oo, data.frame(t(par)))
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
