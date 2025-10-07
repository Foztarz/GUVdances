#FOR A 'CLEAN' RUN, PRESS ctrl+shift+F10 to RESTART Rstudio
graphics.off()
# Details ---------------------------------------------------------------
#       AUTHOR:	James Foster              DATE: 2025 10 07
#     MODIFIED:	James Foster              DATE: 2025 10 07
#
#  DESCRIPTION: Attempt build a no predictors mixture model using BRMS and Stan.
#               Modified from GUV_CircMod_mixture.R, which is currently converging poorly
#               
#       INPUTS: 
#               
#      OUTPUTS: Plots and test statistics
#
#	   CHANGES: - 
#
#   REFERENCES: Sayin S, Couzin-Fuchs E, Petelski I, Günzel Y, Salahshour M, 
#               Lee CY, Graving JM, Li L, Deussen O, Sword GA, et al. (2025)
#               The behavioral mechanisms governing collective motion in swarming locusts.
#               Science.
#               387(6737):995–791
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
#- Select only UV-dim +
#- Simulate bimodal
#- Subset test +
#- Eliminate bimodality in fmu1 posterior
#- Find priors for mu2
#- Find priors for zmu
#- Full dataset test


# Set up workspace --------------------------------------------------------
angle_unit = 'degrees'
angle_rot = 'clock'

#Check the operating system and assign a logical flag (T or F)
sys_win = Sys.info()[['sysname']] == 'Windows'
#find file location
here_path = tryCatch(expr = #look in the folder containing this file: sys.frame(1)$ofile
                       {file.path(dirname(sys.frame(1)$ofile))},
                     error = function(e)
                     {#if that fails, try to find the "Documents" folder
                       file.path(Sys.getenv('HOME'))
                     }
)

tryCatch(expr = #try to load functions from the folder containing this file
           {
             source(file = file.path(here_path,
                                     'GUV_functions.R',
                                     fsep = if(sys_win){'\\'}else{'/'}
                                     )
                    )
           },
         error = function(e)
         {#if that fails, ask the user
          
             path_functions =  
               if(sys_win){#choose.files is only available on Windows
                   message('\n\nPlease select the "GUV_functions.R" file\n\n')
                   Sys.sleep(0.5)#goes too fast for the user to see the message on some computers
                   choose.files(
                     default = file.path(here_path, "*.R"),#For some reason this is not possible in the "root" user
                     caption = 'Please select the "GUV_functions.R" file'
                   )
               }else{
                 message('\n\nPlease select the "GUV_functions.R" file\n\n')
                 Sys.sleep(0.5)#goes too fast for the user to see the message on some computers
                 file.choose(new=F)
               }
           #show the user the path they have selected
           if(is.null(path_functions) | !length(path_functions))
           {stop('No file selected.')}else
           {print(path_functions)}
           source(path_functions)
         }
)

path_file = file.path(here_path, "1Data/colour_dance_reorg.csv")
if(file.exists(path_file))
{
  print(path_file)
}else
{
  # set path to file
  if(sys_win){#choose.files is only available on Windows
    message('\n\nPlease select the ".csv" file\n\n')
    Sys.sleep(0.5)#goes too fast for the user to see the message on some computers
    path_file = choose.files(
      default = file.path(here_path, "*.csv"),#For some reason this is not possible in the "root" user
      caption = 'Please select the "colour_dance_reorg.csv" file'
    )
  }else{
    message('\n\nPlease select the "colour_dance_reorg.csv" file\n\n')
    Sys.sleep(0.5)#goes too fast for the user to see the message on some computers
    path_file = file.choose(new=F)
  }
  #show the user the path they have selected
  if(is.null(path_file) | !length(path_file))
  {stop('No file selected.')}else
  {print(path_file)}
}
# Custom family --------------------------------------------------------
#In the "unwrap" family, all variables 
unwrap_von_mises = custom_family(
  "unwrap_von_mises", dpars = c("mu", "kappa"),
  links = c('identity',#brms cannot accept custom link functions, do via nl instead
            "softplus"), 
  lb = c(-pi, 0),
  ub = c(pi, NA),
  type = "real"
)

### Stan variables ---------------------------------------------------------


### Functions ------------------------------------------------------------

stan_unwrap_fun = stanvar(scode = "
  real unwrap_von_mises_lpdf(real y, real mu, real kappa) {
    return von_mises_lpdf(y | mod_circular(mu), kappa);
  }
  real unwrap_von_mises_rng(real mu, real kappa) {
    return von_mises_rng( mod_circular(mu) , kappa);
  }
  real unwrap_von_mises_vect_lpdf(vector y, real mu, real kappa) {
    real tmp = 0;
    for(i in 1:size(y))
    {
    tmp = tmp + unwrap_von_mises_lpdf(y[i] | mu, kappa);
    }
      return tmp;
  }
",
                          block = 'functions') + 
  stanvar(scode = "
  real mod_circular(real y) {
    return fmod(y + pi(), 2*pi()) - pi();
  }

  real inv_mod_circular(real y) {
    return mod_circular(y);
  }
",
          block = 'functions') 


### Parameters -----------------------------------------------------------


#concentration of random effects on fixed effect on mean angle
kappamu_param = stanvar(scode = "
real kappamu1;
real kappamu2;
                           ",
                            block = "parameters") + 
  stanvar(scode = "
real kappa_mu1 = log1p_exp(kappamu1);
real kappa_mu2 = log1p_exp(kappamu1+kappamu2);
          ", 
          block = 'genquant')
#concentration of random effects on fixed effect on mean angle
uni_kappamu_param = stanvar(scode = "
real kappamu1;
                           ",
                            block = "parameters") + 
  stanvar(scode = "
real kappa_mu1 = log1p_exp(kappamu1);
          ", 
          block = 'genquant')


stanvars_nopred = stan_unwrap_fun + kappamu_param
stanvars_nopred_uni = stan_unwrap_fun + uni_kappamu_param


# Input Variables ----------------------------------------------------------
all_plots = FALSE # to speed up

# Read in the data and format ---------------------------------------------

#select the reorganised data
cd = read.table(file = path_file, 
                header = T, 
                sep  = ',')
if(all_plots)
{
  View(cd)
}
cd = within(cd,
            {
              ID = as.factor(ID) # beedance identifier as a factor
              date = as.factor(date) # date as a factor
              signed_angle = Mod360.180(bearing)  # bearing between -180 and 180
              angle = circular(rad(signed_angle),# bearing between -pi and pi
                               rotation = 'clock') # circular format suppresses later warnings
            }
)

u_id = with(cd, unique(ID)) # unique beedances
length(u_id)#169 beedances

# Fit model ---------------------------------------------------------------

## Simpler variable names ------------------------------------------------
cd = within(cd,
            {
              BR = brightn
              CL = colour
              DT = date
              RN = run
            }
)

## Make a subset for testing with N individuals ----
cd_subs = subset(x = cd,
                 subset = ID %in%
                   sample(u_id,
                          size =  10, # just 10 randomly chosen indivs
                          replace = FALSE) &
                   BR %in% 'l' & # just UV-dim
                   CL %in% 'u'
)

#the model is not taking so long to fit, stick with whole dataset
cd_subs = cd

cd_subs = within(cd_subs,
                 {
                 angle = as.numeric(angle)#remove circular formatting for predictions  
                 }
                 )

## Formula ---------------------------------------------------------------


formula_uni = bf(#modulus may not be necessary, included in lpd function
  formula = angle ~ mu,
  nlf(mu ~ fmu1 + zmu1), 
  fmu1 ~ 1, # all conditions
  zmu1 ~ 0 + ID, # mean angle combines fixed and random effects
  nlf(kappa ~ k1),
  k1 ~ 1 + (1 |ID), #for kappa this occurs in linear space, and BRMS can set it up automatically
  family = unwrap_von_mises,#mod mu, kappa via the softplus
  nl = TRUE)#to accept user-defined extra parameters (zmu) we need to treat the formula as nonlinear

formula_mix = bf(#modulus may not be necessary, included in lpd function
  formula = angle ~ mu,
  nlf(mu1 ~ fmu1 + zmu1), 
  nlf(mu2 ~ fmu1 + fmu2 + zmu1 + zmu2), #assume it is the same, but with some small increment
  fmu1 ~ 1, # all conditions
  zmu1 ~ 0 + ID, # mean angle combines fixed and random effects
  fmu2 ~ 1, #just UV-dim
  zmu2 ~ 0 + ID, # mean angle combines fixed and random effects
  nlf(kappa1 ~ k1),
  nlf(kappa2 ~ k1), #for 2nd mean, this is specific to individual and condition combination
  k1 ~ 1 + (1 |ID), #for kappa this occurs in linear space, and BRMS can set it up automatically
  theta1 ~ 1 + (1 |ID), #for mixture weighting, this is specific to individual and condition combination
  family = mixture(unwrap_von_mises,#mod mu, kappa via the softplus
                   unwrap_von_mises
  ),#mixture, specified by the log ratio of theta1 : theta2
  nl = TRUE)#to accept user-defined extra parameters (zmu) we need to treat the formula as nonlinear


## Priors ----------------------------------------------------------------

#priors for mu
pr_mu_uni = 
  prior(normal(0,pi()/2), class = b,  nlpar = 'fmu1', coef = 'Intercept') + # closer to 0°
  set_prior(paste("target +=", 
                  'unwrap_von_mises_vect_lpdf(b_zmu1 | 0, log1p_exp(kappamu1))',
                  '+ normal_lpdf(b_zmu1 | 0, 2*pi())'# additional prior to keep estimates from walking around the circle
  ),
  check = FALSE) +
  set_prior("target += normal_lpdf(kappamu1 | 3.0, 3.0)", #prior to higher values, indiv differences should be small
            check = FALSE) 
pr_mu_mix = 
  pr_mu_uni + #same as unimodal, plus priors for 2nd mean
  prior(normal(-pi(),pi()/2), class = b, nlpar = 'fmu2', coef = 'Intercept') + # closer to 180°
  set_prior(paste("target +=", 
                  'unwrap_von_mises_vect_lpdf(b_zmu1 | 0, log1p_exp(kappamu1))',
                  '+ normal_lpdf(b_zmu1 | 0, 2*pi())'# additional prior to keep estimates from walking around the circle
  ),
  check = FALSE) +
  set_prior("target += normal_lpdf(kappamu1 | 3.0, 3.0)", #prior to higher values, indiv differences should be small
            check = FALSE) + 
  set_prior(paste("target +=", 
                  'unwrap_von_mises_vect_lpdf(b_zmu2 | 0, log1p_exp(kappamu1 + kappamu2))',
                  '+ normal_lpdf(b_zmu2 | 0, 2*pi())'# additional prior to keep estimates from walking around the circle
  ),
  check = FALSE) +
  set_prior("target += normal_lpdf(kappamu2 | 0.0, 1.0)", #prior to higher values, indiv differences should be small
            check = FALSE)
#priors for kappa
pr_kappa_uni = 
  prior(normal(3.0,3.0), class = b, nlpar = 'k1', coef = 'Intercept') + # bias to oriented
  prior(student_t(3, 0, 3.0), class = sd, nlpar = 'k1')  # weak bias no differences
pr_kappa_mix = pr_kappa_uni  #identical
#priors for theta (mixture weight)
pr_theta_mix = 
  prior(normal(1.0,0.25), class = Intercept, dpar = 'theta1') + # force to mu1 as primary #20251007 approx 30% of data are at mu2, so setting to 1.0
  prior(student_t(3, 0, 0.5), class = sd, dpar = 'theta1') 

#all unimodal priors
pr_uni = pr_mu_uni + pr_kappa_uni
#all mixture priors
pr_mix = pr_mu_mix + pr_kappa_mix + pr_theta_mix


## Save Stancode -----------------------------------------------------------


sc_mix = make_stancode(formula = formula_mix,
                       data = cd_subs,
                       prior = pr_mix,
                       stanvars = stanvars_nopred)

write.table(x = sc_mix,
            file = file.path(dirname(path_file),
                             'sc_nopred_mix.stan'),
            quote = FALSE,
            col.names = FALSE,
            row.names = FALSE)

## Run model -------------------------------------------------------------
wup = 500
sam = 200


### Unimodal -------------------------------------------------------------
#very long compile time
system.time(
  {
    
    np_uni = brm(formula = formula_uni,
                   data = cd_subs,
                   prior = pr_uni,
                   stanvars =stanvars_nopred_uni,
                   warmup = wup,#may be necessary 
                   iter = wup+sam, #doesn't take a lot of runs
                   chains = 4, # 4 chains in parallel
                   cores = 4, # on 4 CPUs
                   control = list(adapt_delta = 0.90),
                   # init = 0, # could ease initialisation
                   # threads = 4, # on 4 CPUs
                   # open_progress = TRUE) # make a progress bar
                   refresh = 100, # echo chain progress every n iterations
                   silent = 1, # echo some Stan messages
                   backend = 'cmdstanr')
    
  }
)

### Bimodal -------------------------------------------------------------
#warning, takes nearly 15 minutes for 700 iterations!
system.time(
  {
    
    np_mix = brm(formula = formula_mix,
                   data = cd_subs,
                   prior = pr_mix,
                   stanvars =stanvars_nopred,
                   warmup = wup,#may be necessary 
                   iter = wup+sam, #doesn't take a lot of runs
                   chains = 4, # 4 chains in parallel
                   cores = 4, # on 4 CPUs
                   control = list(adapt_delta = 0.90),
                   # init = 0, # could ease initialisation
                   # threads = 4, # on 4 CPUs
                   # open_progress = TRUE) # make a progress bar
                   refresh = 100, # echo chain progress every n iterations
                   silent = 1, # echo some Stan messages
                   backend = 'cmdstanr')
    
  }
)


## Model comparison ------------------------------------------------------
#better predictions should justify fitting a mixture model 
loo_uni = loo(np_uni)
loo_mix = loo(np_mix)
lc_unimix = loo_compare(loo_uni, loo_mix)
print(lc_unimix)
#the mixture model has higher predictive power
pnorm(q = lc_unimix[2,1], sd = lc_unimix[2,2],lower.tail = FALSE)

## Plot coefficients -----------------------------------------------------

#check divergences
# nt_mix = nuts_params(np_mix)
# bayesplot::mcmc_scatter(full_mix,
#                       np = nt_mix,
#                       regex_pars = c('Intercept_theta1', '^b_fmu1_Intercept')
#                       )

#main effects means
#weighting (logistic scaled)
plot(np_mix,
     variable = 'Intercept_theta1',
     transform = plogis)
plot(np_mix,
     variable = '^b_theta1',
     regex = TRUE) 
#primary mu
plot(np_mix,
     variable = '^b_fmu1', # all effects have similar names
     regex = TRUE,
     nvariables = 4,
     transform = unwrap_circular_deg) 
#secondary mu
plot(np_mix,
     variable = '^b_fmu2', # all effects have similar names
     regex = TRUE,
     nvariables = 5,
     transform = unwrap_circular_deg) 
#primary kappa
plot(np_mix,
     variable = '^b_k1',
     regex = TRUE)#main effects means converge well
plot(np_mix,
     variable = '^sd_ID__theta1',
     regex = TRUE)
plot(np_mix,
     variable = '^sd_ID__k1',
     regex = TRUE)
plot(np_mix,
     variable = '^kappamu1',
     regex = TRUE)
plot(np_mix,
     variable = '^kappamu2',
     regex = TRUE)

#Many random effects of mu
plot(np_mix,
     variable = '^b_zmu1',
     regex = TRUE,
     ask = FALSE,
     transform = unwrap_circular_deg)
plot(np_mix,
     variable = '^b_zmu2',
     regex = TRUE,
     ask = FALSE,
     transform = unwrap_circular_deg)

## Summarise coefficients ------------------------------------------------
#extract the medians of all parameters
sm_mix = summary(np_mix, robust = TRUE)
print(sm_mix$fixed[c('theta1_Intercept',
                     'fmu1_Intercept',
                     'fmu2_Intercept',
                     'k1_Intercept'),],
      digits = 3)

UnwrapRhats(np_mix, variable = 'fmu')
UnwrapRhats(np_mix)
#find the names of the fixed effects
rn_sm_mix = with(sm_mix, rownames(fixed) )
#Investigate the 
sm_mix$spec_pars
print(sm_mix$fixed[grepl(pattern = 'mu', x = rn_sm_mix ),], digits = 3)
print(sm_mix$fixed[grepl(pattern = 'theta', x = rn_sm_mix ),], digits = 3)
print(sm_mix$fixed[grepl(pattern = 'k1', x = rn_sm_mix ),], digits = 3)



