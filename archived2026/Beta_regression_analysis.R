#FOR A 'CLEAN' RUN, PRESS ctrl+shift+F10 to RESTART Rstudio
graphics.off()
# Details ---------------------------------------------------------------
#       AUTHOR:	James Foster              DATE: 2024 04 03
#     MODIFIED:	James Foster              DATE: 2024 04 10
#
#  DESCRIPTION: Load dance angles, calculate mean vectors and fit beta regression
#               to the mean vectors via brms.
#               
#       INPUTS: 
#               
#      OUTPUTS: Plots and test statistics
#
#	   CHANGES: - Comparisons w/ emmeans
#             - Adapt delta selection for max random effects model
#             - Organise data by colour and brightness
#             - Set and inspect priors
#
#   REFERENCES: Bürkner, P.-C., 2018. 
#               Advanced Bayesian Multilevel Modeling with the R Package brms.
#               The R Journal 10, 395–411. 
#               https://doi.org/10.32614/RJ-2018-017
# 
#               Heiss, A. 2021. 
#               A Guide to Modeling Proportions with Bayesian Beta and 
#               Zero-Inflated Beta Regression Models. 
#               https://doi.org/10.59350/7p1a4-0tw75.
#
#    EXAMPLES:  
#
# 
#TODO   ---------------------------------------------
#TODO   
#- Load data  +
#- Extract results  +
#- Switch to marginal_effects  +
#- Narrower RE prior  +
#- Get modelling consistent +
#- Eliminate divergent transitions
#- Handwrite hypothesis tests for all contrasts
# Useful functions --------------------------------------------------------

# . Load packages ----------------------------------------------------------
#needs installing before first use (in Rstudio, see automatic message)
suppressMessages(#these are disturbing users unnecessarily
  {
    require(circular)#package for handling circular data
    require(CircStats)#package for circular hypothesis tests
    require(cmdstanr)#package for Bayesian modelling via Stan
    require(brms)#package for preparing Stan models
    require(extraDistr)#package for unusual distributions
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


# . Prior plotting --------------------------------------------------------

#plot normal priors
LinesNorm = function(ii, #this value 
                     sq, #the full sequence
                     pr, #prior parameters
                     lwd = 1, #line width
                     param = 'mu', #parameter type
                     xx = seq(from = 0, to = 1, length.out = 1e3), #sequence to plot across
                     ...)
{
  switch(EXPR = param,
         mu = 
           lines(x = xx,
                 y = dnorm(qlogis(xx), mean = ii, sd = 1),
                 col = rgb(0.5-ii/(2.1*max(sq)), 
                           0, 
                           0.5+ii/(2.1*max(sq)),
                           alpha = with(pr, dnorm(ii, mean = mu, sd = sd)/dnorm(0, sd = sd) )),
                 lwd = lwd,
                 ...
           ),
         # https://github.com/boboppie/kruschke-doing_bayesian_data_analysis/blob/master/2e/DBDA2E-utilities.R
         # betaABfromMeanKappa = function( mean , kappa ) {
         #   if ( mean <=0 | mean >= 1) stop("must have 0 < mean < 1")
         #   if ( kappa <=0 ) stop("kappa must be > 0")
         #   a = mean * kappa
         #   b = ( 1.0 - mean ) * kappa
         #   return( list( a=a , b=b ) )
         # }
         phi = 
           lines(x = xx,
                 y = dbeta(x = xx,
                           shape1 = plogis(0)*exp(ii),
                           shape2 = (1-plogis(0)) * exp(ii)),
                 col = rgb(0.5-ii/(2.1*max(sq)), 
                           0, 
                           0.5+ii/(2.1*max(sq)), 
                           alpha = with(pr, dnorm(ii, mean = mu, sd = sd)/dnorm(0, sd = sd)) ),
                 lwd = lwd,
                 ...
           )
  )
  
}


#plot student-T priors
LinesStuT = function(ii,#with value 
                     sq,  #the sequence of values
                     pr,#prior parameters
                     lwd = 1, #line width
                     param = 'mu', #parameter type
                     mu0 = 0, #Intercept mu parameter
                     sd0 = 1, #Intercept standard deviation of mu
                     phi0 = 1, #Intercept phi parameter 
                     xx = seq(from = 0, to = 1, length.out = 1e3),
                     ...)
{
  switch(EXPR = param,
         mu = 
           lines(x = xx,
                 y = dnorm(qlogis(xx), mean = mu0+ii, sd = sd0),
                 col = rgb(0.5-ii/(2.1*max(sq)), 
                           0, 
                           0.5+ii/(2.1*max(sq)),
                           alpha = with(pr, dstudent_t(x = ii,df = df,mu = mu, sigma =  sd)/
                                          dstudent_t(x = mu,df = df,mu = mu, sigma =  sd)) ),
                 lwd = lwd,
                 ...
           ),
         phi = 
           lines(x = xx,
                 y = dbeta(x = xx,
                           shape1 = plogis(mu0)*exp(phi0+ii),
                           shape2 = (1-plogis(mu0)) * exp(phi0+ii)),
                 col = rgb(0.5-ii/(2.1*max(sq)), 
                           0, 
                           0.5+ii/(2.1*max(sq)),
                           alpha = with(pr, dstudent_t(x = ii,df = df,mu = mu, sigma =  sd)/
                                          dstudent_t(x = mu,df = df,mu = mu, sigma =  sd)) ),
                 lwd = lwd,
                 ...
           )
  )
  
}



# . System parameters -----------------------------------------------------

#Check the operating system and assign a logical flag (T or F)
sys_win = Sys.info()[['sysname']] == 'Windows'
#User profile instead of home directory
if(sys_win){
  #get rid of all the backslashes
  ltp = gsub('\\\\', '/', Sys.getenv('USERPROFILE'))#Why does windows have to make this so difficult
}else{#Root directory should be the "HOME" directory on a Mac (or Linux?)
  ltp  =  Sys.getenv('HOME')#Life was easier on Mac
}

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

cd = read.table(file = path_file, 
                header = T, 
                sep  = ',')
View(cd)

cd = within(cd,
            {
              ID = as.factor(bee) # beedance identifier as a factor
              date = as.factor(date) # date as a factor
              signed_angle = Mod360.180(bearing)  # bearing between -180 and 180
              angle = circular(rad(signed_angle),# bearing between -pi and pi
                               rotation = 'clock') # circular format suppresses later warnings
            }
            )

#remove ul2
cd = subset(cd, 
            !(light_type %in% 'ul2') )

u_id = with(cd, unique(ID)) # unique beedances
length(u_id)#169 beedances

#check which individuals experienced both conditions (not necessary w/ brms?)
idx_both  = with(cd, 
               {
                 unique(subset(ID, light_type == 'ul')) %in% 
                   unique(subset(ID, light_type == 'uh'))
               }
              )
#extract the names
u_id_both = u_id[idx_both]


# Reorganise data ---------------------------------------------------------
cd = within(cd,
            {
            col_bright = do.call(what = rbind,
                                 args = strsplit(x = light_type,
                                                 split = '')
                                 )
            colour = col_bright[,1]
            brightn = col_bright[,2]
            }
            )



# Calculate mean vectors --------------------------------------------------


#calculate mean vectors
mean_vectors = aggregate(angle~ID*brightn*colour,
                         data = cd,
                         FUN = rho.circular
                        )
#correct names
mean_vectors = within(mean_vectors,
                       {mean_vector = angle; rm(angle)} # anlge now indicates a mean vector, not an angle
)
#plot mean vectors
boxplot(mean_vector~brightn*colour,
        data = mean_vectors,
        ylim = c(0,1),
        col = adjustcolor('blue', alpha.f = 0.2))
abline(h = c(0,1))
#calculate log kappa
mean_vectors = within(mean_vectors,
                       {
                         kappa <- circular::A1inv(mean_vector)
                         logkappa <- log(kappa)
                       }
)


# Fit a Bayesian Beta distributed model -----------------------------------


## Set up formulas --------------------------------------------------------


#set up model formula

### Simple formula with random effects on intercepts ---------------------


frm_bmod = bf(mean_vector~brightn*colour + # light types affect accuracy
                            (1|ID), # each individual is more or less oriented
              phi ~ brightn*colour + # light types affect variance in accuracy
                      (1|ID), # some individuals vary more in accuracy than others
              family = brms::Beta()) 

### Maximal formula with all random effects ---------------------
#set up model formula with indiv effects on light type
frm_bmod_maxRE = bf(mean_vector~brightn*colour + # light types affect accuracy
                            (1 + brightn*colour|ID), # each individual is more or less oriented
              phi ~ brightn*colour + # light types affect variance in accuracy
                      (1 + brightn*colour|ID), # some individuals vary more in accuracy than others
              family = brms::Beta()) 

### Simpler formula without interactions in random effects ---------------------
#set up model formula with indiv effects on light type
frm_bmod_nointerRE = bf(mean_vector~brightn*colour + # light types affect accuracy
                            (1 + brightn+colour|ID), # each individual is more or less oriented
              phi ~ brightn*colour + # light types affect variance in accuracy
                      (1 + brightn+colour|ID), # some individuals vary more in accuracy than others
              family = brms::Beta()) 


## Set priors ------------------------------------------------------------


#find out the priors needed
priors_bmod = get_prior(formula = frm_bmod, 
                        data = mean_vectors)

# unique(priors_bmod$class)
# [1] "b"         "Intercept" "sd" 

### Suggest priors -------------------------------------------------------

#suggest prior values
pr_mu_Intercept = list(dist = 'normal',
             mu = 0,
             sd = 1)
pr_phi_Intercept = list(dist = 'normal',
             mu = 1,
             sd = 1)
pr_mu_b = list(dist = 'normal',
             mu = 0,
             sd = 1)
pr_phi_b = list(dist = 'normal',
             mu = 0,
             sd = 1)
pr_mu_sd = list(dist = 'student_t',
             df = 3,
             mu = 0,
             sd = 0.5)
pr_phi_sd = list(dist = 'student_t',
                 df = 3,
                 mu = 0,
                 sd = 0.5)
#set up inspection sequences
#potential data sequence
xseq = seq(from = 0, 
           to = 1, 
           length.out = 1e3)
#normal parameter sequences
pr_mu_Intercept_seq = with(pr_mu_Intercept,
                           {
                             seq(from = qnorm(p =  0.05/2,mean = mu, sd = sd),
                                  to = qnorm(p = 1-0.05/2,mean = mu, sd = sd),
                                  length.out = 20)
                           }
                          )
pr_mu_b_seq = with(pr_mu_b,
                           {
                             seq(from = qnorm(p =  0.05/2,mean = mu, sd = sd),
                                  to = qnorm(p = 1-0.05/2,mean = mu, sd = sd),
                                  length.out = 20)
                           }
                          )
pr_phi_Intercept_seq = with(pr_phi_Intercept,
                           {
                             seq(from = qnorm(p =  0.05/2,mean = mu, sd = sd),
                                  to = qnorm(p = 1-0.05/2,mean = mu, sd = sd),
                                  length.out = 20)
                           }
                          )
pr_phi_b_seq = with(pr_phi_b,
                           {
                             seq(from = qnorm(p =  0.05/2,mean = mu, sd = sd),
                                  to = qnorm(p = 1-0.05/2,mean = mu, sd = sd),
                                  length.out = 20)
                           }
                          )
#t-distributed parameter sequence 
pr_mu_sd_seq = with(pr_mu_sd,
                   {
                     seq(from = qstudent_t(p =  0.05/2,df = df,mu = mu, sigma =  sd),
                         to = qstudent_t(p = 1-(0.05/2),df = df,mu = mu, sigma =  sd),
                         length.out = 20)
                   }
)
pr_phi_sd_seq = with(pr_phi_sd,
                   {
                     seq(from = qstudent_t(p =  0.05/2,df = df,mu = mu, sigma =  sd),
                         to = qstudent_t(p = 1-(0.05/2),df = df,mu = mu, sigma =  sd),
                         length.out = 20)
                   }
)


### Inspect priors -------------------------------------------------------

#plot mu expected range
plot(x = xseq,
     y = dnorm(qlogis(xseq), mean = 0, sd = 1),
     type = 'l',
     lwd = 2,
     xlab = 'Mu parameter',
     ylab = 'Prior probability density'
     )
invisible({
sapply(X = pr_mu_Intercept_seq,
       FUN = LinesNorm,
       pr = pr_mu_Intercept,
       param = 'mu',
       sq = pr_mu_Intercept_seq,
       lwd = 1)
})


#look at phi based on beta probability density


#aim for concave probability distributions, set mean phi close to 2 ≈ exp(1)
plot(x = xseq,
     y = dbeta(x = xseq, 
               shape1 = plogis(0)*exp(0+1), 
               shape2 = (1 - plogis(0)) * exp(0+1)),
     type = 'l',
     lwd = 2,
     xlab = 'Fitted mean',
     ylab = 'Prior probability density',
     ylim = c(0,4)
     )
invisible({
  sapply(X = pr_phi_Intercept_seq,
         FUN = LinesNorm,
         pr = pr_phi_Intercept,
         param = 'phi',
         sq = pr_phi_Intercept_seq,
         lwd = 1)
})


#look at t-distributed random effects for mu
plot(x = xseq,
     y = dnorm(qlogis(xseq), mean = 0, sd = 1),
     type = 'l',
     lwd = 2,
     xlab = 'Mu parameter',
     ylab = 'Prior probability density'
)
invisible({
  sapply(X = pr_mu_sd_seq,
         FUN = LinesStuT,
         pr = pr_mu_sd,
         param = 'mu',
         sq = pr_mu_sd_seq,
         mu0 = pr_mu_Intercept$mu,
         sd0 = pr_mu_Intercept$sd,
         lwd = 1)
})

#look at t-distributed random effects on phi
#aim for concave probability distributions, set mean phi close to 2 ≈ exp(1)
plot(x = xseq,
     y = dbeta(x = xseq, 
               shape1 = plogis(0)*exp(0+1), 
               shape2 = (1 - plogis(0)) * exp(0+1)),
     type = 'l',
     lwd = 2,
     xlab = 'Fitted mean',
     ylab = 'Prior probability density',
     ylim = c(0,3)
)
invisible({
  sapply(X = pr_phi_sd_seq,
         FUN = LinesStuT,
         pr = pr_phi_sd,
         param = 'phi',
         sq = pr_phi_sd_seq,
         mu0 = pr_mu_Intercept$mu,
         phi0 = pr_phi_Intercept$mu,
         lwd = 1)
})


### Assign priors --------------------------------------------------------


#set priors on:
# - all coefficients to a normal(0,1) = 95% probability of values 0.12–0.87
# - all grouping paramters to a student_t(df = 3, mu = 0, sigma = 2.5) = values < 5.88
# - all other parameters (?) to a normal(0,10)= 95% probability of values 0–1 +- 1e-7

pr_expression = expression({
  prior = sapply(class, 
                 switch,
                 #the general prior for coefficients (should centre on 0)
                 b = with(pr_mu_b,
                          {paste0(dist,'(',mu,',',sd,')')}),
                 #the general prior for standard deviations (should be long-tailed)
                 sd = with(pr_mu_sd,
                           {paste0(dist,'(',df,',',mu,',',sd,')')}),
                 cor = 'lkj(1)', # The only supported prior for correlation matrices is the 'lkj' prior
                 'normal(0,1)' # with this many params need more informative priors
  )
  #specific prior for phi coefficients (log scaled)
  prior = ifelse(test = class %in% 'b' &
                   dpar %in% 'phi',
                 yes = with(pr_phi_b,
                            {paste0(dist,'(',mu,',',sd,')')}),
                 no = prior
  )
  #specific prior for the intercept of phi (distribution is convex for phi > exp(1))
  prior = ifelse(test = class %in% 'Intercept' &
                   dpar %in% 'phi' |
                   class %in% 'b' & coef %in% '' &#this may also be the intercept condition
                   dpar %in% 'phi' ,
                 yes = with(pr_phi_Intercept,
                            {paste0(dist,'(',mu,',',sd,')')}),
                 no = prior
  )
})

# basic model
priors_bmod = within(priors_bmod, expr = eval(pr_expression))

#model with all random effects
priors_bmod_mxRE = within(get_prior(formula = frm_bmod_maxRE, 
                                    data = mean_vectors), 
                          expr = eval(pr_expression))

#model with random effects without interactions
priors_bmod_niRE= within(get_prior(formula = frm_bmod_nointerRE, 
                                   data = mean_vectors), 
                         expr = eval(pr_expression))


# Run the models -----------------------------------------------------------


## Rescale data ----------------------------------------------------------


with(mean_vectors, hist(mean_vector, breaks = 1e3, xlim = c(0.9, 1.0)))
abline(v = plogis(5.2), col = 'red')#cutoff for mu Intercept, avoid mu = 1?

#rescale mean vectors equal to 1
resc_mvs = within(mean_vectors,
                  {
                  mean_vector = ifelse(test = mean_vector == 1,
                                       yes = plogis(5.2),
                                       no = mean_vector)
                  }
                  )


## Check that stancode can be made ---------------------------------------
frm_list = list(frm_bmod, frm_bmod_nointerRE, frm_bmod_maxRE)
pr_list = list(priors_bmod, priors_bmod_niRE, priors_bmod_mxRE)
dt_list = list(resc_mvs,resc_mvs,resc_mvs)
fm_list = list(Beta(link = 'logit',
                    link_phi = 'log'),
               Beta(link = 'logit',
                   link_phi = 'log'),
               Beta(link = 'logit',
                link_phi = 'log')
               )

sc_test = mapply(FUN = make_stancode,
                 formula = frm_list,
                 prior = pr_list,
                 data = dt_list,
                 family = fm_list
                )

## Check that priors are appropriate -------------------------------------

system.time( #takes about 5 minutes
  {
pr_test = mapply(FUN = brm,
                 SIMPLIFY = FALSE,
                 sample_prior = 'only', #sample only the prior
                 formula = frm_list,
                 prior = pr_list,
                 data = dt_list,
                 family = fm_list,
                 iter = 1e3, # only small number necessary
                 control = list( adapt_delta = 0.95 ), # closer to 1.0 means higher resolution sampling
                 cores = 4,
                 backend = 'cmdstanr',
                 silent = 2, #don't print lots of iteration information
                 refresh = 0
                )
  }
)

lapply(X = pr_test,
       FUN = conditional_effects,
       spaghetti = TRUE,
       ask = FALSE)

## Run the models --------------------------------------------------------
system.time( #takes about 10 minutes
  {
md_list = mapply(FUN = brm,
                 SIMPLIFY = FALSE,
                 formula = frm_list,
                 prior = pr_list,
                 data = dt_list,
                 family = fm_list,
                 iter = 2e3, # only small number necessary
                 control = list( adapt_delta = 0.99 ), # closer to 1.0 means higher resolution sampling
                 cores = 4,
                 backend = 'cmdstanr',
                 silent = 2, #don't print lots of iteration information
                 refresh = 0
)
  }
)

#I do not know where these warnings are coming from, there must be a rounding error
## Exception: beta_lpdf: Second shape parameter is 0
#  Probably related:
## 1 of 4 chains had an E-BFMI less than 0.2.

    # ### Simpler model --------------------------------------------------------
    # #run and time the model
    # system.time(expr = 
    #               {
    #                 bmod_test = brm(formula = frm_bmod ,
    #                                 data = resc_mvs,
    #                                 iter = 5e3, # about 2500 iterations/ minute
    #                                 family = Beta(),
    #                                 prior = priors_bmod,#use default
    #                                 control = list( adapt_delta = 0.95 ), # closer to 1.0 means higher resolution sampling
    #                                 cores = 4,
    #                                 backend = 'cmdstanr',
    #                                 init = 1,
    #                                 silent = 2, #don't print lots of iteration information
    #                                 refresh = 0,
    #                 )
    #               }
    # )
    # #takes about 3 minutes
    # 
    # #run and time the model
    # bmod_test_maxRE_pr = brm(formula = frm_bmod_maxRE ,
    #                       data = resc_mvs,
    #                       iter = 1e3, # about 95 iterations/ minute
    #                       sample_prior = 'only', 
    #                       family = Beta(link = 'logit',
    #                                     link_phi = 'log'),
    #                       prior = priors_bmod_mxRE,#use default
    #                       control = list( adapt_delta = 0.99), # closer to 1.0 means higher resolution sampling
    #                       cores = 4,
    #                       backend = 'cmdstanr',
    #                       init = '0',
    #                       silent = 2, #don't print lots of iteration information
    #                       refresh = 0
    # )
    # 
    # 
    # system.time(expr = 
    #               {
    #                 bmod_test_maxRE = brm(formula = frm_bmod_maxRE ,
    #                                 data = resc_mvs,
    #                                 iter = 2e3, # about 95 iterations/ minute
    #                                 family = Beta(link = 'logit',
    #                                               link_phi = 'log'),
    #                                 prior = priors_bmod_mxRE,#use default
    #                                 control = list( adapt_delta = 0.99), # closer to 1.0 means higher resolution sampling
    #                                 cores = 4,
    #                                 backend = 'cmdstanr',
    #                                 init = '0',
    #                                 silent = 2, #don't print lots of iteration information
    #                                 refresh = 0
    #                 )
    #               }
    # )
    # #takes about 7 minutes
    # #Warning: 3 of 4 chains had an E-BFMI less than 0.2.
    # # Chain 2 Informational Message: The current Metropolis proposal is about to be rejected because of the following issue:
    # # Chain 2 Exception: beta_lpdf: Second shape parameter is 0, but must be positive finite! 
    # 
    # 
    # 
    # #run and time the model
    # system.time(expr = 
    #               {
    #                 bmod_test_nointerRE = brm(formula = frm_bmod_nointerRE ,
    #                                       data = resc_mvs,
    #                                       iter = 2e3, # about 95 iterations/ minute
    #                                       family = Beta(link = 'logit',
    #                                                     link_phi = 'log'),
    #                                       prior = priors_bmod_niRE,#use default
    #                                       control = list( adapt_delta = 0.99), # closer to 1.0 means higher resolution sampling
    #                                       cores = 4,
    #                                       backend = 'cmdstanr',
    #                                       init = '0',
    #                                       silent = 2, #don't print lots of iteration information
    #                                       refresh = 0
    #                 )
    #               }
    # )

# Check model convergence -------------------------------------------------
invisible({
lapply(X = md_list,
       FUN = plot,
       ask = FALSE)

sm_list = lapply(X = md_list,
       FUN = summary)
})
print(sm_list)
    # 
    # plot(bmod_test, ask = FALSE)
    # #
    # #
    # #
    # 
    # summary(bmod_test)
    # 
    # plot(bmod_test_nointerRE, ask = FALSE)
    # #
    # #
    # #
    # #
    # #
    # #
    # 
    # summary(bmod_test_nointerRE)
    # 
    # plot(bmod_test_maxRE, ask = FALSE)
    # #
    # #
    # #
    # #
    # #
    # #
    # 
    # summary(bmod_test_maxRE)

# Compare models ----------------------------------------------------------
names(md_list) = c('REintercepts','REslopes','REinteractions')
loo_list = lapply(X = md_list,
                  FUN = loo)
bmod_lc = loo_compare(loo_list)
print(bmod_lc)
    # bmod_loo = loo(bmod_test)
    # bmod_loo_maxRE = loo(bmod_test_maxRE)
    # bmod_loo_nointerRE = loo(bmod_test_nointerRE)
    # bmod_lc = loo_compare(bmod_loo,bmod_loo_nointerRE, bmod_loo_maxRE)
    # print(bmod_lc)

    # mod_chosen = get(x = dimnames(bmod_lc)[[1]][1])
mod_chosen = with(md_list, 
                  {get(x = dimnames(bmod_lc)[[1]][1])} #the model with the highest predictive density
                  )

# Extract model predictions -----------------------------------------------

cond_eff = conditional_effects(mod_chosen)
plot(x = cond_eff, 
     points = TRUE, 
     point_args = list(width = 0.3, 
                       col = adjustcolor('blue',
                                         alpha.f = 0.2),
                       width = 0.3
                       ),
     ask = FALSE
    )

# Calculate model contrasts -----------------------------------------------


  ##. Marginal effects version -
  # require(marginaleffects)
  # marginaleffects::avg_comparisons(mod_chosen,
  #                                  hypothesis = 'pairwise')


# . emmeans version -------------------------------------------------------
require(emmeans)
em = emmeans::emmeans(mod_chosen,
                      specs = pairwise~brightn*colour,
                      regrid = "response" )
em_con = emmeans::contrast(em,
                           method = "revpairwise")
em_sum = summary(em_con)
print(within(em_sum, 
             {
               difference = ifelse(sign(lower.HPD) == sign(upper.HPD),
                                          yes = '*',#estimated contrast does not overlap with zero
                                          no = '') # estimated contrast does overlap with zero difference
             }
             ) 
      )

##  contrast estimate lower.HPD upper.HPD difference
## l g - h g  -0.2811   -0.3599   -0.2117 *         
## h u - h g   0.0420    0.0141    0.0722 *         
## h u - l g   0.3237    0.2526    0.3943 *         
## l u - l g   0.0343   -0.0432    0.1120           
## l u - h g  -0.2481   -0.2955   -0.1983 *         
## l u - h u  -0.2899   -0.3288   -0.2561 *    

# . Handwritten version ---------------------------------------------------


all_draws = prepare_predictions(mod_chosen)
#collect coefficients
bmod_qnt = fixef(mod_chosen, robust = TRUE)
bmod_coef = cbind(all_draws$dpars$mu$fe$b,
                  all_draws$dpars$phi$fe$b)

#collect and estimate p values from samples
empirical.p = function(x, tails = 2){ifelse(test = median(x)>0,
                                            yes = mean(x<0)*tails, 
                                            no = mean(x>0)*tails)
}

#collect and estimate standard error from samples
sErr = function(x){sd(x) / sqrt(sum(!is.na(x)))}

#collect and estimate z scores from samples
LOR.z = function(lor, Efun = sErr){median(lor) / Efun(lor)}

#collect and estimate z values from samples
P_z = function(z, tails = 2)
{tails*(1-pnorm(abs(z)))}

#two tailed p-values where CI do not overlap with 0 (more conservative)
p_coef2 = apply(X = bmod_coef,
                MARGIN = 2,
                FUN = empirical.p
)

q_coef = apply(X = bmod_coef,
               MARGIN = 2,
               FUN = quantile,
               probs = c(0,1)+0.5*c(1,-1)*0.05
)
names(p_coef2) = sub(names(p_coef2), pattern = 'b_', replacement = '')

z_coef = apply(X = bmod_coef,
               MARGIN = 2,
               FUN = LOR.z,
               Efun = sd # for this sample size, sd would make more sense
)

p_gr8r_z = sapply(X = z_coef,
                  FUN = P_z,
                  tails = 2 # more conservative
)

# Report contrasts --------------------------------------------------------
#collect in data frame with original estimates
test_statistics = data.frame(coefficient = names(p_coef2),
                             estimate = unlist(bmod_qnt[names(p_coef2), 'Estimate']),
                             lower_quantile = unlist(bmod_qnt[names(p_coef2),'Q2.5']),
                             upper_quantile = unlist(bmod_qnt[names(p_coef2),'Q97.5']),
                             `prop(H0 draws) 2-tailed` = p_coef2,
                             `two-tailed signif` = ifelse(p_coef2 <0.05, yes = '*', no = ''),
                             `non-zero CI` = ifelse(
                               apply(X = bmod_qnt[names(p_coef2),c('Q2.5','Q97.5')],
                                     MARGIN = 1,
                                     FUN = function(x){abs(sum(sign(unlist(x))))}), 
                               yes = '*', 
                               no = ''),
                             `Odds-ratio` = ifelse(test = bmod_qnt[names(p_coef2), 'Estimate'] >=0,
                                                   yes = exp(unlist(bmod_qnt[names(p_coef2), 'Estimate'])),
                                                   no = -1/exp(unlist(bmod_qnt[names(p_coef2), 'Estimate']))
                             ),
                             `OR-lower` = ifelse(test = bmod_qnt[names(p_coef2), 'Q2.5'] >=0,
                                                 yes = exp(unlist(bmod_qnt[names(p_coef2), 'Q2.5'])),
                                                 no = -1/exp(unlist(bmod_qnt[names(p_coef2), 'Q2.5']))
                             ),
                             `OR-upper` = ifelse(test = bmod_qnt[names(p_coef2), 'Q97.5'] >=0,
                                                 yes = exp(unlist(bmod_qnt[names(p_coef2), 'Q97.5'])),
                                                 no = -1/exp(unlist(bmod_qnt[names(p_coef2), 'Q97.5']))
                             ),
                             `z-score` = z_coef,
                             `p-greater-z` = p_gr8r_z
)

#print the results for the user                     
with(test_statistics,
     {
       print(
         data.frame(coefficient = coefficient,
                   z = round(unlist(z.score), 3),
                    p = round(unlist(p.greater.z), 4)
                   )
       )
     }
)

##           coefficient      z      p
## 1            Intercept 12.854 0.0000
## 2             brightnl -7.989 0.0000
## 3              colouru  3.537 0.0004
## 4     brightnl:colouru -2.055 0.0399
## 5        phi_Intercept  7.426 0.0000
## 6         phi_brightnl -2.122 0.0338
## 7          phi_colouru  5.939 0.0000
## 8 phi_brightnl:colouru -0.400 0.6890