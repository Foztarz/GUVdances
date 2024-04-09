#FOR A 'CLEAN' RUN, PRESS ctrl+shift+F10 to RESTART Rstudio
graphics.off()
# Details ---------------------------------------------------------------
#       AUTHOR:	James Foster              DATE: 2024 04 03
#     MODIFIED:	James Foster              DATE: 2024 04 03
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
#             - 
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
#- Narrower RE prior
#- Eliminate divergent transitions
#- Handwrite hypothesis tests for all contrasts
#- Get modelling consistent 
# Useful functions --------------------------------------------------------

# . Load packages ----------------------------------------------------------
#needs installing before first use (in Rstudio, see automatic message)
suppressMessages(#these are disturbing users unnecessarily
  {
    require(circular)#package for handling circular data
    require(CircStats)#package for circular hypothesis tests
    require(cmdstanr)#package for Bayesian modelling via Stan
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

#set up model formula
frm_bmod = bf(mean_vector~brightn*colour + # light types affect accuracy
                            (1|ID), # each individual is more or less oriented
              phi ~ brightn*colour + # light types affect variance in accuracy
                      (1|ID), # some individuals vary more in accuracy than others
              family = brms::Beta()) 
#set up model formula with indiv effects on light type
frm_bmod_maxRE = bf(mean_vector~brightn*colour + # light types affect accuracy
                            (1 + brightn*colour|ID), # each individual is more or less oriented
              phi ~ brightn*colour + # light types affect variance in accuracy
                      (1 + brightn*colour|ID), # some individuals vary more in accuracy than others
              family = brms::Beta()) 

#find out the priors needed
priors_bmod = get_prior(formula = frm_bmod, 
                        data = mean_vectors)

#set priors on:
# - all coefficients to a normal(0,1) = 95% probability of values 0.12–0.87
# - all grouping paramters to a student_t(df = 3, mu = 0, sigma = 2.5) = values < 5.88
# - all other parameters (?) to a normal(0,10)= 95% probability of values 0–1 +- 1e-7
priors_bmod = within(priors_bmod, 
                     {
                       prior = sapply(class, 
                                      switch,
                                      b = 'normal(0,1)',
                                      sd = 'student_t(3, 0, 2.5)',
                                      'normal(0,10)'
                       )
                     }
)
priors_bmod_mxRE= within(get_prior(formula = frm_bmod_maxRE, 
                                    data = mean_vectors), 
                     {
                       prior = sapply(class, 
                                      switch,
                                      b = 'normal(0,1)',
                                      sd = 'student_t(3, 0, 2.5)',
                                      cor = 'lkj(1)', # The only supported prior for correlation matrices is the 'lkj' prior
                                      'normal(0,3)' # with this many params need more informative priors
                       )
                     }
)


# Run the models -----------------------------------------------------------

#rescale mean vectors equal to 1
resc_mvs = within(mean_vectors,
                  {
                  mean_vector = ifelse(test = mean_vector == 1,
                                       yes = mean_vector - 1e-8,
                                       no = mean_vector)
                  }
                  )

#run and time the model
system.time(expr = 
              {
                bmod_test = brm(formula = frm_bmod ,
                                data = resc_mvs,
                                iter = 5e3, # about 2500 iterations/ minute
                                family = Beta(),
                                prior = priors_bmod,#use default
                                control = list( adapt_delta = 0.95 ), # closer to 1.0 means higher resolution sampling
                                cores = 4,
                                backend = 'cmdstanr',
                                init = '0',
                                silent = 2, #don't print lots of iteration information
                                refresh = 0
                )
              }
)
#takes about 5 minutes

#run and time the model
system.time(expr = 
              {
                bmod_test_maxRE = brm(formula = frm_bmod_maxRE ,
                                data = resc_mvs,
                                iter = 1e3, # about 95 iterations/ minute
                                family = Beta(),
                                prior = priors_bmod_mxRE,#use default
                                control = list( adapt_delta = 0.95), # closer to 1.0 means higher resolution sampling
                                cores = 4,
                                backend = 'cmdstanr',
                                init = '0',
                                silent = 2, #don't print lots of iteration information
                                refresh = 0
                )
              }
)
#takes about 15 minutes
#Warning: 3 of 4 chains had an E-BFMI less than 0.2.


# Check model convergence -------------------------------------------------

plot(bmod_test, ask = FALSE)
#
#
#

summary(bmod_test)

plot(bmod_test_maxRE, ask = FALSE)
#
#
#
#
#
#

summary(bmod_test_maxRE)

# Compare models ----------------------------------------------------------

bmod_loo = loo(bmod_test)
bmod_loo_maxRE = loo(bmod_test_maxRE)
# bmod_looic = bmod_loo$estimates[3,]
# bmod_looic_maxRE = bmod_loo_maxRE$estimates[3,]

# all_ic = rbind(bmod_looic, bmod_looic_maxRE)
# print(all_ic[order(all_ic[,2],decreasing = TRUE),])#bmod is better, but not by >1SE
bmod_lc = loo_compare(bmod_loo, bmod_loo_maxRE)
print(bmod_lc)

mod_chosen = get(x = dimnames(bmod_lc)[[1]][1])
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

##    contrast estimate lower.HPD upper.HPD difference
##  l g - h g -0.22378  -0.32262   -0.1694 *         
##  h u - h g  0.00369  -0.00853    0.0274           
##  h u - l g  0.23129   0.16242    0.3277 *         
##  l u - h g -0.28590  -0.32270   -0.2396 *         
##  l u - l g -0.05479  -0.12318    0.0421           
##  l u - h u -0.29176  -0.31901   -0.2460 *   

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
## 1            Intercept 16.888 0.0000
## 2             brightnl -9.718 0.0000
## 3              colouru  0.419 0.6750
## 4     brightnl:colouru -1.579 0.1144
## 5        phi_Intercept  6.462 0.0000
## 6         phi_brightnl  0.381 0.7030
## 7          phi_colouru  0.581 0.5611
## 8 phi_brightnl:colouru  0.292 0.7705