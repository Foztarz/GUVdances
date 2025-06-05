#FOR A 'CLEAN' RUN, PRESS ctrl+shift+F10 to RESTART Rstudio
graphics.off()
# Details ---------------------------------------------------------------
#       AUTHOR:	James Foster              DATE: 2025 05 30
#     MODIFIED:	James Foster              DATE: 2025 06 03
#
#  DESCRIPTION: Load data and plot
#               
#       INPUTS: 
#               
#      OUTPUTS: csv
#
#	   CHANGES: - Using Beta(mu, phi) distribution to fit
#             - Comparisons of light sources
#
#   REFERENCES: Edrich, W., Neumeyer, C. and von Helversen, O. (1979).
#               “Anti-sun orientation” of bees with regard to a field of ultraviolet light. 
#               J. Comp. Physiol. 134, 151–157.
#
#               Stavenga, D. G. (2010).
#               On visual pigment templates and the spectral shape of invertebrate rhodopsins and metarhodopsins.
#               J. Comp. Physiol. A 196, 869–878.
#
#               Kühn, N. (2017).
#               Lichtintensität wird der Wellenlänge als Orientierungspunkt
#               während des Schwänzeltanzes der Bienen vorgezogen
#               Masters Thesis, Philipps Universität Marburg
#
#    EXAMPLES:  
#
# 
#TODO   ---------------------------------------------
#TODO   
#- Load data  +
#- Plot       +
#- Fit curves +
#- Switch to beta distribution +
#- Comment




## System parameters -----------------------------------------------------

#Check the operating system and assign a logical flag (T or F)
sys_win = Sys.info()[['sysname']] == 'Windows'
#User profile instead of home directory
if(sys_win){
  #get rid of all the backslashes
  ltp = gsub('\\\\', '/', Sys.getenv('USERPROFILE'))#Why does windows have to make this so difficult
}else{#Root directory should be the "HOME" directory on a Mac (or Linux?)
  ltp  =  Sys.getenv('HOME')#Life was easier on Mac
}

## Select file ---------------------------------------------------------

# set path to file
if(sys_win){#choose.files is only available on Windows
  message('\n\nPlease select the "Edrich1979_brightness-accuracy.csv" file\n\n')
  Sys.sleep(0.5)#goes too fast for the user to see the message on some computers
  path_file  <- choose.files(
    default = file.path(ltp,'Documents', "*.csv"),#For some reason this is not possible in the "root" user
    caption = 'Please select the "Edrich1979_brightness-accuracy.csv" file'
  )
}else{
  message('\n\nPlease select the ".Edrich1979_brightness-accuracy.csv" file\n\n')
  Sys.sleep(0.5)#goes too fast for the user to see the message on some computers
  path_file <- file.choose(new=F)
}
#show the user the path they have selected
if(is.null(path_file) | !length(path_file))
{stop('No file selected.')}else
{print(path_file)}


# Read in the data and format ---------------------------------------------

ed = read.table(file = path_file, 
                header = T, 
                sep  = ',')
View(ed)

with(ed,
     {
     plot(x = log10(intensity),
          y = accuracy,
          xlim = c(9,14),
          ylim = c(0,1),
          pch = 21,
          cex = 1.5,
          col = 'black',
          bg = sapply(paste(wavelength),
                       FUN = switch,
                       `354` = 'purple',
                       `429` = 'blue',
                       `535` = 'darkgreen',
                       'red')
          )
     abline(h = c(0,1))
     }
     )


# Fit models --------------------------------------------------------------
#prepare for model fitting
ed = within(ed,
            {
              logit_accuracy = qlogis(accuracy)
              log10_intensity = log10(intensity)
              colour = sapply(paste(wavelength),
                              FUN = switch,
                              `354` = 'UV',
                              `429` = 'blue',
                              `535` = '0green',
                              NA)
            }
            )


# Subset data by wavelength -----------------------------------------------

ed_uv = subset(ed, wavelength <400)
ed_bl = subset(ed, wavelength == 429)
ed_gr = subset(ed, wavelength >500)



#fit simple GLM
lm_uv = with(ed_uv,
             lm(logit_accuracy ~ log10_intensity) 
             )

lm_bl = with(ed_bl,
             lm(logit_accuracy ~ log10_intensity) 
             )

lm_gr = with(ed_gr,
             lm(logit_accuracy ~ log10_intensity) 
             )

xl = seq(from = 9,
         to = 14,
         length.out = 1e3)

nd = data.frame(log10_intensity = xl)
nd = merge(x = ed, y = nd, all.y = TRUE)

par(mfrow = c(1,1),
    mar = c(4.5, 4.5, 4.5, 4.5))
with(ed,
     {
       plot(x = log10_intensity,
            y = logit_accuracy,
            xlim = c(9,14),
            ylim = c(-5,5),
            pch = 21,
            cex = 1.5,
            col = 'black',
            bg = sapply(paste(wavelength),
                        FUN = switch,
                        `354` = 'purple',
                        `429` = 'blue',
                        `535` = 'darkgreen',
                        'red')
       )
       abline(h = c(0,1))
       abline(lm_uv,
              col = 'magenta',
              lwd = 3)
       abline(lm_bl,
              col = 'cyan',
              lwd = 3)
       abline(lm_gr,
              col = 'green',
              lwd = 3)
     }
)

with(ed,
     {
       plot(x = log10_intensity,
            y = accuracy,
            xlim = c(9,14),
            ylim = c(0,1),
            pch = 21,
            cex = 1.5,
            col = 'black',
            bg = sapply(paste(wavelength),
                        FUN = switch,
                        `354` = 'purple',
                        `429` = 'blue',
                        `535` = 'darkgreen',
                        'red')
       )
       abline(h = c(0,1))
       lines(x = xl,
             y = plogis(xl*coef(lm_uv)[2]+
                          coef(lm_uv)[1]),
              col = 'magenta',
              lwd = 3)
       lines(x = xl,
             y = plogis(xl*coef(lm_bl)[2]+
                          coef(lm_bl)[1]),
              col = 'cyan',
              lwd = 3)
       lines(x = xl,
             y = plogis(xl*coef(lm_gr)[2]+
                          coef(lm_gr)[1]),
              col = 'green',
              lwd = 3)
     }
)


#rough estimate psychometric
ed_uv = within(ed_uv,
               {
               rel_logit_acc = qlogis(
                               (accuracy - min(accuracy))/
                                max(accuracy - min(accuracy))
                                     )
               }
              )

ed_bl = within(ed_bl,
               {
               rel_logit_acc = qlogis(
                               (accuracy - min(accuracy))/
                                max(accuracy - min(accuracy))
                                     )
               }
              )

ed_gr = within(ed_gr,
               {
               rel_logit_acc = qlogis(
                               (accuracy - min(accuracy))/
                                max(accuracy - min(accuracy))
                                     )
               }
              )

relm_uv = with(subset(ed_uv, is.finite(rel_logit_acc)),
             lm(rel_logit_acc ~ log10_intensity) 
)

relm_bl = with(subset(ed_bl, is.finite(rel_logit_acc)),
             lm(rel_logit_acc ~ log10_intensity) 
)

relm_gr = with(subset(ed_gr, is.finite(rel_logit_acc)),
             lm(rel_logit_acc ~ log10_intensity) 
)



with(ed,
     {
       plot(x = log10_intensity,
            y = accuracy,
            xlim = c(9,14),
            ylim = c(0,1),
            pch = 21,
            cex = 1.5,
            col = 'black',
            bg = sapply(paste(wavelength),
                        FUN = switch,
                        `354` = 'purple',
                        `429` = 'blue',
                        `535` = 'darkgreen',
                        'red')
       )
       abline(h = c(0,1))
     }
)
with(ed_uv,
     {
      lines(x = xl,
           y = diff(range(accuracy))*
             plogis(xl*coef(relm_uv)[2]+
                        coef(relm_uv)[1]) +
             min(accuracy),
           col = 'magenta',
           lwd = 3)
     }
)
with(ed_bl,
     {
      lines(x = xl,
           y = diff(range(accuracy))*
             plogis(xl*coef(relm_bl)[2]+
                        coef(relm_bl)[1]) +
             min(accuracy),
           col = 'cyan',
           lwd = 3)
     }
)
with(ed_gr,
     {
      lines(x = xl,
           y = diff(range(accuracy))*
             plogis(xl*coef(relm_gr)[2]+
                        coef(relm_gr)[1]) +
             min(accuracy),
           col = 'green',
           lwd = 3)
     }
)



# Psychometric version ----------------------------------------------------

#load packages
require(cmdstanr)
require(brms)

#set up model fit
formula_nl = bf(
  #set up a formula for the curve as a whole,
  #including parameters found in the data (correct_incorrect, stimulus)
  #and parameters that we wish to estimate (baseline, lapse rate, inflection point, width).
  #Most of these are subject to further fixed (type) and random (animal) effects,
  #these need to be defined for each parameter.
  
  #Two parameters need special transformations
  #To keep the output between 0 and 1, additional effects of lapse rate
  #will be added on the "logit" scale (Lapse = inv_logit(LogitLapse)).
  #To avoid curve widths of 0, we can assume a positive slope (≥0)
  #and add additional effects to width on a log scale (Width = exp(LogWidth))
  formula = accuracy ~ 
    inv_logit(LogitBase) + (1 - inv_logit(LogitLapse) - inv_logit(LogitBase) ) *#curve region
    inv_logit( 4.39*(log10_intensity - Inflex) / exp(LogWidth) ) , #inflection-width curve
  # for each of these parameters, we can set up a separate formula 
  # that describes how to predict them from the data 
  #Base rate of correct choices: "
  LogitBase ~ 1, #Base rate of correct choices: "~ 1" gives the instruction "estimate the mean across all data"
  #Lapse rate on a log(odds) scale:
  LogitLapse ~ colour, #this is similar to the formula in our LMM example
  #inflection point of the initial curve:
  Inflex ~ colour, #N.B. this is similar to the intercept, so it does not include effects of stimulus level
  #log 80% width of the curve:
  LogWidth ~ colour, #N.B. this is similar to the slope, so all of its effects depend on stimulus level
  # family = bernoulli("identity"),
  nl = TRUE)#the joint distribution for these parameters is undefined, and therefore the parameters themselves are "nonlinear"

# #N.B. for plotting we will need to also definte "inv_logit" in our the R environment
inv_logit = inv_logit_scaled # this can be found in brms, but with an additional specifier

#set up priors
prior_nl = get_prior(formula = formula_nl,
                     data = ed)

print(
  with(prior_nl,
       paste0(class, '_', nlpar, '_', coef))
  )

# . . Base rate prior -----------------------------------------------------
#for the baseline, we will use a logit-normal distribution with a bias towards 0.5

#inspect the prior distribution

#set the prior distribution
prior_nl = within(prior_nl, 
                  { prior[
                    class %in% 'b' & #just the fixed effects
                      nlpar %in% 'LogitBase' &
                      coef %in% 'Intercept' 
                  ] = 'normal(0,1)' #a normal distribution centred on plogis(0) = 0.5
                  })
#this prior is automatically bounded between 0 and 1, 
# . . Lapse rate priors ---------------------------------------------------

# set the prior distribution
prior_nl = within(prior_nl, 
                  { prior[
                    class %in% 'b' & #just the fixed effects
                      nlpar %in% 'LogitLapse' &
                      coef %in% 'Intercept' 
                  ] = 'normal(-3,3)' #a normal distribution centred on -3
                  })

#for all other fixed effects on lapse rate, we'll suggest values around 0 (no effect)
prior_nl = within(prior_nl, 
                  { prior[
                    class %in% 'b' & #just the fixed effects
                      nlpar %in% 'LogitLapse' &
                      coef %in% 'colourblue' 
                  ] = 'normal(0,3)' #a normal distribution:mean 0, sd 3
                  })
prior_nl = within(prior_nl, 
                  { prior[
                    class %in% 'b' & #just the fixed effects
                      nlpar %in% 'LogitLapse' &
                      coef %in% 'colourUV' 
                  ] = 'normal(0,3)' #a normal distribution:mean 0, sd 3
                  })

# . . Inflection point priors ---------------------------------------------

#set the priors
prior_nl = within(prior_nl, 
                  { prior[
                    class %in% 'b' & #just the fixed effects
                      nlpar %in% 'Inflex' &
                      coef %in% 'Intercept'
                  ] = 'normal(3,3)' #a normal distribution:mean 3, sd 3
                  })
#for all coefficients, we'll suggest values around 0 (no effect)
prior_nl = within(prior_nl, 
                  { prior[
                    class %in% 'b' & #just the fixed effects
                      nlpar %in% 'Inflex' &
                      coef %in% 'colourblue' 
                  ] = 'normal(0,3)' #a normal distribution:mean 0, sd 3
                  })

prior_nl = within(prior_nl, 
                  { prior[
                    class %in% 'b' & #just the fixed effects
                      nlpar %in% 'Inflex' &
                      coef %in% 'colourUV' 
                  ] = 'normal(0,3)' #a normal distribution:mean 0, sd 3
                  })

# . . Rise region width priors --------------------------------------------

#set prior distribution
prior_nl = within(prior_nl, 
                  { prior[
                    class %in% 'b' & #just the fixed effects
                      nlpar %in% 'LogWidth' &
                      coef %in% 'Intercept' 
                  ] = 'normal(1,3)' #a normal distribution:mean 3, sd 3
                  })
#for all coefficients, we'll suggest values around 0 (no effect)
prior_nl = within(prior_nl, 
                  { prior[
                    class %in% 'b' & #just the fixed effects
                      nlpar %in% 'LogWidth' &
                      coef %in% 'colourblue'
                  ] = 'normal(0,3)' #a normal distribution:mean 0, sd 3
                  })
prior_nl = within(prior_nl, 
                  { prior[
                    class %in% 'b' & #just the fixed effects
                      nlpar %in% 'LogWidth' &
                      coef %in% 'colourUV'
                  ] = 'normal(0,3)' #a normal distribution:mean 0, sd 3
                  })
#inspect assigned priors
print(prior_nl)
##                 prior class       coef group resp dpar      nlpar lb ub  source
##  student_t(3, 0, 2.5) sigma                                        0    default
##                (flat)     b                                Inflex       default
##           normal(0,3)     b colourblue                     Inflex       default
##           normal(0,3)     b   colourUV                     Inflex       default
##           normal(3,3)     b  Intercept                     Inflex       default
##                (flat)     b                             LogitBase       default
##           normal(0,1)     b  Intercept                  LogitBase       default
##           (flat)     b                            LogitLapse       default
##           normal(0,3)     b colourblue                 LogitLapse       default
##           normal(0,3)     b   colourUV                 LogitLapse       default
##           normal(-3,3)     b  Intercept                 LogitLapse       default
##                (flat)     b                              LogWidth       default
##           normal(0,3)     b colourblue                   LogWidth       default
##           normal(0,3)     b   colourUV                   LogWidth       default
##           normal(1,3)     b  Intercept                   LogWidth       default

#double check that the prior distribution is viable by first setting up a short dummy run
# Dummy run
system.time(
  {
    dummy_fit = brm( formula = formula_nl, # using our nonlinear formula
                     data = ed, # our data
                     prior = prior_nl, # our priors 
                     sample_prior = 'only', #ignore the data to check the influence of the priors
                     iter = 300, # short run for 300 iterations
                     chains = 4, # 4 chains in parallel
                     cores = 4, # on 4 CPUs
                     refresh = 0, # don't echo chain progress
                     control = list(adapt_delta = 0.9), #finer sampling
                     backend = 'cmdstanr') # use cmdstanr (other compilers broken)
  }
)
# On my computer this takes <60s, each chain running for <1 seconds (mainly compile time)


#the default plot shows the values estimated for each parameter
# in each chain for each iteration
#fixed effects
plot(dummy_fit, 
     nvariables = 10,
     variable = "^b_", 
     regex = TRUE)


#double check that the prior distribution is viable by first setting up a short dummy run
# Dummy run
system.time(
  {
    nl_fit = brm( formula = formula_nl, # using our nonlinear formula
                     data = ed, # our data
                     prior = prior_nl, # our priors 
                     iter = 4000, # long run for 2000 iterations
                     chains = 4, # 4 chains in parallel
                     cores = 4, # on 4 CPUs
                     refresh = 0, # don't echo chain progress
                     control = list(adapt_delta = 0.9), #finer sampling
                     backend = 'cmdstanr') # use cmdstanr (other compilers broken)
  }
)
# On my computer this takes <60s, each chain running for <1 seconds (mainly compile time)


#the default plot shows the values estimated for each parameter
# in each chain for each iteration
#fixed effects
plot(nl_fit, 
     nvariables = 10,
     variable = "^b_", 
     regex = TRUE)

plot(
  conditional_effects(x = nl_fit, 
                      spaghetti = TRUE, 
                      ndraws = 2e2,
                      effects = 'log10_intensity')
)

#and effects of colour (mean correct accuracy just slightly above baseline)
plot(
  conditional_effects(x = nl_fit, 
                      effects = 'colour')
)

#summary of parameter estimates
full_sm = summary(nl_fit,
                  robust = TRUE)#use the median estimate

# . . Check parameter estimates -----------------------------------------------

#We can also check if these estimates match with our expectations.
#Extract fixed effects estimates
full_fix = full_sm$fixed
#extract rownames
full_fix_rn = rownames(full_fix)

# stimulus levels to add to the dataset to be able to measure this effect.
full_estimates = full_fix$Estimate
names(full_estimates) = full_fix_rn

# to predict all effects, we can use the 'posterior_epred' method
#by default 100 predictions per continuous variable (but fewer for their interactions)
system.time(
  {
    full_cond =brms::conditional_effects(nl_fit, 
                                         method = 'posterior_epred', # posterior epred not working
                                         cores =  parallel::detectCores()-1,
                                         effects = c('log10_intensity:colour')
    )
  }
)#takes <2 seconds
#extract predictions
pred_data = full_cond$`log10_intensity:colour`
#plot each stimulus type
with(ed,
     {
       plot(x = log10_intensity,
            y = accuracy,
            main = 
    'Edrich et al. 1979 Fig. 2. Dance precision
    as a function of light intensity & wavelength',
            xlab = 'log10(intensity) (photons/cm2/s)',
            ylab = 'accuracy (r)',
            ylim = c(0,1),
            xlim = c(9,14),
            pch = 21,
            cex = 1.5,
            col = 'black',
            bg = sapply(paste(wavelength),
                        FUN = switch,
                        `354` = 'purple',
                        `429` = 'blue',
                        `535` = 'darkgreen',
                        'red')
       )
       abline(h = c(0,1))
     }
)
with(data.frame(t(full_estimates)),
     {
  abline(h = c(0,1,
               plogis(LogitBase_Intercept),
               1- plogis(LogitLapse_Intercept) ), 
         lty = c(1,1,3,3)
         )
     }
)

#plot total prediction intervals
with(subset(pred_data, colour == '0green'), 
     {
       polygon(x = c(sort(log10_intensity), rev(sort(log10_intensity))), 
               y = c(lower__[order(log10_intensity)],
                     rev(upper__[order(log10_intensity)])
               ), 
               col = adjustcolor('green', alpha.f = 25/256),
               border = NA,
               lwd = 0.1
       )
     }
)
with(subset(pred_data, colour == 'blue'), 
     {
       polygon(x = c(sort(log10_intensity), rev(sort(log10_intensity))), 
               y = c(lower__[order(log10_intensity)],
                     rev(upper__[order(log10_intensity)])
               ), 
               col = adjustcolor('cyan', alpha.f = 25/256),
               border = NA,
               lwd = 0.1
       )
     }
)
with(subset(pred_data, colour == 'UV'), 
     {
       polygon(x = c(sort(log10_intensity), rev(sort(log10_intensity))), 
               y = c(lower__[order(log10_intensity)],
                     rev(upper__[order(log10_intensity)])
               ), 
               col = adjustcolor('magenta', alpha.f = 25/256),
               border = NA,
               lwd = 0.1
       )
     }
)

#plot the median prediction lines
with(subset(pred_data, colour == '0green'), 
     lines(x = sort(log10_intensity),
           y = estimate__[order(log10_intensity)], 
           col = 'darkgreen',
           lwd = 3)
)
with(subset(pred_data, colour == 'blue'), 
     lines(x = sort(log10_intensity),
           y = estimate__[order(log10_intensity)], 
           col = 'darkblue',
           lwd = 3)
)
with(subset(pred_data, colour == 'UV'), 
     lines(x = sort(log10_intensity),
           y = estimate__[order(log10_intensity)], 
           col = 'purple4',
           lwd = 3)
)

legend(x = 'bottomright',
       legend = paste(sort(unique(ed$wavelength),
                     decreasing = TRUE), 'nm'),
       col = c('darkgreen','blue', 'purple'),
       pch = c(20, 20, 20),
       lwd = 3)
# Beta Psychometric version ----------------------------------------------------


#set up model fit
formula_beta = bf(
  #set up a formula for the curve as a whole,
  #including parameters found in the data (correct_incorrect, stimulus)
  #and parameters that we wish to estimate (baseline, lapse rate, inflection point, width).
  #Most of these are subject to further fixed (type) and random (animal) effects,
  #these need to be defined for each parameter.
  
  #Two parameters need special transformations
  #To keep the output between 0 and 1, additional effects of lapse rate
  #will be added on the "logit" scale (Lapse = inv_logit(LogitLapse)).
  #To avoid curve widths of 0, we can assume a positive slope (≥0)
  #and add additional effects to width on a log scale (Width = exp(LogWidth))
  formula = accuracy ~ 
    inv_logit(LogitBase) + (1 - inv_logit(LogitLapse) - inv_logit(LogitBase) ) *#curve region
    inv_logit( 4.39*(log10_intensity - Inflex) / exp(LogWidth) ) , #inflection-width curve
  # for each of these parameters, we can set up a separate formula 
  # that describes how to predict them from the data 
  #Base rate of correct choices: "
  LogitBase ~ 1, #Base rate of correct choices: "~ 1" gives the instruction "estimate the mean across all data"
  #Lapse rate on a log(odds) scale:
  LogitLapse ~ colour, #this is similar to the formula in our LMM example
  #inflection point of the initial curve:
  Inflex ~ colour, #N.B. this is similar to the intercept, so it does not include effects of stimulus level
  #log 80% width of the curve:
  LogWidth ~ colour, #N.B. this is similar to the slope, so all of its effects depend on stimulus level
  family = Beta(link = "identity", link_phi = 'log'),
  nl = TRUE)#the joint distribution for these parameters is undefined, and therefore the parameters themselves are "nonlinear"


#set up priors
prior_beta = get_prior(formula = formula_beta,
                     data = ed)

print(
  with(prior_beta,
       paste0(class, '_', nlpar, '_', coef))
  )

# . . Base rate prior -----------------------------------------------------
#for the baseline, we will use a logit-normal distribution with a bias towards 0.5

#N.B. this seems to require slightly tighter priors

#set the prior distribution
prior_beta = within(prior_beta, 
                  { prior[
                    class %in% 'b' & #just the fixed effects
                      nlpar %in% 'LogitBase' &
                      coef %in% 'Intercept' 
                  ] = 'normal(-1,1)' #a normal distribution centred on plogis(-1) = 0.27
                  })
#this prior is automatically bounded between 0 and 1, 
# . . Lapse rate priors ---------------------------------------------------

# set the prior distribution
prior_beta = within(prior_beta, 
                  { prior[
                    class %in% 'b' & #just the fixed effects
                      nlpar %in% 'LogitLapse' &
                      coef %in% 'Intercept' 
                  ] = 'normal(-3,2)' #a normal distribution centred on -3
                  })

#for all other fixed effects on lapse rate, we'll suggest values around 0 (no effect)
prior_beta = within(prior_beta, 
                  { prior[
                    class %in% 'b' & #just the fixed effects
                      nlpar %in% 'LogitLapse' &
                      coef %in% 'colourblue' 
                  ] = 'normal(0,2)' #a normal distribution:mean 0, sd 3
                  })
prior_beta = within(prior_beta, 
                  { prior[
                    class %in% 'b' & #just the fixed effects
                      nlpar %in% 'LogitLapse' &
                      coef %in% 'colourUV' 
                  ] = 'normal(0,2)' #a normal distribution:mean 0, sd 3
                  })

# . . Inflection point priors ---------------------------------------------

#set the priors
prior_beta = within(prior_beta, 
                  { prior[
                    class %in% 'b' & #just the fixed effects
                      nlpar %in% 'Inflex' &
                      coef %in% 'Intercept'
                  ] = 'normal(11,2)' #a normal distribution:mean 11, sd 2
                  })
#for all coefficients, we'll suggest values around 0 (no effect)
prior_beta = within(prior_beta, 
                  { prior[
                    class %in% 'b' & #just the fixed effects
                      nlpar %in% 'Inflex' &
                      coef %in% 'colourblue' 
                  ] = 'normal(0,2)' #a normal distribution:mean 0, sd 2
                  })

prior_beta = within(prior_beta, 
                  { prior[
                    class %in% 'b' & #just the fixed effects
                      nlpar %in% 'Inflex' &
                      coef %in% 'colourUV' 
                  ] = 'normal(0,2)' #a normal distribution:mean 0, sd 2
                  })

# . . Rise region width priors --------------------------------------------

#set prior distribution
prior_beta = within(prior_beta, 
                  { prior[
                    class %in% 'b' & #just the fixed effects
                      nlpar %in% 'LogWidth' &
                      coef %in% 'Intercept' 
                  ] = 'normal(1,2)' #a normal distribution:mean 1, sd 3
                  })
#for all coefficients, we'll suggest values around 0 (no effect)
prior_beta = within(prior_beta, 
                  { prior[
                    class %in% 'b' & #just the fixed effects
                      nlpar %in% 'LogWidth' &
                      coef %in% 'colourblue'
                  ] = 'normal(0,2)' #a normal distribution:mean 0, sd 3
                  })
prior_beta = within(prior_beta, 
                  { prior[
                    class %in% 'b' & #just the fixed effects
                      nlpar %in% 'LogWidth' &
                      coef %in% 'colourUV'
                  ] = 'normal(0,2)' #a normal distribution:mean 0, sd 3
                  })
#inspect assigned priors
print(prior_beta)
##                 prior class       coef group resp dpar      nlpar lb ub  source
##  student_t(3, 0, 2.5) sigma                                        0    default
##                (flat)     b                                Inflex       default
##           normal(0,3)     b colourblue                     Inflex       default
##           normal(0,3)     b   colourUV                     Inflex       default
##           normal(3,3)     b  Intercept                     Inflex       default
##                (flat)     b                             LogitBase       default
##           normal(0,1)     b  Intercept                  LogitBase       default
##           (flat)     b                            LogitLapse       default
##           normal(0,3)     b colourblue                 LogitLapse       default
##           normal(0,3)     b   colourUV                 LogitLapse       default
##           normal(-3,3)     b  Intercept                 LogitLapse       default
##                (flat)     b                              LogWidth       default
##           normal(0,3)     b colourblue                   LogWidth       default
##           normal(0,3)     b   colourUV                   LogWidth       default
##           normal(1,3)     b  Intercept                   LogWidth       default

#double check that the prior distribution is viable by first setting up a short dummy run
# Dummy run
system.time(
  {
    dummy_beta = brm( formula = formula_beta, # using our nonlinear formula
                     data = ed, # our data
                     prior = prior_beta, # our priors 
                     sample_prior = 'only', #ignore the data to check the influence of the priors
                     iter = 300, # short run for 300 iterations
                     chains = 4, # 4 chains in parallel
                     cores = 4, # on 4 CPUs
                     refresh = 0, # don't echo chain progress
                     control = list(adapt_delta = 0.9), #finer sampling
                     backend = 'cmdstanr') # use cmdstanr (other compilers broken)
  }
)
# On my computer this takes <60s, each chain running for <1 seconds (mainly compile time)


#the default plot shows the values estimated for each parameter
# in each chain for each iteration
#fixed effects
plot(dummy_beta, 
     nvariables = 10,
     variable = "^b_", 
     regex = TRUE)


#double check that the prior distribution is viable by first setting up a short dummy run
# Dummy run
system.time(
  {
    beta_fit = brm( formula = formula_beta, # using our nonlinear formula
                     data = ed, # our data
                     prior = prior_beta, # our priors 
                     iter = 4000, # long run for 2000 iterations
                     chains = 4, # 4 chains in parallel
                     cores = 4, # on 4 CPUs
                     refresh = 0, # don't echo chain progress
                     control = list(adapt_delta = 0.95), #finer sampling
                     backend = 'cmdstanr') # use cmdstanr (other compilers broken)
  }
)
# On my computer this takes <60s, each chain running for <1 seconds (mainly compile time)


#the default plot shows the values estimated for each parameter
# in each chain for each iteration
#fixed effects
plot(beta_fit, 
     nvariables = 10,
     variable = "^b_", 
     regex = TRUE)

plot(
  conditional_effects(x = beta_fit, 
                      spaghetti = TRUE, 
                      ndraws = 2e2,
                      effects = 'log10_intensity')
)

#and effects of colour (mean correct accuracy just slightly above baseline)
plot(
  conditional_effects(x = beta_fit, 
                      effects = 'colour')
)

#summary of parameter estimates
full_sm_beta = summary(beta_fit,
                  robust = TRUE)#use the median estimate

# . . Check parameter estimates -----------------------------------------------

#We can also check if these estimates match with our expectations.
#Extract fixed effects estimates
full_fix_beta = full_sm_beta$fixed
#extract rownames
full_fix_rn_beta = rownames(full_fix_beta)

# stimulus levels to add to the dataset to be able to measure this effect.
full_estimates_beta = full_fix_beta$Estimate
names(full_estimates_beta) = full_fix_rn_beta

# to predict all effects, we can use the 'posterior_epred' method
#by default 100 predictions per continuous variable (but fewer for their interactions)
system.time(
  {
    full_cond_beta =brms::conditional_effects(beta_fit, 
                                         method = 'posterior_epred', # posterior epred not working
                                         cores =  parallel::detectCores()-1,
                                         effects = c('log10_intensity:colour')
    )
  }
)#takes <2 seconds
#extract predictions
pred_data_beta = full_cond_beta$`log10_intensity:colour`
#plot each stimulus type
with(ed,
     {
       plot(x = log10_intensity,
            y = accuracy,
            main = 
    'Edrich et al. 1979 Fig. 2. Dance precision
    as a function of light intensity & wavelength',
            xlab = 'log10(intensity) (photons/cm2/s)',
            ylab = 'accuracy (r)',
            ylim = c(0,1),
            xlim = c(9,14),
            pch = 21,
            cex = 1.5,
            col = 'black',
            bg = sapply(paste(wavelength),
                        FUN = switch,
                        `354` = 'purple',
                        `429` = 'blue',
                        `535` = 'darkgreen',
                        'red')
       )
       abline(h = c(0,1))
     }
)
with(data.frame(t(full_estimates_beta)),
     {
  abline(h = c(0,1,
               plogis(LogitBase_Intercept),
               1- plogis(LogitLapse_Intercept) ), 
         lty = c(1,1,3,3)
         )
     }
)

#plot total prediction intervals
with(subset(pred_data_beta, colour == '0green'), 
     {
       polygon(x = c(sort(log10_intensity), rev(sort(log10_intensity))), 
               y = c(lower__[order(log10_intensity)],
                     rev(upper__[order(log10_intensity)])
               ), 
               col = adjustcolor('green', alpha.f = 25/256),
               border = NA,
               lwd = 0.1
       )
     }
)
with(subset(pred_data_beta, colour == 'blue'), 
     {
       polygon(x = c(sort(log10_intensity), rev(sort(log10_intensity))), 
               y = c(lower__[order(log10_intensity)],
                     rev(upper__[order(log10_intensity)])
               ), 
               col = adjustcolor('cyan', alpha.f = 25/256),
               border = NA,
               lwd = 0.1
       )
     }
)
with(subset(pred_data_beta, colour == 'UV'), 
     {
       polygon(x = c(sort(log10_intensity), rev(sort(log10_intensity))), 
               y = c(lower__[order(log10_intensity)],
                     rev(upper__[order(log10_intensity)])
               ), 
               col = adjustcolor('magenta', alpha.f = 25/256),
               border = NA,
               lwd = 0.1
       )
     }
)

#plot the median prediction lines
with(subset(pred_data_beta, colour == '0green'), 
     lines(x = sort(log10_intensity),
           y = estimate__[order(log10_intensity)], 
           col = 'darkgreen',
           lwd = 3)
)
with(subset(pred_data_beta, colour == 'blue'), 
     lines(x = sort(log10_intensity),
           y = estimate__[order(log10_intensity)], 
           col = 'darkblue',
           lwd = 3)
)
with(subset(pred_data_beta, colour == 'UV'), 
     lines(x = sort(log10_intensity),
           y = estimate__[order(log10_intensity)], 
           col = 'purple4',
           lwd = 3)
)

legend(x = 'bottomright',
       legend = paste(sort(unique(ed$wavelength),
                     decreasing = TRUE), 'nm'),
       col = c('darkgreen','blue', 'purple'),
       pch = c(20, 20, 20),
       lwd = 3)


# Compare with model without different maxima -----------------------------


#set up model fit
formula_beta_equal = bf(
  #set up a formula for the curve as a whole,
  #including parameters found in the data (correct_incorrect, stimulus)
  #and parameters that we wish to estimate (baseline, lapse rate, inflection point, width).
  #Most of these are subject to further fixed (type) and random (animal) effects,
  #these need to be defined for each parameter.
  
  #Two parameters need special transformations
  #To keep the output between 0 and 1, additional effects of lapse rate
  #will be added on the "logit" scale (Lapse = inv_logit(LogitLapse)).
  #To avoid curve widths of 0, we can assume a positive slope (≥0)
  #and add additional effects to width on a log scale (Width = exp(LogWidth))
  formula = accuracy ~ 
    inv_logit(LogitBase) + (1 - inv_logit(LogitLapse) - inv_logit(LogitBase) ) *#curve region
    inv_logit( 4.39*(log10_intensity - Inflex) / exp(LogWidth) ) , #inflection-width curve
  # for each of these parameters, we can set up a separate formula 
  # that describes how to predict them from the data 
  #Base rate of correct choices: "
  LogitBase ~ 1, #Base rate of correct choices: "~ 1" gives the instruction "estimate the mean across all data"
  #Lapse rate on a log(odds) scale:
  LogitLapse ~ 1, #this is similar to the formula in our LMM example
  #inflection point of the initial curve:
  Inflex ~ colour, #N.B. this is similar to the intercept, so it does not include effects of stimulus level
  #log 80% width of the curve:
  LogWidth ~ colour, #N.B. this is similar to the slope, so all of its effects depend on stimulus level
  family = Beta(link = "identity", link_phi = 'log'),
  nl = TRUE)#the joint distribution for these parameters is undefined, and therefore the parameters themselves are "nonlinear"


#set up priors
prior_beta_equal = get_prior(formula = formula_beta_equal,
                       data = ed)

print(
  with(prior_beta_equal,
       paste0(class, '_', nlpar, '_', coef))
)

# . . Base rate prior -----------------------------------------------------
#for the baseline, we will use a logit-normal distribution with a bias towards 0.5

#N.B. this seems to require even tighter priors

#set the prior distribution
prior_beta_equal = within(prior_beta_equal, 
                    { prior[
                      class %in% 'b' & #just the fixed effects
                        nlpar %in% 'LogitBase' &
                        coef %in% 'Intercept' 
                    ] = 'normal(-1,1)' #a normal distribution centred on plogis(-1) = 0.27
                    })
#this prior is automatically bounded between 0 and 1, 
# . . Lapse rate priors ---------------------------------------------------

# set the prior distribution
prior_beta_equal = within(prior_beta_equal, 
                    { prior[
                      class %in% 'b' & #just the fixed effects
                        nlpar %in% 'LogitLapse' &
                        coef %in% 'Intercept' 
                    ] = 'normal(-3.5,1.5)' #a normal distribution centred on -3
                    })

# . . Inflection point priors ---------------------------------------------

#set the priors
prior_beta_equal = within(prior_beta_equal, 
                    { prior[
                      class %in% 'b' & #just the fixed effects
                        nlpar %in% 'Inflex' &
                        coef %in% 'Intercept'
                    ] = 'normal(11,2)' #a normal distribution:mean 11, sd 2
                    })
#for all coefficients, we'll suggest values around 0 (no effect)
prior_beta_equal = within(prior_beta_equal, 
                    { prior[
                      class %in% 'b' & #just the fixed effects
                        nlpar %in% 'Inflex' &
                        coef %in% 'colourblue' 
                    ] = 'normal(0,2)' #a normal distribution:mean 0, sd 2
                    })

prior_beta_equal = within(prior_beta_equal, 
                    { prior[
                      class %in% 'b' & #just the fixed effects
                        nlpar %in% 'Inflex' &
                        coef %in% 'colourUV' 
                    ] = 'normal(0,2)' #a normal distribution:mean 0, sd 2
                    })

# . . Rise region width priors --------------------------------------------

#set prior distribution
prior_beta_equal = within(prior_beta_equal, 
                    { prior[
                      class %in% 'b' & #just the fixed effects
                        nlpar %in% 'LogWidth' &
                        coef %in% 'Intercept' 
                    ] = 'normal(1,1.5)' #a normal distribution:mean 1, sd 3
                    })
#for all coefficients, we'll suggest values around 0 (no effect)
prior_beta_equal = within(prior_beta_equal, 
                    { prior[
                      class %in% 'b' & #just the fixed effects
                        nlpar %in% 'LogWidth' &
                        coef %in% 'colourblue'
                    ] = 'normal(0,1.5)' #a normal distribution:mean 0, sd 3
                    })
prior_beta_equal = within(prior_beta_equal, 
                    { prior[
                      class %in% 'b' & #just the fixed effects
                        nlpar %in% 'LogWidth' &
                        coef %in% 'colourUV'
                    ] = 'normal(0,1.5)' #a normal distribution:mean 0, sd 3
                    })
#inspect assigned priors
print(prior_beta_equal)
##                 prior class       coef group resp dpar      nlpar lb ub  source
##  student_t(3, 0, 2.5) sigma                                        0    default
##                (flat)     b                                Inflex       default
##           normal(0,3)     b colourblue                     Inflex       default
##           normal(0,3)     b   colourUV                     Inflex       default
##           normal(3,3)     b  Intercept                     Inflex       default
##                (flat)     b                             LogitBase       default
##           normal(0,1)     b  Intercept                  LogitBase       default
##           (flat)     b                            LogitLapse       default
##           normal(0,3)     b colourblue                 LogitLapse       default
##           normal(0,3)     b   colourUV                 LogitLapse       default
##           normal(-3,3)     b  Intercept                 LogitLapse       default
##                (flat)     b                              LogWidth       default
##           normal(0,3)     b colourblue                   LogWidth       default
##           normal(0,3)     b   colourUV                   LogWidth       default
##           normal(1,3)     b  Intercept                   LogWidth       default

#double check that the prior distribution is viable by first setting up a short dummy run
# Dummy run
system.time(
  {
    dummy_beta_equal = brm( formula = formula_beta_equal, # using our nonlinear formula
                      data = ed, # our data
                      prior = prior_beta_equal, # our priors 
                      sample_prior = 'only', #ignore the data to check the influence of the priors
                      iter = 300, # short run for 300 iterations
                      chains = 4, # 4 chains in parallel
                      cores = 4, # on 4 CPUs
                      refresh = 0, # don't echo chain progress
                      control = list(adapt_delta = 0.95), #finer sampling
                      backend = 'cmdstanr') # use cmdstanr (other compilers broken)
  }
)
# On my computer this takes <60s, each chain running for <1 seconds (mainly compile time)


#the default plot shows the values estimated for each parameter
# in each chain for each iteration
#fixed effects
plot(dummy_beta_equal, 
     nvariables = 10,
     variable = "^b_", 
     regex = TRUE)


#double check that the prior distribution is viable by first setting up a short dummy run
# Dummy run
system.time(
  {
    beta_equal_fit = brm( formula = formula_beta_equal, # using our nonlinear formula
                    data = ed, # our data
                    prior = prior_beta_equal, # our priors 
                    iter = 4000, # long run for 2000 iterations
                    chains = 4, # 4 chains in parallel
                    cores = 4, # on 4 CPUs
                    refresh = 0, # don't echo chain progress
                    control = list(adapt_delta = 0.95), #finer sampling
                    backend = 'cmdstanr') # use cmdstanr (other compilers broken)
  }
)
# On my computer this takes <60s, each chain running for <1 seconds (mainly compile time)


#the default plot shows the values estimated for each parameter
# in each chain for each iteration
#fixed effects
plot(beta_equal_fit, 
     nvariables = 10,
     variable = "^b_", 
     regex = TRUE)

plot(
  conditional_effects(x = beta_equal_fit, 
                      spaghetti = TRUE, 
                      ndraws = 2e2,
                      effects = 'log10_intensity')
)

#and effects of colour (mean correct accuracy just slightly above baseline)
plot(
  conditional_effects(x = beta_equal_fit, 
                      effects = 'colour')
)

#summary of parameter estimates
full_sm_beta_equal = summary(beta_equal_fit,
                       robust = TRUE)#use the median estimate

# . . Check parameter estimates -----------------------------------------------

#We can also check if these estimates match with our expectations.
#Extract fixed effects estimates
full_fix_beta_equal = full_sm_beta_equal$fixed
#extract rownames
full_fix_rn_beta_equal = rownames(full_fix_beta_equal)

# stimulus levels to add to the dataset to be able to measure this effect.
full_estimates_beta_equal = full_fix_beta_equal$Estimate
names(full_estimates_beta_equal) = full_fix_rn_beta_equal

# to predict all effects, we can use the 'posterior_epred' method
#by default 100 predictions per continuous variable (but fewer for their interactions)
system.time(
  {
    full_cond_beta_equal =brms::conditional_effects(beta_equal_fit, 
                                              method = 'posterior_epred', # posterior epred not working
                                              cores =  parallel::detectCores()-1,
                                              effects = c('log10_intensity:colour')
    )
  }
)#takes <2 seconds
#extract predictions
pred_data_beta_equal = full_cond_beta_equal$`log10_intensity:colour`
#plot each stimulus type
with(ed,
     {
       plot(x = log10_intensity,
            y = accuracy,
            main = 
'Edrich et al. 1979 Fig. 2. Dance precision
as a function of light intensity & wavelength
—equal maxima',
            xlab = 'log10(intensity) (photons/cm2/s)',
            ylab = 'accuracy (r)',
            ylim = c(0,1),
            xlim = c(9,14),
            pch = 21,
            cex = 1.5,
            col = 'black',
            bg = sapply(paste(wavelength),
                        FUN = switch,
                        `354` = 'purple',
                        `429` = 'blue',
                        `535` = 'darkgreen',
                        'red')
       )
       abline(h = c(0,1))
     }
)
with(data.frame(t(full_estimates_beta_equal)),
     {
       abline(h = c(0,1,
                    plogis(LogitBase_Intercept),
                    1- plogis(LogitLapse_Intercept) ), 
              lty = c(1,1,3,3)
       )
     }
)

#plot total prediction intervals
with(subset(pred_data_beta_equal, colour == '0green'), 
     {
       polygon(x = c(sort(log10_intensity), rev(sort(log10_intensity))), 
               y = c(lower__[order(log10_intensity)],
                     rev(upper__[order(log10_intensity)])
               ), 
               col = adjustcolor('green', alpha.f = 25/256),
               border = NA,
               lwd = 0.1
       )
     }
)
with(subset(pred_data_beta_equal, colour == 'blue'), 
     {
       polygon(x = c(sort(log10_intensity), rev(sort(log10_intensity))), 
               y = c(lower__[order(log10_intensity)],
                     rev(upper__[order(log10_intensity)])
               ), 
               col = adjustcolor('cyan', alpha.f = 25/256),
               border = NA,
               lwd = 0.1
       )
     }
)
with(subset(pred_data_beta_equal, colour == 'UV'), 
     {
       polygon(x = c(sort(log10_intensity), rev(sort(log10_intensity))), 
               y = c(lower__[order(log10_intensity)],
                     rev(upper__[order(log10_intensity)])
               ), 
               col = adjustcolor('magenta', alpha.f = 25/256),
               border = NA,
               lwd = 0.1
       )
     }
)

#plot the median prediction lines
with(subset(pred_data_beta_equal, colour == '0green'), 
     lines(x = sort(log10_intensity),
           y = estimate__[order(log10_intensity)], 
           col = 'darkgreen',
           lwd = 3)
)
with(subset(pred_data_beta_equal, colour == 'blue'), 
     lines(x = sort(log10_intensity),
           y = estimate__[order(log10_intensity)], 
           col = 'darkblue',
           lwd = 3)
)
with(subset(pred_data_beta_equal, colour == 'UV'), 
     lines(x = sort(log10_intensity),
           y = estimate__[order(log10_intensity)], 
           col = 'purple4',
           lwd = 3)
)

legend(x = 'bottomright',
       legend = paste(sort(unique(ed$wavelength),
                           decreasing = TRUE), 'nm'),
       col = c('darkgreen','blue', 'purple'),
       pch = c(20, 20, 20),
       lwd = 3)


# Model comparison --------------------------------------------------------
#How robust is the model to changes in the data structure?
#Would the model make good predictions refitted the model but dropped one datapoint, would it still give good predictions?
# calculate the Leave-One-Out (LOO) cross validation metric for the model
loo_normal = loo(nl_fit)#calculate for full model 
loo_beta = loo(beta_fit)#calculate for beta distributed model
loo_equal = loo(beta_equal_fit)#calculate for equal-lapse model

loo_compare(loo_normal,
            loo_beta,
            loo_equal)


# Inspect curve equations for use elsewhere -------------------------------

#Green light curve
g_curve = 
with(data.frame(t(full_estimates_beta)),
     {
       list(Base = plogis(LogitBase_Intercept),
            Maxima = 1- plogis(LogitLapse_Intercept), 
            Inflection = Inflex_Intercept,  
            Width = exp(LogWidth_Intercept)  
       )
     }
)

#UV light curve
u_curve = 
with(data.frame(t(full_estimates_beta)),
     {
       list(Base = plogis(LogitBase_Intercept),
            Maxima = 1- plogis(LogitLapse_Intercept + LogitLapse_colourUV), 
            Inflection = Inflex_Intercept + Inflex_colourUV,  
            Width = exp(LogWidth_Intercept + LogWidth_colourUV)  
       )
     }
)

print(cbind(g_curve, u_curve),digits = 3)
##           g_curve u_curve
## Base       0.142   0.142  
## Maxima     0.951   0.816  
## Inflection 11.2    10.4   
## Width      1.23    0.643  

# Comparing UV to green
#accuracy reaches 50% at nearly 1 log unit lower intensity for UV
print(full_fix_beta['Inflex_colourUV',c(1,3:4)], digits = 3)
##                 Estimate   l-95% CI   u-95% CI
## Inflex_colourUV   -0.783   -0.965   -0.567

#But maximum accuracy is lower by 1.5 log units for UV
print(full_fix_beta['LogitLapse_colourUV',c(1,3:4)], digits = 3)
##                         Estimate l-95% CI u-95% CI
## LogitLapse_colourUV     1.49    0.508     2.37

#while there is a noticeable difference in width,
#this overlaps with 0 on a log scale
print(full_fix_beta['LogWidth_colourUV',c(1,3:4)], digits = 3)
##                   Estimate l-95% CI u-95% CI
## LogWidth_colourUV   -0.645    -2.34    0.321

log10(4e12)
log10(4e14)

## colour brightn Min.    1st Qu. Median  Mean    3rd Qu. Max.
## g       h      0.2527  0.8145  0.9058  0.8487  0.9691  1.0000
## u       h      0.9107  0.9286  0.9571  0.9520  0.9648  0.9873
## g       l      0.1324  0.5281  0.7348  0.6531  0.8581  0.9847
## u       l      0.0767  0.4977  0.6664  0.6461  0.8146  0.9154

gh_meanQ = list(mn = 0.2527,
                  q1 = 0.8145,
             q2 = 0.9058,
             q3 = 0.9691,
             mx = 1.0000)

gl_meanQ = list(mn = 0.1324,
                q1 = 0.5281,
             q2 = 0.7348,
             q3 = 0.8581,
             mx = 0.9847)

uh_meanQ = list(mn = 0.9107,
                  q1 = 0.9286,
             q2 = 0.9571,
             q3 = 0.9648,
             mx = 0.9873)

ul_meanQ = list(mn = 0.0767,
                q1 = 0.4977,
             q2 = 0.6664,
             q3 = 0.8146,
             mx = 0.9154)

xxl = seq(from = 9, to = 16,
          length.out = 1e3)

#plot each stimulus type
ed_ug = subset(ed, wavelength != 429)

with(ed_ug,
     {
       plot(x = log10_intensity,
            y = accuracy,
            main = 
'Comparison of data from Edrich et al. 1979
with Kühn, 2017',
            xlab = 'log10(intensity) (photons/cm2/s)',
            ylab = 'accuracy (r)',
            ylim = c(0,1),
            xlim = c(9,16),
            pch = 21,
            cex = 1.5,
            col = 'black',
            bg = sapply(paste(wavelength),
                        FUN = switch,
                        `354` = 'purple',
                        `429` = 'blue',
                        `535` = 'darkgreen',
                        'red')
       )
       abline(h = c(0,1))
       abline(v = log10(4*c(1e12,1e14)) )
     }
)

#plot the median prediction lines
with(g_curve, 
     lines(x = xxl,
           y = Base + (Maxima  - Base) * 
             plogis(4.39*(xxl - Inflection)/Width), 
           col = 'darkgreen',
           lwd = 3,
           lty = 3)
)
with(u_curve, 
     lines(x = xxl,
           y = Base + (Maxima  - Base) * 
             plogis(4.39*(xxl - Inflection)/Width), 
           col = 'purple4',
           lwd = 3,
           lty = 3)
)

#plot observed data from Kühn
with(gh_meanQ,
     {
       polygon(x = 0.1*c(-1,1,1,-1)+log10(4e14)-0.1,
               y = c(q1, q1, q3, q3),
               col = adjustcolor('green', alpha.f = 25/256),
               border = 'green'
               )
       lines(x =  0.1*c(-1,1)+log10(4e14)-0.1,
             y = c(q2, q2),
             lwd = 5,
             lend = 'butt',
             col = 'darkgreen')
     }
     )

with(gl_meanQ,
     {
       polygon(x = 0.1*c(-1,1,1,-1)+log10(4e12)-0.1,
               y = c(q1, q1, q3, q3),
               col = adjustcolor('green', alpha.f = 25/256),
               border = 'green'
               )
       lines(x =  0.1*c(-1,1)+log10(4e12)-0.1,
             y = c(q2, q2),
             lwd = 5,
             lend = 'butt',
             col = 'darkgreen')
     }
     )
with(uh_meanQ,
     {
       polygon(x = 0.1*c(-1,1,1,-1)+log10(4e14)+0.1,
               y = c(q1, q1, q3, q3),
               col = adjustcolor('magenta', alpha.f = 25/256),
               border = 'magenta'
               )
       lines(x =  0.1*c(-1,1)+log10(4e14)+0.1,
             y = c(q2, q2),
             lwd = 5,
             lend = 'butt',
             col = 'purple')
     }
     )

with(ul_meanQ,
     {
       polygon(x = 0.1*c(-1,1,1,-1)+log10(4e12)+0.1,
               y = c(q1, q1, q3, q3),
               col = adjustcolor('magenta', alpha.f = 25/256),
               border = 'magenta'
               )
       lines(x =  0.1*c(-1,1)+log10(4e12)+0.1,
             y = c(q2, q2),
             lwd = 5,
             lend = 'butt',
             col = 'purple')
     }
     )


legend(x = 'bottomright',
       legend = c(
                   paste('Edrich', 
                        sort(unique(ed_ug$wavelength),
                           decreasing = TRUE), 'nm'),
                   paste('  Kühn', c(528, 365), 'nm' ) 
       ),
       col = c('darkgreen', 'purple', 'green', 'magenta'),
       pch = c(20, 20, 15, 15),
       lty = c(3,3,NA, NA),
       lwd = 3)


# Compare with visual pigment template ------------------------------------

# Make a spline template for a visual pigment
StavengaSpline = function(spec_range = c(300, 700), #bounds of spectrum in nanometers
                           lambda_max,#peak sensitivity
                           a_type = 'a1'){#pigment type, only a1 available currently
  wlns =  seq(from = min(spec_range),
              to = max(spec_range), 
              length.out = 1e3) #
  #Stavenga, D. G. (2010). On visual pigment templates and the spectral shape of invertebrate rhodopsins and metarhodopsins. Journal of Comparative Physiology A: Neuroethology, Sensory, Neural, and Behavioral Physiology, 196(11), 869–878. doi:10.1007/s00359-010-0568-7
  # modified lognormal
  Mlognorm  = function(wl,l_max, a0, a1)
  {
    x = log10(wl/l_max)
    return(
      exp(-a0*x^2 *
            (1+a1*x+3*a1^2*x^2)
          )
      )
  }
  if(a_type == 'a1')
  {
    #alpha band
    a_band = Mlognorm(wl = wlns,
                      l_max = lambda_max, 
                      a0 = 380, 
                      a1 = 6.09)
    #beta band
    b_band = 0.29*Mlognorm(wl = wlns, 
                           l_max = 340, 
                           a0 = 247, 
                           a1 = 3.59)
    #gamma band
    g_band = 1.99*Mlognorm(wl = wlns, 
                           l_max = 276, 
                           a0 = 647, 
                           a1 = 23.4)
  }else
  {stop('a2 and a3 pigments not yet implemented')}
  # N.B. Stavenga normalises to max(a.band), I normalise to function max
  r_stav  = (a_band + b_band + g_band)/
              max(a_band + b_band + g_band)
  return(    smooth.spline(x = wlns,
                           y =  r_stav)    )
}#StavengaSpline <- function(spec.range, lambda.max)


g_photopig = predict(StavengaSpline(lambda_max = 526))
u_photopig = predict(StavengaSpline(lambda_max = 344))

g_Kühn_g = predict(StavengaSpline(lambda_max = 526), x = 528)
g_Kühn_u = predict(StavengaSpline(lambda_max = 526), x = 365)
u_Kühn_g = predict(StavengaSpline(lambda_max = 344), x = 528)
u_Kühn_u = predict(StavengaSpline(lambda_max = 344), x = 365)
g_Edrich_g = predict(StavengaSpline(lambda_max = 526), x = 535)
g_Edrich_u = predict(StavengaSpline(lambda_max = 526), x = 354)
u_Edrich_g = predict(StavengaSpline(lambda_max = 344), x = 535)
u_Edrich_u = predict(StavengaSpline(lambda_max = 344), x = 354)

ww = seq(from = 300,
         to = 700,
         length.out = 1e3)
fwhm_sigma = 2 * sqrt(2 * log(2))

#The UV LED has a full-width at half max of ~20nm
rough_UV_LED = dnorm(x = ww,
                     mean = 365,
                     sd = 20/fwhm_sigma
                     )/dnorm(0, sd = 20/fwhm_sigma)
#The green LED has a full-width at half max of ~40nm
rough_G_LED = dnorm(x = ww,
                     mean = 528,
                     sd = 40/fwhm_sigma
                    )/dnorm(0, sd = 40/fwhm_sigma)
#The UV filter has a full-width at half max of ~80nm
rough_UV_filter = dnorm(x = ww,
                     mean = 354,
                     sd = 80/fwhm_sigma
                     )/dnorm(0, sd = 80/fwhm_sigma)
#The green filter has a full-width at half max of ~90nm
rough_G_filter = dnorm(x = ww,
                     mean = 528,
                     sd = 90/fwhm_sigma
                    )/dnorm(0, sd = 90/fwhm_sigma)

ug_ratio = g_photopig$y/
            u_photopig$y
#somehow this produces negative values?
ug_ratio[ug_ratio <0] = 0
#and absurd outliers
ug_ratio[ug_ratio > mean(ug_ratio) + 6*sd(ug_ratio)] =
  max(ug_ratio[ug_ratio < mean(ug_ratio) + 6*sd(ug_ratio)])
ug_ratio[ww > 550] = predict(smooth.spline(x = ww,
                                           y = ug_ratio,
                                           spar = 0.5)
                             )$y[ww>550]

plot(NULL,
     xlim = c(300, 700),
     ylim = c(0,1),
     xlab = 'wavelength (nm)',
     ylab = 'relative absoption | emission | transmission',
     main = 'comparison of light sources\nused by Edrich & Kühn')
abline(h = c(0,1))
lines(u_photopig,
      lwd = 3,
      col = 'purple')
lines(g_photopig,
      lwd = 3,
      col = 'darkgreen')
polygon(x = c(ww, rev(range(ww))),
        y = c(rough_UV_LED, 0, 0),
        col = adjustcolor('magenta', alpha.f = 25/255),
        border = 'magenta')
polygon(x = c(ww, rev(range(ww))),
        y = c(rough_G_LED, 0, 0),
        col = adjustcolor('green', alpha.f = 25/255),
        border = 'green')
polygon(x = c(ww, rev(range(ww))),
        y = c(rough_UV_filter, 0, 0),
        col = adjustcolor('magenta4', alpha.f = 25/255),
        border = 'magenta4',
        lty = 3)
polygon(x = c(ww, rev(range(ww))),
        y = c(rough_G_filter, 0, 0),
        col = adjustcolor('green3', alpha.f = 25/255),
        border = 'green3',
        lty = 3)
points(x = c(354, 365, 528, 535),
       y = c(0,0,0,0),
       pch = c(4,3,3,4),
       lwd = 3,
       col = c('magenta4', 'magenta', 'green2', 'green4'))
legend(x = 'topright',
       legend = c('UV opsin',
                  'green opsin',
                  'Schott UG1 filter',
                  'Schott VG9 filter',
                  'UV LED',
                  'green LED'),
       col = c('purple',
               'darkgreen',
               'magenta4',
               'green4',
               'magenta',
               'green2'),
       pch = c(NA, NA, 4, 4, 3, 3),
       lty = c(1, 1, 3, 3, 1, 1),
       lwd = 2,
       cex = 0.5
       )

plot(NULL,
     xlim = c(300, 700),
     ylim = c(0.1,1e23),
     xlab = 'wavelength (nm)',
     ylab = 'G/UV relative sensitivity',
     main = 'Comparison of light sources\nwith predicted relative sensitivity',
     log = 'y')
abline(h = 1,
       lwd = 3,
       lty = 2)
lines(x = ww,
      y = ug_ratio,
      lwd = 5,
      col = 'black')
abline(v = c(354, 365, 528, 535))
par(new = TRUE)
plot(NULL,
     xlim = c(300, 700),
     ylim = c(0,1),
     xlab = '',
     ylab = '',
      axes = F)
axis(side = 4, at = 0:4/4)
polygon(x = c(ww, rev(range(ww))),
        y = c(rough_UV_LED, 0, 0),
        col = adjustcolor('magenta', alpha.f = 25/255),
        border = 'magenta')
polygon(x = c(ww, rev(range(ww))),
        y = c(rough_G_LED, 0, 0),
        col = adjustcolor('green', alpha.f = 25/255),
        border = 'green')
polygon(x = c(ww, rev(range(ww))),
        y = c(rough_UV_filter, 0, 0),
        col = adjustcolor('magenta4', alpha.f = 25/255),
        border = 'magenta4',
        lty = 3)
polygon(x = c(ww, rev(range(ww))),
        y = c(rough_G_filter, 0, 0),
        col = adjustcolor('green3', alpha.f = 25/255),
        border = 'green3',
        lty = 3)
