#FOR A 'CLEAN' RUN, PRESS ctrl+shift+F10 to RESTART Rstudio
graphics.off()
# Details ---------------------------------------------------------------
#       AUTHOR:	James Foster              DATE: 2025 05 30
#     MODIFIED:	James Foster              DATE: 2025 05 30
#
#  DESCRIPTION: Load data and plot
#               
#       INPUTS: 
#               
#      OUTPUTS: csv
#
#	   CHANGES: - 
#
#   REFERENCES: Edrich, W., Neumeyer, C. and von Helversen, O. (1979).
#               “Anti-sun orientation” of bees with regard to a field of ultraviolet light. 
#               J. Comp. Physiol. 134, 151–157.
#
#
#    EXAMPLES:  
#
# 
#TODO   ---------------------------------------------
#TODO   
#- Load data  
#- Plot
#- Fit curves




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
       pch = c(20, 20, 20))
