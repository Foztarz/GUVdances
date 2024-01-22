require(circular)
raw_data <- data.frame(bearings = round(c(
  sapply(X = c(3,3,3,3,3,5,5,5,5,5),#concentration parameter "kappa"
         FUN = function(k, ...)
         {
           suppressWarnings(
             rvonmises(n = 10, 
                       mu = runif(n=1,min = -pi, max = pi),
                       kappa = k,
                       ...)
           )
         },
         control.circular = list(template = 'geographics',
                                 units = 'degrees',
                                 mod = '2pi')
  )
)),
condition = sort(rep(c('A','B'), 50)), #two conditions with 50 angles
beetle = sort(rep(1:5, 10)) #5 beetles with 10 angles per condition
)
View(raw_data)#inspect data
#calculate mean vectors
mean_vectors = aggregate(bearings~beetle*condition,
                         data = raw_data,
                         FUN = function(i, ...){
                           suppressWarnings(
                             {rho.circular(x = as.circular(i, ...))}
                           )
                         },
                         control.circular = list(template = 'geographics',
                                                 units = 'degrees',
                                                 mod = '2pi')
)
#correct names
mean_vectors <- within(mean_vectors,
                       {mean_vector = bearings; rm(bearings)} # bearing now indicates a mean vector, not an angle
)
#plot mean vectors
boxplot(mean_vector~condition, data = mean_vectors, ylim = c(0,1))
#calculate log kappa
mean_vectors <- within(mean_vectors,
                       {
                         kappa <- circular::A1inv(mean_vector)
                         logkappa <- log(kappa)
                       }
)
#plot log kappa
boxplot(logkappa~condition, data = mean_vectors)# variances differ
#check for normal distribution
shap_test <- aggregate(logkappa~condition, data = mean_vectors, FUN = shapiro.test)
#they are at least normally distributed (p>0.05)
within(shap_test, {p <- logkappa; rm(logkappa)})


# Try a linear model with log kappa ---------------------------------------
require(lme4)
lm_test = lmer(formula = logkappa~condition+(1|beetle), #cannot account for different variances
               data = mean_vectors)
lm_null = lmer(formula = logkappa~1+(1|beetle), 
               data = mean_vectors)
anova(lm_test, lm_null) # model with condition is better, but not significantly


# Fit beta distributed model ----------------------------------------------
require(glmmTMB)

glm_test = glmmTMB::glmmTMB(formula = mean_vector~condition+(1|beetle),
                            dispformula = ~ condition, # cannot account for indiv effects on dispersion
                            family = beta_family(link = "logit"),
                            data = mean_vectors 
                            # control = glmmTMBControl(optimizer = nloptr::bobyqa)
                            )

glm_null = glmmTMB::glmmTMB(formula = mean_vector~1+(1|beetle),
                            dispformula = ~ 1,
                            family = beta_family(link = "logit"),
                            data = mean_vectors#, 
                            )
#compare all linear models
all_aic = AIC(lm_test, lm_null, glm_test, glm_null)
print(all_aic[order(all_aic[,2],decreasing = TRUE),])#glm is objectively better
# df       AIC
# lm_test   4  19.94452
# lm_null   3  19.57929
# glm_null  3 -24.43961
# glm_test  5 -34.39136
anova(glm_test, glm_null) # significant effect of condition on a beta scale
summary(glm_test) #no significant effect of condition on dispersion


# Fit a Bayesian Beta distributed model -----------------------------------
require(brms)

frm_bmod = bf(mean_vector~condition+(1|beetle),
              phi ~ condition+(1|beetle),
              family = brms::Beta()) 
priors_bmod = get_prior(formula = frm_bmod, 
                        data = mean_vectors)
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

system.time(expr = 
  {
bmod_test = brm(formula = frm_bmod ,
                 data = mean_vectors,
                iter = 5e3, # about 3 minutes
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

plot(bmod_test)
summary(bmod_test)
cond_eff = conditional_effects(bmod_test)
plot(cond_eff)

bmod_loo = loo(bmod_test)
bmod_looic = bmod_loo$estimates[3,]

all_ic = rbind(all_aic, bmod_se = rev(bmod_looic))
print(all_ic[order(all_ic[,2],decreasing = TRUE),])#bmod is better, but not by >1SE
