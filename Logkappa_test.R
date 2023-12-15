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
