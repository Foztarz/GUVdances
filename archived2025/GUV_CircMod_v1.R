#FOR A 'CLEAN' RUN, PRESS ctrl+shift+F10 to RESTART Rstudio
graphics.off()
# Details ---------------------------------------------------------------
#       AUTHOR:	James Foster              DATE: 2024 12 05
#     MODIFIED:	James Foster              DATE: 2024 12 20
#
#  DESCRIPTION: Attempt to run a two-way interaction model on the GUV dances data
#               using the circular modulo modelling method devel. by Jake Graving.
#               
#       INPUTS: 
#               
#      OUTPUTS: Plots and test statistics
#
#	   CHANGES: - Investigation of MLE fit to each individual
#             - CircMLE plots
#
#   REFERENCES: Sayin, S., Graving, J., et al. in revision
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
#- Read in data +
#- Inspect data +
#- Set up priors for simple model +
#- Generate Stan code for simple model  +
#- Set up Stan parameters for simple model  +
#- Set up parameters for random effect model +
#- Random effects interactions model (really!) +
#- Try larger subset +
#- Extract predictions +
#- Correct zmu predictions
#- Plot predictions 1/2
#- Try longer run for full model
#- Simplify model
# - Model comparison
# - Vectorise zkappa


# Set up workspace --------------------------------------------------------


## Load packages ----------------------------------------------------------
#needs installing before first use (in Rstudio, see automatic message)
suppressMessages(#these are disturbing users unnecessarily
  {
    # require(CircStats)#package for circular hypothesis tests #Is this being used at all?
    require(circular)#package for handling circular data
    require(CircMLE)#package for circular mixture models
    require(brms)#package for preparing Stan models
  }
)


## Plot spacing function -------------------------------------------------
#generates the same spacing as R's default barplot function
BarSpacer = function(n, 
                     spa = 0.2,
                     wdt = 1.0)
{
  seq(from = spa+1-wdt/2,
      to = n*(1+spa)-wdt/2,
      length.out = n )
}

## General functions -----------------------------------------------------


#convert angles to signed angles in (-180, 180)
Mod360.180 = function(x)
{#use atan2 to convert any angle to the range (-180,180)
  deg(
    atan2(y = sin(rad(x)),
          x = cos(rad(x))
    )
  )
}

#the radian equivalent
mod_circular = function(x)
{
  atan2(y = sin(x),
        x = cos(x))
}

#the unwrapped (no discontinuities)
unwrap_circular = function(x)
{
  mux = mean.circular(x = circular(x = x, template = 'none'))
  centx = atan2(y = sin(x - mux),
                x = cos(x  - mux))
  unwrx = centx + mux
}

#degree version
unwrap_circular_deg = function(x)
{
  mux = mean.circular(x = circular(x = x, template = 'none'))
  centx = atan2(y = sin(x - mux),
                x = cos(x  - mux))
  unwrx = centx + mux
  return(deg(unwrx))
}

#invert the softplus link
#https://en.wikipedia.org/wiki/Softplus
#we are using this as our _inverse_ link function for kappa,
#maps almost 1:1 but keeps values >0 for low estimates
softplus = function(x)
{
  log(exp(x)+1) 
}
#this would return our kappa estimates back to the original scale
inv_softplus = function(x)
{
  log(exp(x)-1) 
}

#convert inv_softplus scaled kappa to mean vector estimate
Softpl_to_meanvec = function(x)
{
  circular::A1(
    softplus(x)
  )
}

#convert circular to normalised
NormCirc = function(x,
                    plusmean = TRUE)
{
  mn = mean.circular(x) * as.numeric(plusmean)
  return(mod_circular(x - mn) + mn)
}
## von Mises model inspection functions ----------------------------------

#histograms on a vertical axis
#any data, but plotted as a histogram on a vertical rather tahn horizontal axis
VertHist = function(data, # numerical data vector
                    breaks = 1e2,
                    ylab = 'data',
                    xlab = 'density',
                    ylim = NULL,
                    main = '',
                    col = 'gray',
                    border = NA,
                    ...)
{
  hst = hist(x = data, # calculate the histogram but don't plot it
             breaks = breaks, # user defined breaks
             plot = FALSE)
  with(hst,
       {
         plot(x = NULL, #open an empty plot
              xlim = c(0, max(density)),
              ylim = if(is.null(ylim)){range(mids)}else{ylim},
              xlab = xlab,
              ylab = ylab,
              main = main)
         #plot each bar
         for(i in 1:length(mids))
         {
           rect(xleft = 0,
                xright = density[i],
                ybottom = breaks[i], 
                ytop = breaks[i + 1],
                col = col,
                border = border,
                ...
           )
         }
       }
  )
}

#function to construct the transformation list
#transforms variable estimates in the model output
#by default circular estimates are converted to the 360° interval around the mean
#and kappa estimates are converted to softplus and then to mean vector estimates
Make_vmtrans = function(mod,
                        angfun = NormCirc,
                        kapfun = Softpl_to_meanvec)
{
  fx_nm = names(fixef(mod)[,1])
  nm_ang = fx_nm[grepl(pattern = 'mu', x = fx_nm)]
  nm_kap = fx_nm[grepl(pattern = 'kappa', x = fx_nm)]
  tlst = c(replicate(n = length(nm_ang),
                     expr = angfun,
                     simplify = FALSE),
           replicate(n = length(nm_kap),
                     expr = kapfun,
                     simplify = FALSE) )
  names(tlst) = paste0('b_', c(nm_ang, nm_kap))
  return(tlst)
}

#plot the variable estimates from a von Mises model
PlotVMfix = function(mod,
                     angfun = NormCirc,
                     kapfun = Softpl_to_meanvec)
{
  plot(mod,
       nvariables = 10,
       variable = "^b_", 
       transformations = Make_vmtrans(mod, angfun, kapfun),
       regex = TRUE)
}

## Stan variables ---------------------------------------------------------


### Functions ------------------------------------------------------------


#set up the modulo function
mod_circular_fun = stanvar(scode = "
  real mod_circular(real y) {
    return atan2(sin(y), cos(y));
  }
",
                           block = 'functions')
#set up a von Mises PDF that converts to modulo (adapted from BRMS default)
von_mises3_fun = stanvar(scode = "
real von_mises3_lpdf(real y, real mu, real kappa) {
     if (kappa < 100) {
       return von_mises_lpdf(mod_circular(y) | mu, kappa);
     } else {
       return normal_lpdf(mod_circular(y) | mu, sqrt(1 / kappa));
     }
   }
",
                         block = 'functions')
#a circular mean that might be useful
meancirc_fun = stanvar(scode = "
   // calculate the mean angle of a circular distribution (in radians)
   real mean_circular(vector y){
      int N = size(y);
      real sumsin = 0;
      real sumcos = 0;
     for (n in 1:N){
       sumsin += sin(y[n]);
       sumcos += cos(y[n]);
     }
     sumsin = sumsin/N;
     sumcos = sumcos/N;
     return(atan2(  sumsin, sumcos) );
   }
",
                       block = 'functions')
#could unwrap estimate; not especially useful as only one estimate per iteration
unwrapcirc_fun = stanvar(scode = "
   // calculate the unwrapped version (no discontinuities)
vector unwrap_circular(vector y)
{
  int N = size(y);
  real mux = mean_circular(y);
  vector[N] centx = y - mux;
  for(n in 1:N)
  {
  centx[n] = mod_circular(centx[n]);
  }
  vector[N] unwrx = centx + mux;
  return(unwrx);
}
",
                         block = 'functions')

stan_mvm_fun = meancirc_fun + 
  mod_circular_fun + 
  unwrapcirc_fun +
  von_mises3_fun

### Parameters -----------------------------------------------------------


#generate modulo outputs of fixed effects
# mu_gen = stanvar(scode = "
# real mu_circ = mod_circular(b_fmu[1]);
# real mu_offs = mod_circular(b_fmu[2]);
# ",
#                  block = 'genquant')
# #or maybe?
mu_gen = stanvar(scode = "
vector [M_1] mu_circ; //modulo circular estimate
for (i in 1:size(b_fmu)){
mu_circ[i] = mod_circular(b_fmu[i]);
}
",
                 block = 'genquant')

#random effects on mean angle
  # zmu_var = stanvar(scode = "
  # vector[K_zmu] zmu_id;  // regression coefficients;
  # for (i in 1:K_zmu){
  #   zmu_id[i] = mod_circular(b_zmu[i]); //each in modulus format
  # }
  # ", 
  # block = 'genquant')

#concentration of random effects on mean angle
#superceded by  slopes?
  # zkappa_var = stanvar(scode = "
  # real zkappa;", 
  #                      block = "parameters") + #define in the parameters block
  #   stanvar(scode = "
  # real kappa_id = log1p_exp(zkappa);
  #           ", 
  #           block = 'genquant') #inverse softplus (estimated on softplus scale)

# vector[K_zmu] zmu_id;  // regression coefficients;
zmu_var_slope = stanvar(scode = "
vector[N_1] zmu_id_condition1;  // condition on coefficients per indiv;
vector[N_1] zmu_id_condition2;  // condition on coefficients per indiv;
vector[N_1] zmu_id_condition3;  // condition on coefficients per indiv;
vector[N_1] zmu_id_condition4;  // condition on coefficients per indiv;
for (i in 1:N_1){
zmu_id_condition1[i] = mod_circular(b_zmu[i]); //each in modulus format
zmu_id_condition2[i] = mod_circular(b_zmu[i+N_1]); //each in modulus format
zmu_id_condition3[i] = mod_circular(b_zmu[i+N_1*2]); //each in modulus format
zmu_id_condition4[i] = mod_circular(b_zmu[i+N_1*3]); //each in modulus format
}
          ", 
block = 'genquant')
# zmu_id[i] = mod_circular(b_zmu[i]); //each in modulus format

#concentration of random effects on fixed effect on mean angle
zkappa_var_slope = stanvar(scode = "
real zkappa1;
real zkappa2;
real zkappa3;
real zkappa4;
                           ",
                           block = "parameters") + 
  stanvar(scode = "
real kappa_id_condition1 = log1p_exp(zkappa1);
real kappa_id_condition2 = log1p_exp(zkappa1+zkappa2);
real kappa_id_condition3 = log1p_exp(zkappa1+zkappa3);
real kappa_id_condition4 = log1p_exp(zkappa1+zkappa4);
          ", 
          block = 'genquant')



stanvars_intercepts = stan_mvm_fun + mu_gen  #+ zmu_var+ zkappa_var
stanvars_slopes = stan_mvm_fun + mu_gen  + zkappa_var_slope + 
  zmu_var_slope #includes intercepts+ zkappa_var

# Input Variables ----------------------------------------------------------

all_plots = FALSE # to speed up
#  .  User input -----------------------------------------------------------



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
    caption = 'Please select the "colour_dance_reorg.csv" file'
  )
}else{
  message('\n\nPlease select the "colour_dance_reorg.csv" file\n\n')
  Sys.sleep(0.5)#goes too fast for the user to see the message on some computers
  path_file <- file.choose(new=F)
}
#show the user the path they have selected
if(is.null(path_file) | !length(path_file))
{stop('No file selected.')}else
{print(path_file)}


# Read in the data and format ---------------------------------------------

#select the reorganised data
cd = read.table(file = path_file, 
                header = T, 
                sep  = ',')
View(cd)

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
#add kappa
mle_estimates = aggregate(angle~ID*brightn*colour,
                         data = cd,
                         FUN = function(x)
                         {with(mle.vonmises(circular(x,
                                                    template = 'none'),
                                            bias = TRUE), c(mu, kappa))}
)
#add to the summary table and 
#calculate inverse softplus kappa
mean_vectors = within(mean_vectors,
                      {
                        mu = deg(mle_estimates$angle[,1])
                        kappa = mle_estimates$angle[,2]
                        iskappa  = inv_softplus(kappa)
                      }
                      )
#plot mean vectors
boxplot(mean_vector~colour*brightn,
        data = mean_vectors,
        ylim = c(0,1),
        col = adjustcolor(c('green2',
                            'purple',
                            'darkgreen',
                            'purple4'), alpha.f = 0.5),
        outline = FALSE)
stripchart(mean_vector~colour*brightn,
           data = mean_vectors,
           add = TRUE,
           vertical = TRUE,
           method = 'jitter',
           col = adjustcolor(c('green2',
                               'purple',
                               'darkgreen',
                               'purple4'), alpha.f = 0.5),
           bg = gray(level = 0,
                     alpha = 0.1),
           pch = 21)

stim = apply(X = expand.grid(c = c('g','u'),b = c('h','l')),
             FUN = paste,
             collapse = '',
             MARGIN = 1)
invisible(
  {
  lapply(X = u_id,
         FUN = function(id)
         {
            
           with(subset(x = mean_vectors,
                      subset = ID %in% id),
                {
           lines(x = 1:length(stim),
                     y = mean_vector[ match(x = stim,
                                            table =  paste0(colour,brightn),
                                            nomatch = NA) ],
                 col = gray(0,0.2)
                 )
                })
         })
})
abline(h = c(0,1))


#plot kappa estimates
boxplot(kappa~colour*brightn,
        data = mean_vectors,
        ylim = c(0.01,500),
        log = 'y',
        col = adjustcolor(c('green2',
                            'purple',
                            'darkgreen',
                            'purple4'), alpha.f = 0.5),
        outline = FALSE)
stripchart(kappa~colour*brightn,
           data = mean_vectors,
           add = TRUE,
           vertical = TRUE,
           method = 'jitter',
           col = adjustcolor(c('green2',
                               'purple',
                               'darkgreen',
                               'purple4'), alpha.f = 0.5),
           bg = gray(level = 0,
                     alpha = 0.1),
           pch = 21)
invisible(
  {
    lapply(X = u_id,
           FUN = function(id)
           {
             
             with(subset(x = mean_vectors,
                         subset = ID %in% id),
                  {
                    lines(x = 1:length(stim),
                          y = kappa[ match(x = stim,
                                                 table =  paste0(colour,brightn),
                                                 nomatch = NA) ],
                          col = gray(0,0.2)
                    )
                  })
           })
  })



#plot inverse softplus kappa (model scale)
boxplot(iskappa~colour*brightn,
        data = mean_vectors,
        ylim = c(-10,100),
        col = adjustcolor(c('green2',
                            'purple',
                            'darkgreen',
                            'purple4'), alpha.f = 0.5),
        outline = FALSE)
stripchart(iskappa~colour*brightn,
           data = mean_vectors,
           add = TRUE,
           vertical = TRUE,
           method = 'jitter',
           col = adjustcolor(c('green2',
                               'purple',
                               'darkgreen',
                               'purple4'), alpha.f = 0.5),
           bg = gray(level = 0,
                     alpha = 0.1),
           pch = 21)

invisible(
  {
    lapply(X = u_id,
           FUN = function(id)
           {
             
             with(subset(x = mean_vectors,
                         subset = ID %in% id),
                  {
                    lines(x = 1:length(stim),
                          y = iskappa[ match(x = stim,
                                                 table =  paste0(colour,brightn),
                                                 nomatch = NA) ],
                          col = gray(0,0.2)
                    )
                  })
           })
  })
abline(h = inv_softplus(A1inv(c(0.2,0.7))), # benchmarks could be -0.7 and 2.0
       lty = 3)
mtext(side = 4,
     at = inv_softplus(A1inv(c(0.2,0.7))),
     text = paste(' rho = ', c(0.2,0.7)),
     las = 2,
     cex = 0.4
     )

#plot angles
#by brightness
boxplot(mu~colour*brightn,
        data = mean_vectors,
        ylim = c(-180,180),
        col = NA,
        border = NA,
        outline = FALSE)
stripchart(mu~colour*brightn,
           data = mean_vectors,
           add = TRUE,
           vertical = TRUE,
           method = 'jitter',
           col = adjustcolor(c('green2',
                               'purple',
                               'darkgreen',
                               'purple4'), alpha.f = 0.5),
           bg = gray(level = 0,
                     alpha = 0.1),
           pch = 21)
invisible(
  {
    lapply(X = u_id,
           FUN = function(id)
           {
             
             with(subset(x = mean_vectors,
                         subset = ID %in% id),
                  {
                    lines(x = 1:length(stim),
                          y = mu[ match(x = stim,
                                             table =  paste0(colour,brightn),
                                             nomatch = NA) ],
                          col = gray(0,0.2)
                    )
                  })
           })
  })
abline(h = c(-180, -90, 0, 90,  180),
       lty = c(1,3,1,3,1))
#plot angles
#by colour
boxplot(mu~brightn*colour,
        data = mean_vectors,
        ylim = c(-180,180),
        col = NA,
        border = NA,
        outline = FALSE)
stripchart(mu~brightn*colour,
           data = mean_vectors,
           add = TRUE,
           vertical = TRUE,
           method = 'jitter',
           col = adjustcolor(c('green2',
                               'darkgreen',
                               'purple',
                               'purple4'), alpha.f = 0.5),
           bg = gray(level = 0,
                     alpha = 0.1),
           pch = 21)
invisible(
  {
    lapply(X = u_id,
           FUN = function(id)
           {
             
             with(subset(x = mean_vectors,
                         subset = ID %in% id),
                  {
                    lines(x = 1:length(stim),
                          y = mu[ match(x = sort(stim),#now alphabetical (g before u, h before l)
                                             table =  paste0(colour,brightn),
                                             nomatch = NA) ],
                          col = gray(0,0.2)
                    )
                  })
           })
  })
abline(h = c(-180, -90, 0, 90,  180),
       lty = c(1,3,1,3,1))
#plot angles
#by colour
boxplot(mu~brightn,
        data = subset(mean_vectors, colour == 'u'),
        ylim = c(-180,180)*2,
        col = NA,
        border = NA,
        outline = FALSE)
invisible(
  {
    lapply(X = u_id,
           FUN = function(id)
           {
             
             with(subset(x = mean_vectors,
                         subset = ID %in% id),
                  {
                    points(x = 1:2 + rnorm(1,sd=0.05),
                          y = unwrap_circular_deg(
                                  mu[ match(x = sort(stim)[3:4],#just UV now alphabetical (g before u, h before l)
                                             table =  paste0(colour,brightn),
                                             nomatch = NA) ]
                                  ),
                          type = 'b',
                          col = adjustcolor(gray(0,0.5),offset = c(runif(3),0))
                    )
                  })
           })
  })
abline(h = -4:4*90,
       lty = c(1,3))

#average orientation changes
#uv
mu_diff_u = sapply(X = u_id,
                   FUN = function(id)
                   {
                     
                     with(subset(x = mean_vectors,
                                 subset = ID %in% id & colour %in% 'u'),
                          {
                            deg(
                            mod_circular(
                           rad(mu[brightn %in% 'l'] -  mu[brightn %in% 'h'])
                            )
                            )
                          })
                   })
mu_diff_u = unlist(mu_diff_u)

#green
mu_diff_g = sapply(X = u_id,
                   FUN = function(id)
                   {
                     
                     with(subset(x = mean_vectors,
                                 subset = ID %in% id & colour %in% 'g'),
                          {
                            deg(
                            mod_circular(
                           rad(mu[brightn %in% 'l'] -  mu[brightn %in% 'h'])
                            )
                            )
                          })
                   })
mu_diff_g = unlist(mu_diff_g)

#high green to uv
mu_diff_hgu = sapply(X = u_id,
                   FUN = function(id)
                   {
                     
                     with(subset(x = mean_vectors,
                                 subset = ID %in% id ),
                          {
                            deg(
                            mod_circular(
                           rad(mu[brightn %in% 'h'& colour %in% 'u'] -  
                                 mu[brightn %in% 'h' & colour %in% 'g'])
                            )
                            )
                          })
                   })
mu_diff_hgu = unlist(mu_diff_hgu)

#high green to uv
mu_diff_lgu = sapply(X = u_id,
                     FUN = function(id)
                     {
                       
                       with(subset(x = mean_vectors,
                                   subset = ID %in% id ),
                            {
                              deg(
                                mod_circular(
                                  rad(mu[brightn %in% 'l'& colour %in% 'u'] -  
                                        mu[brightn %in% 'l' & colour %in% 'g'])
                                )
                              )
                            })
                     })
mu_diff_lgu = unlist(mu_diff_lgu)

par(mfrow = c(2,2))
VertHist(unwrap_circular_deg(mu_diff_g),
         breaks = 24,
         ylim = c(-360,360),
         ylab = 'Change in mean angle (unwrapped)',
         col = 'darkgreen',
         main = 'High to low green')
abline(h = -4:4*90,
       lty = c(1,3))

VertHist(unwrap_circular_deg(mu_diff_u),
         ylim = c(-360,360),
         breaks = 24,
         ylab = 'Change in mean angle (unwrapped)',
         col = 'purple3',
         main =  'High to low UV')
abline(h = -4:4*90,
       lty = c(1,3))

VertHist(unwrap_circular_deg(mu_diff_hgu),
         breaks = 24,
         ylim = c(-360,360),
         ylab = 'Change in mean angle (unwrapped)',
         col = 'gray',
         main = 'High green to UV')
abline(h = -4:4*90,
       lty = c(1,3))

VertHist(unwrap_circular_deg(mu_diff_lgu),
         ylim = c(-360,360),
         breaks = 24,
         ylab = 'Change in mean angle (unwrapped)',
         col = 'black',
         main =  'Low green to UV')
abline(h = -4:4*90,
       lty = c(1,3))

cmle_g = circ_mle(data = circular(mu_diff_g, template = 'geographics'))
cmle_u = circ_mle(data = circular(mu_diff_u, template = 'geographics'))
cmle_h = circ_mle(data = circular(mu_diff_hgu, template = 'geographics'))
cmle_l = circ_mle(data = circular(mu_diff_lgu, template = 'geographics'))

par(mfrow = c(2,2), mar = c(0,0,0,0))
par(pty = 's')
plot_circMLE(data = circular(x = unwrap_circular_deg(mu_diff_g),
                             units = 'degrees',
                             rotation = 'clock',
                             modulo = '2pi',
                             zero = pi/2), 
             table = cmle_g,
             col = c('darkgreen',"red", "black", "black"))
mtext(text = 'High to low green')                  

plot_circMLE(data = circular(x = unwrap_circular_deg(mu_diff_u),
                             units = 'degrees',
                             rotation = 'clock',
                             modulo = '2pi',
                             zero = pi/2), 
             table = cmle_u,
             col = c('purple3',"red", "black", "black"))
mtext(text = 'High to low UV')                
plot_circMLE(data = circular(x = unwrap_circular_deg(mu_diff_hgu),
                             units = 'degrees',
                             rotation = 'clock',
                             modulo = '2pi',
                             zero = pi/2), 
             table = cmle_h,
             col = c('lightgray',"red", "black", "black"))
mtext(text = 'High green to UV')                  

plot_circMLE(data = circular(x = unwrap_circular_deg(mu_diff_lgu),
                             units = 'degrees',
                             rotation = 'clock',
                             modulo = '2pi',
                             zero = pi/2), 
             table = cmle_l,
             col = c('darkgray',"red", "black", "black"))
mtext(text =  'Low green to UV')              


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
                                    size =  20,
                                    replace = FALSE)
                )


## Formula ---------------------------------------------------------------

#set up model fit
formula_int_slope = bf(
  formula = angle ~ #set up a formula for the mean angle, modulus to (-pi,pi)
    mod_circular(fmu + zmu), # mean angle combines fixed and random effects
  fmu ~  BR + CL + BR:CL, # fixed effects only change as a function of condition
  zmu ~  0+ID + ID:BR + ID:CL + ID:BR:CL, #random effects change as a function of individual and condition and their interaction
  kappa ~ BR + CL + BR:CL + (1 + BR + CL + BR:CL|ID), #for kappa this occurs in linear space, and BRMS can set it up automatically
  family = von_mises(link = "identity", # the mean angle will be returned as-is
                     link_kappa = 'softplus'),#kappa will be returned via the softplus link https://en.wikipedia.org/wiki/Softplus
  nl = TRUE)#to accept user-defined extra parameters (zmu) we need to treat the formula as nonlinear


# sc = make_stancode(formula = formula_int_slope,
#                    data = cd)
# write.table(x = sc,
#             file = file.path(dirname(path_file),
#                              'sc_CircMod_v1.stan')
#           )
## Priors ----------------------------------------------------------------
prior_int_slope = get_prior(formula = formula_int_slope,
                            data = cd_subs,
                            check = FALSE)
dim(prior_int_slope)
#suggests 695  possible priors!
#they will take this general structure
print(subset(prior_int_slope, nlpar %in% 'fmu')['coef'])
## coef
## ''    
## BRl 
## BRl:CLu
## CLu 
## Intercept 


### assign BRMS default priors -------------------------------------------

#ideally the ones for zkappa should be vectorised like the ones for kappa
# for (n in 1:N) {
#   // add more terms to the linear predictor
#   kappa[n] += r_1_kappa_1[J_1[n]] * Z_1_kappa_1[n] + r_1_kappa_2[J_1[n]] * Z_1_kappa_2[n] + r_1_kappa_3[J_1[n]] * Z_1_kappa_3[n] + r_1_kappa_4[J_1[n]] * Z_1_kappa_4[n];
# }

prior_int_slope = within(prior_int_slope,
                         {
   #fixed effects on mean angle are von Mises distributed (von_mises3 converts estimates to modulo (-pi,pi))                              
   prior[nlpar %in% 'fmu' & coef %in% 'Intercept'] = 'von_mises3(0, 0.1)'# very weak bias to zero (could be no bias?)
   prior[nlpar %in% 'fmu' & class %in% 'b'] = 'von_mises3(0, 1.0)'#moderate bias to zero, no effect
   #random effects on mean angle are von Mises distributed, with a kappa parameter estimated from the data
   #the intercept condition is high intensity green light
   prior[nlpar %in% 'zmu' & coef %in% 'Intercept'] = 'von_mises3(0, log1p_exp(zkappa1))'
   prior[nlpar %in% 'zmu' & class %in% 'b'] = 'von_mises3(0, log1p_exp(zkappa1))'
   prior[nlpar %in% 'zmu' & class %in% 'b' 
         & grepl(pattern = 'BRl', #the random effect of low brightness has a different kappa
                 x = coef)] = 'von_mises3(0, log1p_exp(zkappa1+zkappa2))'
   prior[nlpar %in% 'zmu' & class %in% 'b' 
         & grepl(pattern = 'CLu', #the random effect of UV has a different kappa
                 x = coef)] = 'von_mises3(0, log1p_exp(zkappa1+zkappa3))'
   prior[nlpar %in% 'zmu' & class %in% 'b' 
         & grepl(pattern = 'BRl:CLu', #the random effect of low brightness & UV has a different kappa
                 x = coef)] = 'von_mises3(0, log1p_exp(zkappa1+zkappa4))'
   #fixed effects on kappa are normally distributed on the softplus scale
   prior[dpar %in% 'kappa' & class %in% 'Intercept'] = 'normal(2.0, 5.0)'#weak expectation of kappa around 2 (mean vector around 0.70)
   prior[dpar %in% 'kappa' & class %in% 'b'] = 'normal(0.0, 5.0)'#weak expectation of condition effect around 0
   #random effects on kappa are t-distributed on the softplus scale
   prior[dpar %in% 'kappa' & class %in% 'sd'] = 'student_t(3, 0, 5.0)' #wide prior 
                         }
)


### add extra priors for the random effects mean angles -----------------



# prior_int_slope = prior_int_slope + #random effects kappas are t-distributed on a softplus scale
#   set_prior("target += student_t_lpdf(zkappa1 | 3, 25, 5)", #expect high concentration (low variation) 
#             check = FALSE)+
#   set_prior("target += student_t_lpdf(zkappa1+zkappa2 | 3, 25, 5)", #expect high concentration (low variation) 
#             check = FALSE)+ #random effects kappas are t-distributed on a softplus scale
#   set_prior("target += student_t_lpdf(zkappa1+zkappa3 | 3, 25, 5)", #expect high concentration (low variation) 
#             check = FALSE)+
#   set_prior("target += student_t_lpdf(zkappa1+zkappa4 | 3, 25, 5)", #expect high concentration (low variation) 
#             check = FALSE)

#different strategy, anchor estimates around zkappa 1 but with wide dist
#now attempting very wide dist on zkappa intercept, to include probability mass around zero
prior_int_slope = prior_int_slope + #random effects kappas are t-distributed on a softplus scale
  set_prior("target += student_t_lpdf(zkappa1 | 3, 30, 20)", #expect high concentration (low variation) 
            check = FALSE)+
  set_prior("target += student_t_lpdf(zkappa2 | 3, 0, 5)", #expect high concentration (low variation) 
            check = FALSE)+ #random effects kappas are t-distributed on a softplus scale
  set_prior("target += student_t_lpdf(zkappa3 | 3, 0, 5)", #expect high concentration (low variation) 
            check = FALSE)+
  set_prior("target += student_t_lpdf(zkappa4 | 3, 0, 5)", #expect high concentration (low variation) 
            check = FALSE)



sc = make_stancode(formula = formula_int_slope,
                   data = cd,
                   prior = prior_int_slope,
                   stanvars = stanvars_slopes)
write.table(x = sc,
            file = file.path(dirname(path_file),
                             'sc_CircMod_v1.stan'),
            quote = FALSE,
            col.names = FALSE,
            row.names = FALSE)


# ## Dummy run to check the influence of the priors ------------------
# 
# 
# #double check that the prior distribution is viable by first setting up a short dummy run
# # Dummy run
# #Warning takes 15 min just to compile!
# #TODO work out why this samples less efficiently than with data
# system.time( #currently takes about 60 minutes for 10000 iterations
#   {
#     dummy_int_slope = brm( formula = formula_int_slope, # using our nonlinear formula
#                            data = cd_subs, # our data
#                            prior = prior_int_slope, # our priors 
#                            stanvars = stanvars_slopes,
#                            sample_prior = 'only', #ignore the data to check the influence of the priors
#                            iter = 10000, # can only estimate with enough iterations for params
#                            chains = 4, # 4 chains in parallel
#                            cores = 4, # on 4 CPUs
#                            refresh = 0, # don't echo chain progress
#                            backend = 'cmdstanr') # use cmdstanr (other compilers broken)
#   }
# )
# 
# if(all_plots)
# {
#   plot(dummy_int_slope,
#        variable = 'fmu',
#        regex = TRUE,
#        transform = unwrap_circular_deg)
#   plot(dummy_int_slope,
#        variable = '^kappa_id',
#        regex = TRUE)
#   #samples inefficiently?
#   plot(dummy_int_slope,
#        variable = '^zkappa',
#        regex = TRUE)
#   # plot(dummy_int_slope,
#   #      variable = 'zmu_id',
#   #      transform = unwrap_circular_deg,
#   #      nvariables = 5,
#   #      ask = FALSE)
#   plot(dummy_int_slope,
#        variable = '^zmu_id_condition',
#        transform = unwrap_circular_deg,
#        nvariables = 5,
#        regex = TRUE,
#        ask = FALSE)
#   # plot(dummy_int_slope)
# }


## Subset run --------------------------------------------------------------

# subset run
system.time(#takes less than 35 minutes for 50 individuals
  {
    full_int_slope = brm( formula = formula_int_slope, # using our nonlinear formula
                          data = cd_subs, # our data
                          prior = prior_int_slope, # our priors 
                          stanvars = stanvars_slopes,
                          warmup = 500,#may be necessary 
                          iter = 500 +200, #doesn't take a lot of runs
                          chains = 4, # 4 chains in parallel
                          cores = 4, # on 4 CPUs
                          refresh = 0, # don't echo chain progress
                          backend = 'cmdstanr') # use cmdstanr (other compilers broken)
  }
)
if(all_plots)
{
  #main effects means converge well
  plot(full_int_slope,
       variable = 'fmu',
       regex = TRUE,
       transform = unwrap_circular_deg) 
  #main effects means converge well
  plot(full_int_slope,
       variable = '^b_kappa',
       regex = TRUE)
  #conditional kappa mostly converge except condition 2 
  #TODO check how this is calculated
  plot(full_int_slope,
       variable = '^kappa_id',
       regex = TRUE)
  #random effects mu kappas converge well in inv_softplus space
  plot(full_int_slope,
       variable = '^zkappa',
       regex = TRUE)
  # plot(full_int_slope,
  #      variable = 'zmu_id',
  #      transform = unwrap_circular_deg,
  #      nvariables = 5,
  #      ask = FALSE)
  #individual means converge well, some bimodality in posterior
  plot(full_int_slope,
       variable = '^zmu_id_condition',
       transform = unwrap_circular_deg,
       nvariables = 4,
       regex = TRUE,
       ask = FALSE)
  # plot(full_int_slope)
}

sm_vm = summary(full_int_slope, robust = TRUE)
rn_sm_vm = rownames(sm_vm$fixed)
#fairly good convergence for main effects means
sm_vm$spec_pars
sm_vm$fixed[grepl(pattern = '^kappa', x = rn_sm_vm ),]



# Extract predictions -----------------------------------------------------



## Collect fixed effects predictions -------------------------------------
#Get fixef predictions
sm_vm = summary(full_int_slope, robust = TRUE)
prms_vm = with(sm_vm,
               rbind(fixed, #the fixed effects
                     spec_pars) #generated parameters 
)
est_vm = data.frame(t(t(prms_vm)['Estimate',])) # extract just the estimate
#all draws for circular variables
#circular intercept
full_int_slope_mu_circ_draws = brms::as_draws_df(full_int_slope,
                                                  variable = 'mu_circ') 
                                    # 
uw_mu_circ = apply(X = full_int_slope_mu_circ_draws[1:4],
                  MARGIN = 2,
                  FUN = unwrap_circular
                  )


# . Collect random effects predictions ------------------------------------
Cpal = colorRampPalette(colors = c(2:6,
                                   'seagreen',
                                   # 'salmon',
                                   # 'slategray3',
                                   'orange',
                                   'navajowhite4'
))
id_cols = sample(x = Cpal(n = 20),
                 size = 20,
                 replace = FALSE)
# [1] "#E59B14" "#DF536B" "#3EB0A2" "#BA22C0" "#24B0E5" "#626078" "#8A9630"
# [8] "#8B795E" "#7DB455" "#3ACAE0"

# [1] "#2297E6" "#69923E" "#EE9E0D" "#3DAFA5" "#3FC3DF" "#57A4D9" "#BC8B35" "#9F2F9F" "#72BE53" "#B529C1"
# [11] "#DF536B" "#CD6467" "#61D04F" "#5B6673" "#FFA500" "#884190" "#58C764" "#8B795E" "#DD981A" "#8667CD"
# [21] "#28E2E5" "#447865" "#CD0BBC" "#A9885F" "#6E85D3" "#E1A10C" "#85AC57" "#2E8B57" "#4B8E4A" "#BB7663"
# [31] "#AC8543" "#2B9FD0" "#22A1E5" "#9B7F50" "#24B7E5" "#B61DAD" "#979A5B" "#CD9228" "#4FBF7A" "#A59925"
# [41] "#25C1E5" "#27D7E5" "#C39D18" "#23ACE5" "#725482" "#879631" "#9D48C7" "#34A7BA" "#26CCE5" "#46B78F"

#Get raneff predictions for kappa
# cf_vm_k = coef(full_int_slope)$ID[,,'kappa_Intercept']
cf_vm = coef(full_int_slope)
if_cf = names(coef(full_int_slope)$ID[1,1,])
cf_vm_k = coef(full_int_slope)$ID[,,grepl(pattern = '^kappa', x = if_cf )]
cf_vm_k_gh = cf_vm_k[,,'kappa_Intercept'] #intercept condition
cf_vm_k_gl = cf_vm_k[,,'kappa_BRl']+cf_vm_k_gh # add intercept condition
cf_vm_k_uh = cf_vm_k[,,'kappa_CLu']+cf_vm_k_gh
cf_vm_k_ul = cf_vm_k[,,'kappa_BRl:CLu']+cf_vm_k_gh

#We can collect linear scaled predictions for zmu
#these could potentially be wrong if the estimate is close to +-pi
    # rn = rownames(prms_vm)
    # rn_zmu = rn[grep(pattern = 'zmu',
    #                  x = rn)]
    # colnames(cf_vm) = colnames(prms_vm)[1:4]
    # #these are likely wrong for zmu
    # ran_prms_vm = rbind(cf_vm, prms_vm[rn_zmu,1:4])

#all draws
# full_int_slope_zmu_draws = brms::as_draws_df(full_int_slope,
#                                              variable = 'zmu_id') 
full_int_slope_mu_condition_draws = brms::as_draws_df(full_int_slope,
                                                       variable = 'mu_circ') 
full_int_slope_zmu_condition_draws = brms::as_draws_df(full_int_slope,
                                                       variable = '^b_zmu',
                                                       regex = TRUE) 
zmu_nms = names(full_int_slope_zmu_condition_draws)
#TODO #appears to need unwrapping?
zmu_draws_gh = full_int_slope_zmu_condition_draws[
                          ,grepl(pattern = '^[BRl:CLu].*$', #TODO better regex
                                 x = zmu_nms )]
zmu_draws_gl = full_int_slope_zmu_condition_draws[
                          ,grepl(pattern = 'BRl$', 
                                 x = zmu_nms )]
zmu_draws_uh = full_int_slope_zmu_condition_draws[
                          ,grepl(pattern = 'CLu$',  #TODO better regex
                                 x = zmu_nms )]
zmu_draws_ul = full_int_slope_zmu_condition_draws[
                          ,grepl(pattern = 'BRl:CLu$', 
                                 x = zmu_nms )]

zmu_draws_gh = apply(X = zmu_draws_gh,
                     MARGIN = 2,
                     FUN = unwrap_circular
)
zmu_draws_gl = apply(X = zmu_draws_gl,
                     MARGIN = 2,
                     FUN = unwrap_circular
)
zmu_draws_uh = apply(X = zmu_draws_uh,
                     MARGIN = 2,
                     FUN = unwrap_circular
)
zmu_draws_ul = apply(X = zmu_draws_ul,
                     MARGIN = 2,
                     FUN = unwrap_circular
)

deg_pred_gh = circular(x = deg(zmu_draws_gh +
                                         uw_mu_circ[,1]),
                    type = 'angles',
                    unit = 'degrees',
                    template = 'geographics',
                    modulo = '2pi',
                    zero = pi/2,
                    rotation = 'clock')
deg_pred_gl = circular(x = deg(zmu_draws_gl +
                                         uw_mu_circ[,2]) + 
                                     deg_pred_condition1,
                    type = 'angles',
                    unit = 'degrees',
                    template = 'geographics',
                    modulo = '2pi',
                    zero = pi/2,
                    rotation = 'clock')
deg_pred_uh = circular(x = deg(zmu_draws_uh +
                                         uw_mu_circ[,3]) + 
                                     deg_pred_condition1,
                    type = 'angles',
                    unit = 'degrees',
                    template = 'geographics',
                    modulo = '2pi',
                    zero = pi/2,
                    rotation = 'clock')
deg_pred_ul = circular(x = deg(zmu_draws_ul +
                                         uw_mu_circ[,4]) + 
                                     deg_pred_condition1,
                    type = 'angles',
                    unit = 'degrees',
                    template = 'geographics',
                    modulo = '2pi',
                    zero = pi/2,
                    rotation = 'clock')


# Plot predictions --------------------------------------------------------
angle_unit = 'degrees'
angle_rot = 'clock'

# . plot circular random effects ------------------------------------------
dt_dim = dim(cd_subs)
u_id = with(cd_subs, unique(ID))
n_indiv = length(u_id)
i = which(u_id %in% '2016.08.26.14.39')

par(pty = 's')#sometimes gets skipped? Needs to come first
# par(mar =rep(0,4),
#     mfrow = rep(x = ceiling(sqrt(n_indiv)), 
#                 times = 2) )
with(subset(x = cd_subs,
            subset = ID %in% '2016.08.26.14.39' &
              CL %in% 'g'&
              BR %in% 'h'),
     {
       plot.circular(x = circular(x = angle,
                                  type = 'angles',
                                  unit = angle_unit,
                                  template = 'geographics',
                                  modulo = '2pi',
                                  zero = pi/2,
                                  rotation = angle_rot
       ),
       sep = 2/dt_dim[1],
       stack = TRUE,
       bins = 360/5,
       col = adjustcolor(id_cols[1], alpha.f = 200/256)
       ) 
     }
)
par(new = T)
plot.circular(x = circular(x = NULL, 
                           type = 'angles',
                           unit = 'degrees',
                           template = 'geographics',
                           modulo = '2pi',
                           zero = pi/2,
                           rotation = angle_rot
),
shrink = 1.05,
axes = F
)
points.circular(x = circular(x = deg_pred_gh[,i],
                             type = 'angles',
                             unit = angle_unit,
                             template = 'geographics',
                             modulo = '2pi',
                             zero = pi/2,
                             rotation = angle_rot
),
sep = -1e-3,
stack = TRUE,
bins = 360,
col = adjustcolor(id_cols[i], alpha.f = 1/256)
)

arrows.circular(x = median.circular(circular(x = deg(uw_mu_circ[,1]),
                                             type = 'angles',
                                             unit = angle_unit,
                                             template = 'geographics',
                                             modulo = '2pi',
                                             zero = pi/2,
                                             rotation = angle_rot
)),
y = Softpl_to_meanvec(median(cf_vm_k_gh)),
length =0, 
lwd = 1,
col = adjustcolor(id_cols[i], alpha.f = 1)
)
points.circular(x = circular(x = deg(uw_mu_circ[,1]),
                             type = 'angles',
                             unit = angle_unit,
                             template = 'geographics',
                             modulo = '2pi',
                             zero = pi/2,
                             rotation = angle_rot
),
sep = -1e-3,
stack = TRUE,
bins = 360,
col = adjustcolor(id_cols[i], alpha.f = 1/256)
) 

#these are wrong, because of the unbalanced nature of the dataset!

