#FOR A 'CLEAN' RUN, PRESS ctrl+shift+F10 to RESTART Rstudio
graphics.off()
# Details ---------------------------------------------------------------
#       AUTHOR:	James Foster              DATE: 2025 09 02
#     MODIFIED:	James Foster              DATE: 2025 09 04
#
#  DESCRIPTION: Attempt to run a two-way interaction model on the GUV dances data
#               using the circular modulo modelling method devel. by Jake Graving.
#               Includes a secondary mean to account for bimodality for UV-dim stimulus.
#               Modified from GUV_CircMod_v2.R
#               
#       INPUTS: 
#               
#      OUTPUTS: Plots and test statistics
#
#	   CHANGES: - 
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
#- Remove all unnecessary sections  
#- Make custom family +
#- Generate mod circular variables + 
#- Extract bimodal effects +
#- Left-right bimodal +
#- Think about priors for opposite means
#- Constrain mu slope SD to improve bimodality identifiability  +
#- Try initialising at 0  +
#- Check extraction of ul condition
#- Test von Mises version
#- Sorted version of bimodal
#- Try priors for speed
#- Move functions out of script
#- Test with full dataset


# Set up workspace --------------------------------------------------------
angle_unit = 'degrees'
angle_rot = 'clock'

## Load packages ----------------------------------------------------------
#needs installing before first use (in Rstudio, see automatic message)
suppressMessages(#these are disturbing users unnecessarily
  {
    require(circular)#package for handling circular data
    require(CircMLE)#package for circular mixture models
    require(brms)#package for preparing Stan models
  }
)


## Plot spacing function -------------------------------------------------
#TODO check usage
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
  # atan2(y = sin(x),
  #       x = cos(x))
  ( (x + pi) %% 2*pi ) - pi
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
#any data, but plotted as a histogram on a vertical rather than horizontal axis
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

#TODO check usage
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

PCfun = function(angles,
                 col,
                 shrink = 1.5,
                 title = '',
                 side = 1)
{
  ca = circular(x = angles,
                units = 'degrees',
                rotation = 'clock')
  plot.circular(x = ca,
                col = col,
                stack = TRUE,
                sep = 0.1,
                bins = 355/5,
                units = 'degrees',
                rotation = 'clock',
                zero = pi/2,
                shrink = shrink)
  mtext(text = title,
        side = side,
        line = -2)
  lines(x = c(0,0),
        y = c(-1,1),
        col = 'gray')
  arrows.circular(x = mean.circular(ca),
                  y = rho.circular(ca),
                  zero = pi/2,
                  rotation = 'clock',
                  col = col,
                  length =0.1)
}


#generic mean angle simulator
MeanRvm = function(n, #representative sample size
                   mu = circular(0), #mean (defaults to 0rad)
                   kappa, #kappa required
                   au = 'degrees', #units
                   ar = 'clock') #rotation direction
{
  mean.circular(rvonmises(n = n, 
                          mu = circular(mu, units = au, rotation = ar), 
                          kappa = kappa,
                          control.circular = list(units = au, rotation = ar)))
}


#Simulate confidence intervals for a unimodal or bimodal distribution
#fitted to a vector of "angles"
#Simulate confidence intervals for a unimodal or bimodal distribution
#fitted to a vector of "angles"
CI_vM = function(angles, #vector of angles fitted (used for sample size)
                 m1, #primary mean
                 k1, #primary concentration
                 m2 = NA, #secondary mean (ignored if NULL or NA)
                 k2 = NA, #secondary kappa
                 w1 = 1, #weighting of primary mean
                 force_mu = FALSE, #force median at true mu?
                 n = 1e4, #number of simulations
                 au = 'degrees', 
                 ar = 'clock',
                 calc_q = TRUE,
                 alternative = 'one.sided', #two.sided less conservative
                 interval = 0.95, #confidence interval to calculate
                 speedup_parallel = TRUE
)
{
  if(speedup_parallel) #3x faster
  {
    cl = parallel::makePSOCKcluster(parallel::detectCores()-1)
    parallel::clusterExport(cl = cl, 
                            varlist = c('mean.circular',
                                        'circular',
                                        'rvonmises'),
                            envir = .GlobalEnv
    )
    parallel::clusterExport(cl = cl, 
                            varlist = c('MeanRvm',
                                        'angles',
                                        'm1',
                                        'k1',
                                        'm2',
                                        'k2',
                                        'w1',
                                        'n',
                                        'au',
                                        'ar'),
                            envir = environment()
    )
    #simulate primary mean
    m1_est = 
      parallel::parSapply(cl = cl,
                          X = 1:n,
                          FUN = function(i)
                          {
                            eval.parent(
                              {
                                MeanRvm(n = round(length(angles)*w1), #estimate number of observations at primary mean
                                        mu = m1, 
                                        kappa = k1,
                                        au = au,
                                        ar = ar)
                              }
                            )
                          },
                          simplify = 'array' #return an array of simulated angles
      )
    if(!is.na(m2)) #if there is a valid secondary mean
    {
      m2_est = 
        parallel::parSapply(cl = cl,
                            X = 1:n,
                            FUN = function(i)
                            {
                              eval.parent(
                                {
                                  MeanRvm(n = round(length(angles)*(1-w1)), #estimate number of observations at secondary mean
                                          mu = m2, 
                                          kappa = k2,
                                          au = au,
                                          ar = ar)
                                }
                              )
                            },
                            simplify = 'array' #return an array of simulated angles
        )
    }
    parallel::stopCluster(cl)
  }else
  { #if not using parallel, use the slower version via replicate()
    m1_est = replicate(n = n, 
                       MeanRvm(n = round(length(angles)*w1), 
                               mu = m1, 
                               kappa = k1,
                               au = au,
                               ar = ar)
    )
    if(!is.na(m2))
    {
      m2_est = replicate(n = n, 
                         MeanRvm(n = round(length(angles)*(1-w1)), 
                                 mu = m2, 
                                 kappa = k2,
                                 au = au,
                                 ar = ar)
      )
    }
  }
  return(
    if(calc_q) #calculate quantiles only if requested
    {
      #either two-sided, symmetrical around mean change
      #or one-sided, from zero change towards mean change
      probs1 = switch(alternative,
                      two.sided = sort(c(c(0,1)+c(1,-1)*(1-interval)/2, 0.5)),
                      one.sided = sort(c(c(0,1)+
                                           (if(Mod360.180(m1)>0) #N.B. quantile.circular counts anticlockwise
                                           {c(1,0)}else
                                           {c(0,-1)}
                                           )*(1-interval), 0.5)),
                      sort(c(c(0,1)+ #default to one-sided
                               (if(Mod360.180(m1)>0)
                               {c(1,0)}else
                               {c(0,-1)}
                               )*(1-interval), 0.5))
      )
      if(is.na(m2))
      {
        if(force_mu)
        {
          Mod360.180(
            quantile( Mod360.180(as.numeric(m1_est) - m1),
                      probs = probs1) + m1
          )
        }else
        {
          Mod360.180(
            quantile.circular(x = circular(x = m1_est,
                                           units = au,
                                           rotation = ar),
                              probs = probs1)
          )
        }
      }else
      {
        probs2 = switch(alternative,
                        two.sided = sort(c(c(0,1)+c(1,-1)*(1-interval)/2, 0.5)),
                        one.sided = sort(c(c(0,1)+
                                             (if(Mod360.180(m2)>0)
                                             {c(1,0)}else
                                             {c(0,-1)}
                                             )*(1-interval), 0.5)),
                        sort(c(c(0,1)+ #default to one-sided
                                 (if(Mod360.180(m2)<0)
                                 {c(1,0)}else
                                 {c(0,-1)}
                                 )*(1-interval), 0.5))
        )
        list(m1 = if(force_mu)
        {
          Mod360.180(
            quantile( Mod360.180(as.numeric(m1_est) - m1),
                      probs = probs1) + m1
          )
        }else
        {
          Mod360.180(
            quantile.circular(x = circular(x = m1_est,
                                           units = au,
                                           rotation = ar),
                              probs = probs1)
          )
        },
        m2 = if(force_mu)
        {
          Mod360.180(
            quantile( Mod360.180(as.numeric(m2_est) - m2),
                      probs = probs2) + m2
          )
        }else
        {
          Mod360.180(
            quantile.circular(x = circular(x = m2_est,
                                           units = au,
                                           rotation = ar),
                              probs = probs2)
          )
        }
        )
      }
    }else
    { #if quantiles not requested, return the simulations (mainly for troubleshooting)
      if(is.na(m2))
      {
        m1_est = 
          sapply(X = m1_est, FUN = Mod360.180)
      }else
      {
        list(
          m1_est = 
            sapply(X = m1_est, FUN = Mod360.180),     
          m2_est = 
            sapply(X = m2_est, FUN = Mod360.180),
        )
      }
    }
  )
}


PlotCI_vM = function(ci_vec,
                     col = 'salmon',
                     lwd = 2,
                     radius = 0.95,
                     ...)#passed to lines()
{
  ci_vec = as.numeric(ci_vec)#remove circular formatting!
  #changed on 20250815, plotting issues near median
  angle_seq1.1 = 
    seq(from = ci_vec[1], #lower
        to = ci_vec[1] +
          Mod360.180(ci_vec[2]-ci_vec[1]), #median
        length.out =1e2/2)
  angle_seq1.2 = 
    seq(from = ci_vec[2], #median
        to = ci_vec[2] +
          Mod360.180(ci_vec[3]-ci_vec[2]) , #upper
        length.out =1e2/2)
  lines(x = radius*sin( rad(angle_seq1.1) ),
        y = radius*cos( rad(angle_seq1.1) ),
        col = col,
        lwd = lwd,
        lend = 'butt',
        ...
  )
  lines(x = radius*sin( rad(angle_seq1.2) ),
        y = radius*cos( rad(angle_seq1.2) ),
        col = col,
        lwd = lwd,
        lend = 'butt',
        ...
  )
  if(!is.na(ci_vec[4]))
  {
    #changed on 20250815
    angle_seq2.1 = 
      seq(from = ci_vec[1+3],
          to = ci_vec[1+3] +
            Mod360.180(ci_vec[2+3]-ci_vec[1+3]),
          length.out =1e2/2)
    
    angle_seq2.2 = 
      seq(from = ci_vec[2+3],
          to = ci_vec[2+3] +
            Mod360.180(ci_vec[3+3]-ci_vec[2+3]) ,
          length.out =1e2/2)
    lines(x = radius*sin( rad(angle_seq2.1) ),
          y = radius*cos( rad(angle_seq2.1) ),
          col = col,
          lwd = lwd,
          lend = 'butt',
          ....)
    lines(x = radius*sin( rad(angle_seq2.1) ),
          y = radius*cos( rad(angle_seq2.1) ),
          col = col,
          lwd = lwd,
          lend = 'butt',
          ...)
  }
}

#TODO check usage
#plot circular model estimates
PCestimates = function(angles,
                       col = 'darkblue',
                       adjcol = 1/256,
                       sep = 2/length(unlist(angles)),
                       shrink = 1.5,
                       lw = 3,
                       title = '',
                       titline = -2,
                       side = 1,
                       add_mv = TRUE,
                       add_ci = FALSE,
                       add_cred = TRUE,
                       force_mu = TRUE,
                       interval = 0.95,
                       angle_unit = 'degrees',
                       angle_rot = 'clock'
)
{
  ca = circular(unlist(angles),
                units = angle_unit,
                rot = angle_rot
  )
  plot.circular(x = ca,
                col = adjustcolor(col,
                                  alpha.f = 1/256),
                stack = TRUE,
                sep = sep,
                bins = 355/5,
                units = angle_unit,
                rotation = angle_rot,
                zero = pi/2,
                shrink = shrink)
  mtext(text = title,
        side = side,
        line = titline)
  lines(x = c(0,0),
        y = c(-1,1),
        col = 'gray')
  if(add_mv)
  {
    arrows.circular(x = mean.circular(ca),
                    y = rho.circular(ca),
                    zero = pi/2,
                    rotation = angle_rot,
                    col = col,
                    lwd = lw,
                    length =0.1)
  }
  if(add_ci)
  {
    mlevm = mle.vonmises(ca, bias = TRUE)
    
    #calculate vector of estimates
    
    civ = with(mlevm,
               {
                 CI_vM(angles = ca,
                       m1 = mu,
                       k1 = kappa,
                       alternative = 'two.sided',
                       force_mu = if(kappa == 0 ){TRUE}else{FALSE})
               }
    )
  }
  if(add_cred)
  {
    if(force_mu)
    {
      civ =  Mod360.180(
        quantile( Mod360.180(as.numeric(ca) - mean.circular(ca)),
                  ,
                  probs = sort(c(c(0,1)+c(1,-1)*(1-interval)/2, 0.5))) + mean.circular(ca)
      )
    }else
    {
      civ = Mod360.180(
        quantile.circular(x = ca,
                          probs = sort(c(c(0,1)+c(1,-1)*(1-interval)/2, 0.5)))
      )
    }
  }
  #plot quantiles of estimates
  if(add_ci | add_cred)
  {
    PlotCI_vM(ci_vec = civ,
              col = col, 
              lwd = lw,
              # radius = 1+250*sep*shrink)
              radius = 1+0.5)
  }
}




## Stan variables ---------------------------------------------------------


### Functions ------------------------------------------------------------


  # #set up the modulo function
  # mod_circular_fun = stanvar(scode = "
  #   real mod_circular(real y) {
  #     return atan2(sin(y), cos(y));
  #   }
  # ",
#set up the modulo function
# From Stan manual
#"Warning: These functions can seriously hinder sampling and optimization efficiency for gradient-based methods (e.g., NUTS, HMC, BFGS) if applied to parameters (including transformed parameters and local variables in the transformed parameters or model block). The problem is that they break gradients due to discontinuities coupled with zero gradients elsewhere. They do not hinder sampling when used in the data, transformed data, or generated quantities blocks."
mod_circular_fun = stanvar(scode = "
  real mod_circular(real y) {
    return fmod(y + pi(), 2*pi()) - pi();
  }
",
                           block = 'functions')

#TODO check usage
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

#TODO check usage
#set up a von Mises PDF that converts to modulo (adapted from BRMS default)
# a mixture model of the form PDF = vonmises1*lambda + vonmises2*(1-lambda)
von_misesmix_fun = stanvar(scode = "
real von_misesmix_lpdf(real y, real mu1, real kappa1, real mu2, real kappa2, real lambda) {
     if (kappa1 < 100 && kappa2 < 100) {
       return log( exp(von_mises_lpdf(mod_circular(y) | mu1, kappa1)+log(lambda)) + exp(von_mises_lpdf(mod_circular(y) | mu2, kappa2)+log(1-lambda)) );
     } else {
       return log( exp(normal_lpdf(mod_circular(y) | mu1, sqrt(1/kappa1))+log(lambda)) + exp(normal_lpdf(mod_circular(y) | mu2, 1/sqrt(kappa2))+log(1-lambda)) );
     }
   }
",
                         block = 'functions')

#TODO check usage
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

#TODO check usage
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
  von_mises3_fun +
  von_misesmix_fun

### Parameters -----------------------------------------------------------


#generate modulo outputs of fixed effects
mu_gen = stanvar(scode = "
vector [M_1] mu_circ; //modulo circular estimate
for (i in 1:size(b_fmu)){
mu_circ[i] = mod_circular(b_fmu[i]);
}
",
                 block = 'genquant')

#random effects on mean angle
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
real kappa_id_condition4 = log1p_exp(zkappa1+zkappa2+zkappa3+zkappa4);
          ", 
          block = 'genquant')

#TODO check usage
lambda_mix = stanvar(scode = "
real Intercept_logit_lambda;  // lambda for each individual
vector[K_zmu] b_logit_lambda;  // lambda for each individual
                           ",
                   block = "parameters") +
                   stanvar(scode = "
  for (i in 1:K_zmu) {
    lprior += von_misesmix_lpdf(b_fmu2[4] | 0, log1p_exp(Intercept_kappa + sum(b_kappa)), pi(), log1p_exp(Intercept_kappa + sum(b_kappa)), inv_logit(Intercept_logit_lambda + b_logit_lambda[1]));
                      }       ",
                   block = "tparameters")
  

stanvars_intercepts = stan_mvm_fun + mu_gen  #+ zmu_var+ zkappa_var
stanvars_slopes = stan_mvm_fun + mu_gen  + zkappa_var_slope + 
  zmu_var_slope 
stanvars_mix = stan_mvm_fun + mu_gen  + zkappa_var_slope + 
  zmu_var_slope 


# Input Variables ----------------------------------------------------------
all_plots = FALSE # to speed up

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
                                    size =  10,
                                    replace = FALSE)
                )
# cd_subs = cd

## Custom family --------------------------------------------------------
#In the "unwrap" family, all variables 
unwrap_von_mises = custom_family(
  "unwrap_von_mises", dpars = c("mu", "kappa"),
  links = c('identity',#brms cannot accept custom link functions, do via nl instead
            "softplus"), 
  lb = c(-pi, 0), ub = c(pi, NA),
  type = "real"
)

stan_unwrap_fun = stanvar(scode = "
  real unwrap_von_mises_lpdf(real y, real mu, real kappa) {
    return von_mises_lpdf(y | mod_circular(mu), kappa);
  }
  real unwrap_von_mises_rng(real mu, real kappa) {
    return von_mises_rng( mod_circular(mu) , kappa);
  }
",
  block = 'functions') + 
  stanvar(scode = "
  real inv_mod_circular(real y) {
    return mod_circular(y);
  }
",
  block = 'functions') 

#generate modulo outputs of fixed effects
unwrap_mu_gen = stanvar(scode = "
  vector [Kc_mu1 +1] mu_circ1; //modulo circular estimate
  vector [Kc_mu2 +1] mu_circ2; //modulo circular estimate
  mu_circ1[1] = mod_circular(Intercept_mu1); //Intercept case
  mu_circ2[1] = mod_circular(Intercept_mu2); //Intercept case
  for (i in 1:Kc_mu1){
  mu_circ1[i+1] = mod_circular(b_mu1[i]);
  mu_circ2[i+1] = mod_circular(b_mu2[i]);
  }
",
block = 'genquant')

#post processing
log_lik_unwrap_von_mises  = function(i, prep) {
  mu = brms::get_dpar(prep, "mu", i = i)
  phi = brms::get_dpar(prep, "kappa", i = i)
  y = prep$data$Y[i]
  unwrap_von_mises_lpdf(y, mu, kappa)
}

von_mises_rng = function(mu, kappa)
{
  cc = circular::rvonmises(n = 1,
                      mu = circular(x = mu, template = 'none'),
                      kappa = kappa)
  return(as.numeric(cc))
}
posterior_predict_unwrap_von_mises = function(i, draws, ...) {
  mu = brms::get_dpar(draws, "mu", i = i)
  phi  =  brms::get_dpar(draws, "kappa", i = i)
  von_mises_rng(mu, kappa)
}

posterior_epred_unwrap_von_mises = function(draws, ...)
  {
    mu  =  mod_circular(brms::get_dpar(draws, "mu"))
    kappa  =  brms::get_dpar(draws, "kappa")
  }



## Formula ---------------------------------------------------------------
#N.B. Was limited to just the interaction effects for 2nd distribution, now includes all
#set up model fit
formula_mix = bf(#modulus may not be necessary, included in lpd function
  formula = angle ~ mu,#set up a formula for the mean angles, modulus to (-pi,pi)
              mu1 ~ BR + CL + BR:CL + (1 + BR + CL + BR:CL|ID), # mean angle combines fixed and random effects
              mu2 ~ BR + CL + BR:CL + (1 + BR + CL + BR:CL|ID), # mean angle effect of condition combination
              kappa1 ~ BR + CL + BR:CL + (1 + BR + CL + BR:CL|ID), #for kappa this occurs in linear space, and BRMS can set it up automatically
              kappa2 ~ BR + CL + BR:CL + (1 + BR + CL + BR:CL|ID), #for 2nd mean, this is specific to individual and condition combination
            theta1 ~ BR + CL + BR:CL + (1 + BR + CL + BR:CL|ID), #for mixture weighting, this is specific to individual and condition combination
  family = mixture(unwrap_von_mises,#mod mu, kappa via the softplus
                   unwrap_von_mises
                   ),#mixture, specified by the log ratio of theta1 : theta2
  nl = FALSE)#to accept user-defined extra parameters (zmu) we need to treat the formula as nonlinear
## Priors ----------------------------------------------------------------

#priors for mu
pr_mu_mix = 
  prior(normal(0,pi()/12), class = Intercept, dpar = 'mu1') + # anchor to 0° #TODO loosen
  prior(normal(0,pi()/12), class = Intercept, dpar = 'mu2') + # anchor to 0° #Needs to be tight if theta1 is high
  prior(normal(pi()/3,pi()/4), class = b, dpar = 'mu1') + # weak bias to right turns
  prior(normal(-pi()/3,pi()/4), class = b, dpar = 'mu2') + # weak bias to left turns
  prior(lognormal( log(pi()/9), 0.3), dpar = 'mu1', class = 'sd', group  = 'ID') + #small differences across conditions
  prior(lognormal( log(pi()/9), 0.3), dpar = 'mu2', class = 'sd', group  = 'ID') + #small differences across conditions
  prior(lognormal( log(pi()/15), 0.5), dpar = 'mu2', class = 'sd', coef = 'Intercept', group  = 'ID') + #keep sd under 90°
  prior(lognormal( log(pi()/15), 0.5), dpar = 'mu1', class = 'sd', coef = 'Intercept', group  = 'ID') #keep sd under 90°
#priors for kappa
pr_kappa_mix = 
  prior(normal(5.0,1.0), class = Intercept, dpar = 'kappa1') + # bias to oriented
  prior(normal(5.0,1.0), class = Intercept, dpar = 'kappa2') + # bias to oriented
  prior(normal(0,1.0), class = b, dpar = 'kappa1') + # strong bias to no change
  prior(normal(0,1.0), class = b, dpar = 'kappa2') + # strong bias to no change
  prior(student_t(3,0, 1.0), class = sd, dpar = 'kappa1') + # weak bias to no turn
  prior(student_t(3,0, 1.0), class = sd, dpar = 'kappa2') # weak bias to no turn
#priors for theta (mixture weight)
pr_theta_mix = 
  prior(normal(10,1.0), class = Intercept, dpar = 'theta1') + # bias to mu1 as primary
  prior(normal(0,10), class = b, dpar = 'theta1') + # weak bias to zero
  prior(student_t(3, 0, 0.5), class = sd, dpar = 'theta1') 


pr_mix = pr_mu_mix + pr_kappa_mix + pr_theta_mix


## Save Stancode -----------------------------------------------------------


sc_mix = make_stancode(formula = formula_mix,
                   data = cd_subs,
                   prior = pr_mix,
                   stanvars = stan_unwrap_fun + mod_circular_fun + unwrap_mu_gen)

write.table(x = sc_mix,
            file = file.path(dirname(path_file),
                             'sc_CircMod_mix.stan'),
            quote = FALSE,
            col.names = FALSE,
            row.names = FALSE)

## Run model -------------------------------------------------------------
wup = 1000
sam = 200

#very long compile time
system.time(
  {
    
full_mix = brm(formula = formula_mix,
         data = cd_subs,
         prior = pr_mix,
         stanvars = stan_unwrap_fun + #unwrapped lpdf
                     mod_circular_fun + #relies on circular modulus
                     unwrap_mu_gen, #add modulus to generated quantities
         warmup = wup,#may be necessary 
         iter = wup+sam, #doesn't take a lot of runs
         chains = 4, # 4 chains in parallel
         cores = 4, # on 4 CPUs
         init = 0, # could ease initialisation
         # threads = 4, # on 4 CPUs
         # open_progress = TRUE) # make a progress bar
         refresh = 100, # echo chain progress every n iterations
         silent = 0, # echo Stan messages
         backend = 'cmdstanr')

}
)

## Plot coefficients -----------------------------------------------------

#main effects means
#primary mu
plot(full_mix,
     variable = '^b_mu1', # all effects have similar names
     regex = TRUE,
     nvariables = 4,
     transform = unwrap_circular_deg) 
#secondary mu
plot(full_mix,
     variable = '^b_mu2', # all effects have similar names
     regex = TRUE,
     nvariables = 5,
     transform = unwrap_circular_deg) 
#weighting (logistic scaled)
plot(full_mix,
     variable = '^b_theta1',
     regex = TRUE) 
#primary kappa
plot(full_mix,
     variable = '^b_kappa1',
     regex = TRUE)#main effects means converge well
#secondary kappa
plot(full_mix,
     variable = '^b_kappa2',
     regex = TRUE)#main effects means converge well
#linear SD should be less than 180°
plot(full_mix,
     variable = '^sd_ID__mu1',
     regex = TRUE,
     transform = deg)
plot(full_mix,
     variable = '^sd_ID__mu2',
     regex = TRUE,
     transform = deg)
plot(full_mix,
     variable = '^sd_ID__theta1',
     regex = TRUE)
plot(full_mix,
     variable = '^mu_circ1',
     regex = TRUE,
     transform = unwrap_circular_deg)
plot(full_mix,
     variable = '^mu_circ2',
     regex = TRUE,
     transform = unwrap_circular_deg)


## Summarise coefficients ------------------------------------------------
#extract the medians of all parameters
sm_mix = summary(full_mix, robust = TRUE)
#find the names of the fixed effects
rn_sm_mix = with(sm_mix, rownames(fixed) )
#Investigate the 
sm_mix$spec_pars
sm_mix$fixed[grepl(pattern = '^mu', x = rn_sm_mix ),]
sm_mix$fixed[grepl(pattern = '^theta', x = rn_sm_mix ),]
sm_mix$fixed[grepl(pattern = '^kappa', x = rn_sm_mix ),]





# Extract predictions -----------------------------------------------------
ce = conditional_effects(full_mix,effects = 'BR:CL')#weird output

pe = posterior_predict(full_mix)#doesn't seem to work

## Collect fixed effects predictions -------------------------------------
#TODO make a function
#Get fixef predictions
sm_mix = summary(full_mix, robust = TRUE)
prms_mix = with(sm_mix,
               rbind(fixed, #the fixed effects
                     spec_pars) #generated parameters 
)
est_mix = data.frame(t(t(prms_mix)['Estimate',])) # extract just the estimate
#all draws for circular variables
#circular fixed effects
md_mu1_draws = brms::as_draws_df(full_mix,
                          variable = 'mu1',
                          regex = TRUE)
md_mu2_draws = brms::as_draws_df(full_mix,
                          variable = 'mu2',
                          regex = TRUE)
#kappa effects
md_kappa1_draws = brms::as_draws_df(full_mix,
                               variable = 'kappa1',
                               regex = TRUE)
md_kappa2_draws = brms::as_draws_df(full_mix,
                                 variable = 'kappa2',
                                 regex = TRUE)
#mixture weights
md_theta1_draws = brms::as_draws_df(full_mix,
                               variable = 'theta1',
                               regex = TRUE)

## Invert link functions -------------------------------------------------
#and generate fixef predictions for each draw

#unwrap link
#extract fixef coef for later comparison
uw_mu1 = apply(X = md_mu1_draws[1:4],
                  MARGIN = 2,
                  FUN = unwrap_circular
                  )
uw_mu2 = apply(X = md_mu2_draws[1:4],
                  MARGIN = 2,
                  FUN = unwrap_circular
                  )
#construct fixef predictions
gh_mu1 = apply(X = md_mu1_draws[,1],
                                 MARGIN = 2,
                                 FUN = unwrap_circular
                                 )
gh_mu2 = apply(X = md_mu2_draws[,1],
                                 MARGIN = 2,
                                 FUN = unwrap_circular
                                 )

gl_mu1 = unwrap_circular(apply(X = md_mu1_draws[1:2],
                           FUN = sum,
                           MARGIN = 1))
gl_mu2 = unwrap_circular(apply(X = md_mu2_draws[,1:2],
                           FUN = sum,
                           MARGIN = 1))

uh_mu1 = unwrap_circular(apply(X = md_mu1_draws[,c(1,3)],
                           FUN = sum,
                           MARGIN = 1))
uh_mu2 = unwrap_circular(apply(X = md_mu2_draws[,c(1,3)],
                           FUN = sum,
                           MARGIN = 1))

ul_mu1 = unwrap_circular(apply(X = md_mu1_draws[,1:4],
                           FUN = sum,
                           MARGIN = 1))
ul_mu2 = unwrap_circular(apply(X = md_mu2_draws[,1:4],
                           FUN = sum,
                           MARGIN = 1))


#softplus link
gh_kappa1 = softplus(md_kappa1_draws[,1])
gh_kappa2 = softplus(md_kappa2_draws[,1])
gl_kappa1 = softplus(apply(X =
                            md_kappa1_draws[,1:2],
                          FUN = sum,
                          MARGIN = 1))
gl_kappa2 = softplus(apply(X =
                            md_kappa2_draws[,1:2],
                          FUN = sum,
                          MARGIN = 1))
uh_kappa1 = softplus(apply(X = 
                            md_kappa1_draws[,c(1,3)],
                          FUN = sum,
                          MARGIN = 1))
uh_kappa2 = softplus(apply(X = 
                            md_kappa2_draws[,c(1,3)],
                          FUN = sum,
                          MARGIN = 1))
ul_kappa1 = softplus(apply(X = 
                            md_kappa1_draws[,c(1:4)],
                          FUN = sum,
                          MARGIN = 1))
ul_kappa2 = softplus(apply(X = 
                            md_kappa2_draws[,c(1:4)],
                          FUN = sum,
                          MARGIN = 1))
#inv_logit link
gh_lambda = apply(X = md_theta1_draws[,1], 
                  MARGIN = 2,
                  FUN = plogis)
gl_lambda = plogis(apply(X =
                            md_theta1_draws[,1:2],
                          FUN = sum,
                          MARGIN = 1))
uh_lambda = plogis(apply(X = 
                            md_theta1_draws[,c(1,3)],
                          FUN = sum,
                          MARGIN = 1))
ul_lambda = plogis(apply(X = 
                             md_theta1_draws[,c(1:4)],
                           FUN = sum,
                           MARGIN = 1))

## Collect random effects predictions ------------------------------------
u_id = with(cd_subs, unique(ID))
n_indiv = length(u_id)
Cpal = colorRampPalette(colors = c(2:6,
                                   'seagreen',
                                   'salmon3',
                                   'gray25',
                                   # 'slategray4',
                                   'orange3',
                                   # 'navajowhite4'
                                   'darkred',
                                   'darkblue'
))
id_cols = sample(x = Cpal(n = 20),
                 size = n_indiv,
                 replace = FALSE)

#Get raneff predictions for kappa
cf_mix = coef(full_mix)$ID

cf_mix_k1_gh = cf_mix[,,'kappa1_Intercept'] #intercept condition
cf_mix_k1_gl = cf_mix[,,'kappa1_BRl'] # don't add intercept condition!
cf_mix_k1_uh = cf_mix[,,'kappa1_CLu']
cf_mix_k1_ul = cf_mix[,,'kappa1_BRl:CLu']

cf_mix_k2_gh = cf_mix[,,'kappa2_Intercept'] #intercept condition
cf_mix_k2_gl = cf_mix[,,'kappa2_BRl']+cf_mix_k2_gh # add intercept condition
cf_mix_k2_uh = cf_mix[,,'kappa2_CLu']+cf_mix_k2_gh
cf_mix_k2_ul = apply(X = cf_mix[,,c('kappa2_Intercept',
                                    'kappa2_BRl',
                                    'kappa2_CLu',
                                    'kappa2_BRl:CLu')
                      ],
                      MARGIN = 1,
                      FUN = sum)

#Get raneff predictions for lambda

cf_mix_l_gh = cf_mix[,,'theta1_Intercept'] #intercept condition
cf_mix_l_gl = cf_mix[,,'theta1_BRl']+cf_mix_l_gh # add intercept condition
cf_mix_l_uh = cf_mix[,,'theta1_CLu']+cf_mix_l_gh
cf_mix_l_ul = apply(X = cf_mix[,,c('theta1_Intercept',
                                    'theta1_BRl',
                                    'theta1_CLu',
                                    'theta1_BRl:CLu')
                    ],
                    MARGIN = 1,
                    FUN = sum)


#find mu conditions
nm_mu = names(md_mu1_draws)
nm_ID = grep(pattern = 'r_ID',
             x = nm_mu)
#extract draws for individual
md_mu1_ID_draws = md_mu1_draws[,nm_ID]
md_mu2_ID_draws = md_mu2_draws[,nm_ID]
#find condition names
nm_gh = grep(pattern = ',Intercept]',
             x = names(md_mu1_ID_draws))
nm_gl = grep(pattern = ',BRl]',
             x = names(md_mu1_ID_draws))
nm_uh = grep(pattern = ',CLu]',
             x = names(md_mu1_ID_draws))
nm_ul = grep(pattern = ',BRl:CLu]',
             x = names(md_mu1_ID_draws))

#Sort according to name
zmu_draws1_gh = md_mu1_ID_draws[,nm_gh]
zmu_draws1_gl = md_mu1_ID_draws[
                          ,nm_gl]+zmu_draws1_gh
zmu_draws1_uh = md_mu1_ID_draws[
                          ,nm_uh]+zmu_draws1_gh
zmu_draws1_ul = md_mu1_ID_draws[
                          ,nm_ul]+zmu_draws1_gh+zmu_draws1_gl+zmu_draws1_uh
zmu_draws2_gh = md_mu2_ID_draws[,nm_gh]
zmu_draws2_gl = md_mu2_ID_draws[
                          ,nm_gl]+zmu_draws2_gh
zmu_draws2_uh = md_mu2_ID_draws[
                          ,nm_uh]+zmu_draws2_gh
zmu_draws2_ul = md_mu2_ID_draws[
                          ,nm_ul]+zmu_draws2_gh+zmu_draws2_gl+zmu_draws2_uh
#unwrap estimates
#TODO check if this is necessary
zmu_draws1_gh = apply(X = zmu_draws1_gh,
                     MARGIN = 2,
                     FUN = unwrap_circular
)
zmu_draws1_gl = apply(X = zmu_draws1_gl,
                     MARGIN = 2,
                     FUN = unwrap_circular
)
zmu_draws1_uh = apply(X = zmu_draws1_uh,
                     MARGIN = 2,
                     FUN = unwrap_circular
)
zmu_draws1_ul = apply(X = zmu_draws1_ul,
                     MARGIN = 2,
                     FUN = unwrap_circular
)
zmu_draws2_gh = apply(X = zmu_draws2_gh,
                     MARGIN = 2,
                     FUN = unwrap_circular
)
zmu_draws2_gl = apply(X = zmu_draws2_gl,
                     MARGIN = 2,
                     FUN = unwrap_circular
)
zmu_draws2_uh = apply(X = zmu_draws2_uh,
                     MARGIN = 2,
                     FUN = unwrap_circular
)
zmu_draws2_ul = apply(X = zmu_draws2_ul,
                     MARGIN = 2,
                     FUN = unwrap_circular
)


## Construct circular predictors -----------------------------------------


deg_pred1_gh = circular(x = deg(zmu_draws1_gh + c(gh_mu1)),
                    type = 'angles',
                    unit = 'degrees',
                    template = 'geographics',
                    modulo = '2pi',
                    zero = pi/2,
                    rotation = 'clock')
deg_pred1_gl = circular(x = deg(zmu_draws1_gl + c(gh_mu1)),
                    type = 'angles',
                    unit = 'degrees',
                    template = 'geographics',
                    modulo = '2pi',
                    zero = pi/2,
                    rotation = 'clock')
deg_pred1_uh = circular(x = deg(zmu_draws1_uh + c(uh_mu1)),
                    type = 'angles',
                    unit = 'degrees',
                    template = 'geographics',
                    modulo = '2pi',
                    zero = pi/2,
                    rotation = 'clock')
deg_pred1_ul = circular(x = deg(zmu_draws1_ul + c(ul_mu1)),
                    type = 'angles',
                    unit = 'degrees',
                    template = 'geographics',
                    modulo = '2pi',
                    zero = pi/2,
                    rotation = 'clock')

deg_pred2_gh = circular(x = deg(zmu_draws2_gh + c(gh_mu2)),
                    type = 'angles',
                    unit = 'degrees',
                    template = 'geographics',
                    modulo = '2pi',
                    zero = pi/2,
                    rotation = 'clock')
deg_pred2_gl = circular(x = deg(zmu_draws2_gl + c(gh_mu2)),
                    type = 'angles',
                    unit = 'degrees',
                    template = 'geographics',
                    modulo = '2pi',
                    zero = pi/2,
                    rotation = 'clock')
deg_pred2_uh = circular(x = deg(zmu_draws2_uh + c(uh_mu2)),
                    type = 'angles',
                    unit = 'degrees',
                    template = 'geographics',
                    modulo = '2pi',
                    zero = pi/2,
                    rotation = 'clock')
deg_pred2_ul = circular(x = deg(zmu_draws2_ul + c(ul_mu2)),
                    type = 'angles',
                    unit = 'degrees',
                    template = 'geographics',
                    modulo = '2pi',
                    zero = pi/2,
                    rotation = 'clock')


# Plot predictions --------------------------------------------------------
angle_unit = 'degrees'
angle_rot = 'clock'

## plot circular random effects ------------------------------------------
dt_dim = dim(cd_subs)
u_id = with(cd_subs, unique(ID))
n_indiv = length(u_id)
par(pty = 's')#sometimes gets skipped? Needs to come first
par(mar =rep(0,4),
    mfrow = c(n_indiv,
                4) )
for(id in u_id)
{
  # i = which(u_id %in% '2016.07.12.14.43')
  i = grep(x = colnames(deg_pred1_gh),
           pattern = id)
  #gh
  with(subset(x = cd_subs,
              subset = ID %in% id &
                CL %in% 'g'&
                BR %in% 'h'),
       {
         plot.circular(x = circular(x = deg(angle),
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
         col = 'green3'
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
  points.circular(x = circular(x = deg_pred1_gh[,i],
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
  col = adjustcolor('darkblue', alpha.f = 1/256)
  )
  points.circular(x = circular(x = deg_pred2_gh[,i],
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
  col = adjustcolor('darkred', alpha.f = 1/256)
  )
  #mu1 estimate
  arrows.circular(x = median.circular(circular(x = deg_pred1_gh[,i],
                                               type = 'angles',
                                               unit = angle_unit,
                                               template = 'geographics',
                                               modulo = '2pi',
                                               zero = pi/2,
                                               rotation = angle_rot
  )),
  y = Softpl_to_meanvec(cf_mix_k1_gh[id, 'Estimate']),
  length =0, 
  lwd = 5*plogis(cf_mix_l_gh[id, 'Estimate']),
  col = adjustcolor(id_cols[i], alpha.f = 1)
  )
  #mu2 estimate
  arrows.circular(x = median.circular(circular(x = deg_pred2_gh[,i],
                                               type = 'angles',
                                               unit = angle_unit,
                                               template = 'geographics',
                                               modulo = '2pi',
                                               zero = pi/2,
                                               rotation = angle_rot
  )),
  y = Softpl_to_meanvec(cf_mix_k2_gh[id, 'Estimate']),
  length =0, 
  lwd = 5*plogis(-cf_mix_l_gh[id, 'Estimate']),
  col = adjustcolor(id_cols[i], alpha.f = 1)
  )
  #gl
  with(subset(x = cd_subs,
              subset = ID %in% id &
                CL %in% 'g'&
                BR %in% 'l'),
       {
         plot.circular(x = circular(x = deg(angle),
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
         col = 'darkgreen'
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
  points.circular(x = circular(x = deg_pred1_gl[,i],
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
  col = adjustcolor('darkblue', alpha.f = 1/256)
  )
  points.circular(x = circular(x = deg_pred2_gl[,i],
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
  col = adjustcolor('darkred', alpha.f = 1/256)
  )
  #mu1 estimate
  arrows.circular(x = median.circular(circular(x = deg_pred1_gl[,i],
                                               type = 'angles',
                                               unit = angle_unit,
                                               template = 'geographics',
                                               modulo = '2pi',
                                               zero = pi/2,
                                               rotation = angle_rot
  )),
  y = Softpl_to_meanvec(cf_mix_k1_gl[id, 'Estimate']),
  length =0, 
  lwd = 5*plogis(cf_mix_l_gl[id, 'Estimate']),
  col = adjustcolor(id_cols[i], alpha.f = 1)
  )
  #mu2 estimate
  arrows.circular(x = median.circular(circular(x = deg_pred2_gl[,i],
                                               type = 'angles',
                                               unit = angle_unit,
                                               template = 'geographics',
                                               modulo = '2pi',
                                               zero = pi/2,
                                               rotation = angle_rot
  )),
  y = Softpl_to_meanvec(cf_mix_k2_gl[id, 'Estimate']),
  length =0, 
  lwd = 5*plogis(-cf_mix_l_gl[id, 'Estimate']),
  col = adjustcolor(id_cols[i], alpha.f = 1)
  )
  
  #uh
  with(subset(x = cd_subs,
              subset = ID %in% id &
                CL %in% 'u'&
                BR %in% 'h'),
       {
         plot.circular(x = circular(x = deg(angle),
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
         col = 'magenta'
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
  points.circular(x = circular(x = deg_pred1_uh[,i],
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
  col = adjustcolor('darkblue', alpha.f = 1/256)
  )
  points.circular(x = circular(x = deg_pred2_uh[,i],
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
  col = adjustcolor('darkred', alpha.f = 1/256)
  )
  #mu1 estimate
  arrows.circular(x = median.circular(circular(x = deg_pred1_uh[,i],
                                               type = 'angles',
                                               unit = angle_unit,
                                               template = 'geographics',
                                               modulo = '2pi',
                                               zero = pi/2,
                                               rotation = angle_rot
  )),
  y = Softpl_to_meanvec(cf_mix_k1_uh[id, 'Estimate']),
  length =0, 
  lwd = 5*plogis(cf_mix_l_uh[id, 'Estimate']),
  col = adjustcolor(id_cols[i], alpha.f = 1)
  )
  #mu2 estimate
  arrows.circular(x = median.circular(circular(x = deg_pred2_uh[,i],
                                               type = 'angles',
                                               unit = angle_unit,
                                               template = 'geographics',
                                               modulo = '2pi',
                                               zero = pi/2,
                                               rotation = angle_rot
  )),
  y = Softpl_to_meanvec(cf_mix_k2_uh[id, 'Estimate']),
  length =0, 
  lwd = 5*plogis(-cf_mix_l_uh[id, 'Estimate']),
  col = adjustcolor(id_cols[i], alpha.f = 1)
  )
  #ul
  with(subset(x = cd_subs,
              subset = ID %in% id &
                CL %in% 'u'&
                BR %in% 'l'),
       {
         plot.circular(x = circular(x = deg(angle),
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
         col = 'purple'
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
  points.circular(x = circular(x = deg_pred1_ul[,i],
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
  col = adjustcolor('darkblue', alpha.f = 1/256)
  )
  points.circular(x = circular(x = deg_pred2_ul[,i],
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
  col = adjustcolor('darkred', alpha.f = 1/256)
  )
  #mu1 estimate
  arrows.circular(x = median.circular(circular(x = deg_pred1_ul[,i],
                                               type = 'angles',
                                               unit = angle_unit,
                                               template = 'geographics',
                                               modulo = '2pi',
                                               zero = pi/2,
                                               rotation = angle_rot
  )),
  y = Softpl_to_meanvec(cf_mix_k1_ul[id]),
  length =0, 
  lwd = 5*plogis(cf_mix_l_ul[id]),
  col = adjustcolor(id_cols[i], alpha.f = 1)
  )
  #mu2 estimate
  arrows.circular(x = median.circular(circular(x = deg_pred2_uh[,i],
                                               type = 'angles',
                                               unit = angle_unit,
                                               template = 'geographics',
                                               modulo = '2pi',
                                               zero = pi/2,
                                               rotation = angle_rot
  )),
  y = Softpl_to_meanvec(cf_mix_k2_ul[id]),
  length =0, 
  lwd = 5*plogis(-cf_mix_l_ul[id]),
  col = adjustcolor(id_cols[i], alpha.f = 1)
  )
}


# Plot coefficients -------------------------------------------------------

par(
    mfrow = c(2, 4),
    mar = c(4,4,4,4)
    )

PCestimates(angles =  deg(uw_mu1[,1]),
      col = 'green3',
      title = 'Green Bright (Intercept)',
      titline = 2
)
PCestimates(angles =  deg(uw_mu1[,2]),
      col = 'darkgreen',
      title = 'Green Dim - Green Bright',
      titline = 2
)
PCestimates(angles =  deg(uw_mu1[,4]),
      col = 'purple',
      title = 'UV Dim - Green Bright',
      titline = 2
)
PCestimates(angles =  deg(uw_mu1[,3]),
      col = 'gray40',
      title = 'UV Bright - Green Bright',
      titline = 2
)

VertHist(data = unlist(gh_kappa1),
         ylab = 'kappa',
         breaks = 2e1,
         main = 'Green Bright (Intercept)',
         col = 'green3'
         )

VertHist(data = unlist(gl_kappa1 - gh_kappa1),
         ylab = 'Δkappa',
         breaks = 2e1,
         ylim = c(-1, 1)*8,
         col = 'darkgreen',
         main = 'Green Dim - Green Bright')
abline(h = 0,
       col = 'gray')

VertHist(data = unlist(ul_kappa1 - gh_kappa1),
         ylab = 'Δkappa',
         breaks = 2e1,
         ylim = c(-1, 1)*8,
         col = 'purple',
         main = 'UV Dim - Green Bright')
abline(h = 0,
       col = 'gray')

VertHist(data = unlist(uh_kappa1 - gh_kappa1),
         ylab = 'Δkappa',
         breaks = 2e1,
         ylim = c(-1, 1)*8,
         col = 'gray40',
         main = 'UV Bright - Green Bright')
abline(h = 0,
       col = 'gray')

#TODO update this
      # #rho version
      # par(
      #     mfrow = c(2, 4),
      #     mar = c(4,4,4,4)
      #     )
      # 
      # PCestimates(angles =  deg(uw_mu1[,1]),
      #       col = 'green3',
      #       title = 'Green Bright (Intercept)'
      # )
      # PCestimates(angles =  deg(uw_mu1[,2]),
      #       col = 'darkgreen',
      #       title = 'Green Dim - Green Bright'
      # )
      # PCestimates(angles =  deg(uw_mu1[,4]),
      #       col = 'purple',
      #       title = 'UV Dim - Green Bright'
      # )
      # PCestimates(angles =  deg(uw_mu1[,3]),
      #       col = 'gray40',
      #       title = 'UV Bright - Green Bright'
      # )
      # 
      # VertHist(data = A1(unlist(gh_kappa1)),
      #          ylab = 'rho',
      #          ylim = c(0,1),
      #          breaks = 2e1,
      #          main = 'Green Bright (Intercept)',
      #          col = 'green3'
      #          )
      # 
      # VertHist(data = A1(unlist(gl_kappa1)) - A1(unlist(gh_kappa1)),
      #          ylab = 'Δrho',
      #          breaks = 2e1,
      #          ylim = c(-1, 1)*0.5,
      #          col = 'darkgreen',
      #          main = 'Green Dim - Green Bright')
      # abline(h = 0,
      #        col = 'gray')
      # 
      # VertHist(data = A1(unlist(ul_kappa1)) - A1(unlist(gh_kappa1)),
      #          ylab = 'Δrho',
      #          breaks = 2e1,
      #          ylim = c(-1, 1)*0.5,
      #          col = 'purple',
      #          main = 'UV Dim - Green Bright')
      # abline(h = 0,
      #        col = 'gray')
      # 
      # VertHist(data = A1(unlist(uh_kappa1)) - A1(unlist(gh_kappa1)),
      #          ylab = 'Δrho',
      #          breaks = 2e1,
      #          ylim = c(-1, 1)*0.5,
      #          col = 'gray40',
      #          main = 'UV Bright - Green Bright')
      # abline(h = 0,
      #        col = 'gray')





