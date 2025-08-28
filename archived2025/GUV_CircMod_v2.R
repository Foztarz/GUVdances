#FOR A 'CLEAN' RUN, PRESS ctrl+shift+F10 to RESTART Rstudio
graphics.off()
# Details ---------------------------------------------------------------
#       AUTHOR:	James Foster              DATE: 2025 08 24
#     MODIFIED:	James Foster              DATE: 2025 08 28
#
#  DESCRIPTION: Attempt to run a two-way interaction model on the GUV dances data
#               using the circular modulo modelling method devel. by Jake Graving.
#               Modified from GUV_CircMod_v1.R
#               
#       INPUTS: 
#               
#      OUTPUTS: Plots and test statistics
#
#	   CHANGES: - Implemented fmod link
#             - Switched to normal priors for zkappa
#             - Stacked zkappa (slower, but more logical)
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
#- Remove all unnecessary sections  +
#- Test with real data  +
#- Double-check kappa extractions +
#- Extract coefficients +
#- Try priors for speed
#- Test with full dataset
#- Try normal raneff
#- Try vM raneff
#- Try bimodal fixef


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

lambda_mix = stanvar(scode = "
vector[K_zmu] b_logit_lambda;  // lambda for each individual
                           ",
                   block = "parameters") +
                   stanvar(scode = "
  for (i in 1:K_zmu) {
    lprior += von_misesmix_lpdf(b_fmu2[4] | 0, log1p_exp(Intercept_kappa + sum(b_kappa)), pi(), log1p_exp(Intercept_kappa + sum(b_kappa)), inv_logit(b_logit_lambda[1]));
                      }       ",
                   block = "tparameters")
  

stanvars_intercepts = stan_mvm_fun + mu_gen  #+ zmu_var+ zkappa_var
stanvars_slopes = stan_mvm_fun + mu_gen  + zkappa_var_slope + 
  zmu_var_slope 
stanvars_mix = stan_mvm_fun + mu_gen  + zkappa_var_slope + 
  zmu_var_slope + lambda_mix #includes intercepts+ zkappa_var

#implementation for bimodal effects of condition?
# b_zmu[i]

# #set up a von Mises PDF that converts to modulo (adapted from BRMS default)
# von_mises3_fun = stanvar(scode = "
# real von_misesbim_lpdf(real y, real mu, real kappa, real lambda) {
#      if (kappa < 100) {
#        return von_mises_lpdf(mod_circular(y) | mu, kappa)+log(lambda) + von_mises_lpdf(mod_circular(y) | mu+pi(), kappa)+log(1.0-lambda);
#      } else {
#        return normal_lpdf(mod_circular(y) | mu, sqrt(1 / kappa)) +log(lambda) + normal_lpdf(mod_circular(y) | mu+pi(), sqrt(1 / kappa)) +log(1-lambda) ;
#      }
#    }
# ",
#                          block = 'functions')

# Input Variables ----------------------------------------------------------

all_plots = FALSE # to speed up
##  User input -----------------------------------------------------------



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


## Formula ---------------------------------------------------------------


### Unimodal -------------------------------------------------------------

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


### Bimodal --------------------------------------------------------------

#set up model fit
formula_mix = bf(
  formula = angle ~ #set up a formula for the mean angle, modulus to (-pi,pi)
    mod_circular(fmu + fmu2 + zmu), # mean angle combines fixed and random effects
  fmu ~  BR + CL + BR:CL, # fixed effects only change as a function of condition
  fmu2 ~ 0 + BR:CL, # bimodal fixed effects only for the UV-dim condition
  zmu ~  0+ID + ID:BR + ID:CL + ID:BR:CL, #random effects change as a function of individual and condition and their interaction
  kappa ~ BR + CL + BR:CL + (1 + BR + CL + BR:CL|ID), #for kappa this occurs in linear space, and BRMS can set it up automatically
  family = von_mises(link = "identity", # the mean angle will be returned as-is
                     link_kappa = 'softplus'),#kappa will be returned via the softplus link https://en.wikipedia.org/wiki/Softplus
  nl = TRUE)#to accept user-defined extra parameters (zmu) we need to treat the formula as nonlinear


## Priors ----------------------------------------------------------------
#Unimodal
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

#Bimodal
prior_mix = get_prior(formula = formula_mix,
                            data = cd_subs,
                            check = FALSE)
dim(prior_mix)
#suggests 695  possible priors!
#they will take this general structure
print(subset(prior_mix, nlpar %in% 'fmu2')['coef'])
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

#potentially good priors
# 
# "Intercept": bmb.Prior("Normal", mu=vm_prior[1], sigma=1*np.pi/180), 
# # We expect the effect of condition to be a change of 180° (2 x 90°), most priors find this fairly well
# "Cond": bmb.Prior("Normal", mu=np.pi, sigma=60*np.pi/180), #
# # No expectations about different mean headings by condition, beyond a wider spread
# "catDoLP": bmb.Prior("Normal", mu=0, sigma=30*np.pi/180), #
# # Individual-level effects for μ: #bias to high kappa appears to cause divergent transitions here
# "1|Individual": bmb.Prior(
#   # In these experiments, we have no expectation that beetles would choose the same direction in their 1st trial
#   # Informative prior used expecting SD close to that measured empirically
#   # Using a normal distribution because that converges better
#   "Normal", mu=0, sigma=bmb.Prior("HalfStudentT", nu = 3, sigma = sd_prior/4) #assume zero, but scale by ML estimate
# ),
# # Priors for the κ-model (inside the "kappa" dictionary):
# "kappa": {
#   # Fixed effects for κ:
#   # # Across two trials the MLE for this should range from ≈1.5–2.8
#   #20250327 try higher precision to make individual mean directions easier to estimate
#   "Intercept": bmb.Prior("Normal", mu=inverse_softplus(10), sigma=1.0), 
#   "catDoLP": bmb.Prior("Normal", mu=0, sigma=1.0), 
#   # Individual-level effects for κ:
#   "1|Individual": bmb.Prior(
#     "Normal", mu=0, sigma=bmb.Prior("LogNormal", mu=np.log(0.5), sigma=0.1) #Informative prior, individual differences in concentration are small
#   ),


### Unimodal -------------------------------------------------------------

prior_int_slope = within(prior_int_slope,
                         {
   #fixed effects on mean angle are von Mises distributed (von_mises3 converts estimates to modulo (-pi,pi))                              
   prior[nlpar %in% 'fmu' & coef %in% 'Intercept'] = 'normal(0, 15*pi()/180)'# strong bias to zero
   prior[nlpar %in% 'fmu' & class %in% 'b'] = 'normal(0, pi()/3)'#moderate bias to zero, no effect #too small reduces efficiency!
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
                 x = coef)] = 'von_mises3(0, log1p_exp(zkappa1+zkappa2+zkappa3+zkappa4))'
   #fixed effects on kappa are normally distributed on the softplus scale
   prior[dpar %in% 'kappa' & class %in% 'Intercept'] = 'normal(3.0, 3.0)'#weak expectation of kappa around 3 (mean vector around 0.80)
   prior[dpar %in% 'kappa' & class %in% 'b'] = 'normal(0.0, 2.0)'#expectation of condition effect around 0
   #random effects on kappa are t-distributed on the softplus scale
   prior[dpar %in% 'kappa' & class %in% 'sd'] = 'student_t(3, 0, 2.0)' #narrow prior 
                         }
)


### Bimodal -------------------------------------------------------------

prior_mix = within(prior_mix,
                         {
   #fixed effects on mean angle are von Mises distributed (von_mises3 converts estimates to modulo (-pi,pi))                              
   prior[nlpar %in% 'fmu' & coef %in% 'Intercept'] = 'normal(0, 15*pi()/180)'# strong bias to zero
   prior[nlpar %in% 'fmu' & class %in% 'b'] = 'normal(0, pi()/3)'#moderate bias to zero, no effect #too small reduces efficiency!
   prior[nlpar %in% 'fmu'  & coef %in% c("BRl:CLu")] = 'normal(pi()/2, pi()/3)'#strong bias to the rightward turns
   prior[nlpar %in% 'fmu2' & coef %in% c("BRh:CLg", "BRh:CLu", "BRl:CLg")] = 'normal(0, 1e-3)'#No effect on nearly all conditions
   prior[nlpar %in% 'fmu2' & coef %in% c("BRl:CLu")] = 'von_misesmix(0, log1p_exp(Intercept_kappa + sum(b_kappa)), pi(), log1p_exp(Intercept_kappa + sum(b_kappa)), inv_logit(quantile(b_logit_lambda,0.5)))'#Bimodal effect with change of either 0 or 180°
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
                 x = coef)] = 'von_mises3(0, log1p_exp(zkappa1+zkappa2+zkappa3+zkappa4))'
   #fixed effects on kappa are normally distributed on the softplus scale
   prior[dpar %in% 'kappa' & class %in% 'Intercept'] = 'normal(3.0, 3.0)'#weak expectation of kappa around 3 (mean vector around 0.80)
   prior[dpar %in% 'kappa' & class %in% 'b'] = 'normal(0.0, 2.0)'#expectation of condition effect around 0
   #random effects on kappa are t-distributed on the softplus scale
   prior[dpar %in% 'kappa' & class %in% 'sd'] = 'student_t(3, 0, 2.0)' #narrow prior 
                         }
)


## add extra priors for the random effects mean angles -----------------



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
    # prior_int_slope = prior_int_slope + #random effects kappas are t-distributed on a softplus scale
    #   set_prior("target += student_t_lpdf(zkappa1 | 3, 30, 20)", #expect high concentration (low variation) 
    #             check = FALSE)+
    #   set_prior("target += student_t_lpdf(zkappa2 | 3, 0, 1.0)", #expect high concentration (low variation) 
    #             check = FALSE)+ #random effects kappas are t-distributed on a softplus scale
    #   set_prior("target += student_t_lpdf(zkappa3 | 3, 0, 1.0)", #expect high concentration (low variation) 
    #             check = FALSE)+
    #   set_prior("target += student_t_lpdf(zkappa4 | 3, 0, 0.5)", #expect high concentration (low variation) 
    #             check = FALSE)

#Unimodal
prior_int_slope = prior_int_slope + #random effects kappas are t-distributed on a softplus scale
  set_prior("target += normal_lpdf(zkappa1 | 3, 3)", #expect high concentration (low variation) 
            check = FALSE)+
  set_prior("target += normal_lpdf(zkappa2 | 0, 0.5)", #expect high concentration (low variation) 
            check = FALSE)+ #random effects kappas are t-distributed on a softplus scale
  set_prior("target += normal_lpdf(zkappa3 | 0, 0.5)", #expect high concentration (low variation) 
            check = FALSE)+
  set_prior("target += normal_lpdf(zkappa4 | 0, 0.3)", #expect high concentration (low variation) 
            check = FALSE)

#Bimodal
prior_mix = prior_mix + #random effects kappas are t-distributed on a softplus scale
  set_prior("target += normal_lpdf(zkappa1 | 3, 3)", #expect high concentration (low variation) 
            check = FALSE)+
  set_prior("target += normal_lpdf(zkappa2 | 0, 0.5)", #expect high concentration (low variation) 
            check = FALSE)+ #random effects kappas are t-distributed on a softplus scale
  set_prior("target += normal_lpdf(zkappa3 | 0, 0.5)", #expect high concentration (low variation) 
            check = FALSE)+
  set_prior("target += normal_lpdf(zkappa4 | 0, 0.3)", #expect high concentration (low variation) 
            check = FALSE)+
  set_prior("target += normal_lpdf(logit_lambda | 0, 1)", #expect high concentration (low variation) 
            check = FALSE)


### Check they can be written out ---------------------------------------------------



#unimodal
sc = make_stancode(formula = formula_int_slope,
                   data = cd,
                   prior = prior_int_slope,
                   stanvars = stanvars_slopes)
write.table(x = sc,
            file = file.path(dirname(path_file),
                             'sc_CircMod_v2.stan'),
            quote = FALSE,
            col.names = FALSE,
            row.names = FALSE)

#bimodal
sc_mix = make_stancode(formula = formula_mix,
                   data = cd,
                   prior = prior_mix,
                   stanvars = stanvars_slopes)
write.table(x = sc_mix,
            file = file.path(dirname(path_file),
                             'sc_CircMod_mix.stan'),
            quote = FALSE,
            col.names = FALSE,
            row.names = FALSE)


## Dummy run to check the influence of the priors ------------------


#double check that the prior distribution is viable by first setting up a short dummy run
# Dummy run
#Warning takes a long time to compile!
#TODO work out why this samples less efficiently than with data
# system.time( #currently takes about <1 minutes for 1000 iterations
#   {
#     dummy_int_slope = brm( formula = formula_int_slope, # using our nonlinear formula
#                            data = cd_subs, # our data
#                            prior = prior_int_slope, # our priors
#                            stanvars = stanvars_slopes,
#                            sample_prior = 'only', #ignore the data to check the influence of the priors
#                            iter = 1000, # can only estimate with enough iterations for params
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
# }


# Subset run --------------------------------------------------------------


## Unimodal --------------------------------------------------------------

# subset run
system.time(#takes less than 4 minutes for 10 individuals
  {
    full_int_slope = brm( formula = formula_int_slope, # using our nonlinear formula
                          data = cd_subs, # our data
                          prior = prior_int_slope, # our priors 
                          stanvars = stanvars_slopes,
                          warmup = 1000,#may be necessary 
                          iter = 1000+1000, #doesn't take a lot of runs
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
  plot(full_int_slope,
       variable = '^kappa_id',
       regex = TRUE)
  #random effects mu kappas converge well in inv_softplus space
  plot(full_int_slope,
       variable = '^zkappa',
       regex = TRUE)
  #individual means converge well, some bimodality in posterior
  plot(full_int_slope,
       variable = '^zmu_id_condition',
       transform = unwrap_circular_deg,
       nvariables = 4,
       regex = TRUE,
       ask = FALSE)
  #individual kappas
  plot(full_int_slope,
       variable = '^sd_ID__kappa',
       regex = TRUE)
  # plot(full_int_slope)
}

sm_vm = summary(full_int_slope, robust = TRUE)
rn_sm_vm = rownames(sm_vm$fixed)
#fairly good convergence for main effects means
sm_vm$spec_pars
sm_vm$fixed[grepl(pattern = '^kappa', x = rn_sm_vm ),]

## Bimodal --------------------------------------------------------------

# subset run
system.time(#takes 22 minutes for 10 individuals
  {
    full_mix = brm( formula = formula_mix, # using our nonlinear formula
                          data = cd_subs, # our data
                          prior = prior_mix, # our priors 
                          stanvars = stanvars_mix,
                          warmup = 500,#may be necessary 
                          iter = 500+200, #doesn't take a lot of runs
                          chains = 4, # 4 chains in parallel
                          cores = 4, # on 4 CPUs
                          refresh = 0, # don't echo chain progress
                          backend = 'cmdstanr') # use cmdstanr (other compilers broken)
  }
)
if(all_plots)
{
  #main effects means converge well
  plot(full_mix,
       variable = 'fmu',
       regex = TRUE,
       nvariables = 8,
       transform = unwrap_circular_deg) 
  #weighting parameter
  plot(full_mix,
       variable = 'lambda',
       regex = TRUE,
       transform = plogis) 
  #main effects means converge well
  plot(full_mix,
       variable = '^b_kappa',
       regex = TRUE)
  #conditional kappa mostly converge except condition 2 
  plot(full_mix,
       variable = '^kappa_id',
       regex = TRUE)
  #random effects mu kappas converge well in inv_softplus space
  plot(full_mix,
       variable = '^zkappa',
       regex = TRUE)
  #individual means converge well, some bimodality in posterior
  plot(full_mix,
       variable = '^zmu_id_condition',
       transform = unwrap_circular_deg,
       nvariables = 4,
       regex = TRUE,
       ask = FALSE)
  #individual kappas
  plot(full_mix,
       variable = '^sd_ID__kappa',
       regex = TRUE)
  # plot(full_int_slope)
}

sm_mix = summary(full_mix, robust = TRUE)
rn_sm_mix = rownames(sm_mix$fixed)
#fairly good convergence for main effects means
sm_mix$spec_pars
sm_mix$fixed[grepl(pattern = '^kappa', x = rn_sm_mix ),]



# Extract predictions -----------------------------------------------------



## Collect fixed effects predictions -------------------------------------
#Get fixef predictions
sm_vm = summary(full_int_slope, robust = TRUE)
sm_mix = summary(full_mix, robust = TRUE)

prms_vm = with(sm_vm,
               rbind(fixed, #the fixed effects
                     spec_pars) #generated parameters 
)
prms_mix = with(sm_mix,
               rbind(fixed, #the fixed effects
                     spec_pars) #generated parameters 
)
est_vm = data.frame(t(t(prms_vm)['Estimate',])) # extract just the estimate
est_mix = data.frame(t(t(prms_vm)['Estimate',])) # extract just the estimate
#all draws for circular variables
#circular fixed effects
full_int_slope_mu_circ_draws = brms::as_draws_df(full_int_slope,
                                                  variable = 'mu_circ')
full_mix_mu_circ_draws = brms::as_draws_df(full_mix,
                                                  variable = 'mu_circ')
full_mix_fmu2_draws = brms::as_draws_df(full_mix,
                                                  variable = '^b_fmu2',
                                        regex = TRUE)
#kappa effects
full_int_slope_kappa_draws = brms::as_draws_df(full_int_slope,
                                               variable = 'kappa',
                                               regex = TRUE)
full_mix_kappa_draws = brms::as_draws_df(full_mix,
                                               variable = 'kappa',
                                               regex = TRUE)
#invert link
#unwrap link
uw_mu_circ = apply(X = full_int_slope_mu_circ_draws[1:4],
                  MARGIN = 2,
                  FUN = unwrap_circular
                  )
uw_mu_circ_mix = apply(X = full_mix_mu_circ_draws[1:4],
                  MARGIN = 2,
                  FUN = unwrap_circular
                  )
uw_fmu2_mix = apply(X = full_mix_fmu2_draws[1:4],
                  MARGIN = 2,
                  FUN = unwrap_circular
                  )
#softplus link
gh_kappa = softplus(full_int_slope_kappa_draws[,1])
gl_kappa = softplus(apply(X =
                            full_int_slope_kappa_draws[,1:2],
                          FUN = sum,
                          MARGIN = 1))
uh_kappa = softplus(apply(X = 
                            full_int_slope_kappa_draws[,c(1,3)],
                          FUN = sum,
                          MARGIN = 1))
ul_kappa = softplus(apply(X = 
                            full_int_slope_kappa_draws[,c(1:4)],
                          FUN = sum,
                          MARGIN = 1))

gh_kappa_mix = softplus(full_mix_kappa_draws[,1])
gl_kappa_mix = softplus(apply(X =
                            full_mix_kappa_draws[,1:2],
                          FUN = sum,
                          MARGIN = 1))
uh_kappa_mix = softplus(apply(X = 
                            full_mix_kappa_draws[,c(1,3)],
                          FUN = sum,
                          MARGIN = 1))
ul_kappa_mix = softplus(apply(X = 
                            full_mix_kappa_draws[,c(1:4)],
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
cf_mix = coef(full_mix)

if_cf = names(coef(full_int_slope)$ID[1,1,])
cf_vm_k = coef(full_int_slope)$ID[,,grepl(pattern = '^kappa', x = if_cf )]
cf_mix_k = coef(full_mix)$ID[,,grepl(pattern = '^kappa', x = if_cf )]

cf_vm_k_gh = cf_vm_k[,,'kappa_Intercept'] #intercept condition
cf_vm_k_gl = cf_vm_k[,,'kappa_BRl']+cf_vm_k_gh # add intercept condition
cf_vm_k_uh = cf_vm_k[,,'kappa_CLu']+cf_vm_k_gh
cf_vm_k_ul = cf_vm_k[,,'kappa_BRl:CLu']+cf_vm_k_gh

cf_mix_k_gh = cf_mix_k[,,'kappa_Intercept'] #intercept condition
cf_mix_k_gl = cf_mix_k[,,'kappa_BRl']+cf_mix_k_gh # add intercept condition
cf_mix_k_uh = cf_mix_k[,,'kappa_CLu']+cf_mix_k_gh
cf_mix_k_ul = cf_mix_k[,,'kappa_BRl:CLu']+cf_mix_k_gh

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

full_mix_mu_condition_draws = brms::as_draws_df(full_mix,
                                                       variable = 'mu_circ') 
full_mix_zmu_condition_draws = brms::as_draws_df(full_mix,
                                                       variable = '^b_zmu',
                                                       regex = TRUE) 
full_mix_lambda_draws = brms::as_draws_df(full_mix,
                                                       variable = '^logit_lambda',
                                                       regex = TRUE) 

zmu_nms = names(full_int_slope_zmu_condition_draws)
ul_nms = grepl(pattern = 'BRl:CLu$', 
               x = zmu_nms )
gl_nms = grepl(pattern = 'BRl$', 
               x = zmu_nms ) & !(ul_nms)
uh_nms = grepl(pattern = 'CLu$', 
               x = zmu_nms ) & !(ul_nms)
gh_nms = !(gl_nms | uh_nms | ul_nms |
             grepl(pattern = 'chain|iteration|draw', 
                                            x = zmu_nms ))

#Sort according to name
zmu_draws_gh = full_int_slope_zmu_condition_draws[
                          ,gh_nms]
zmu_draws_gl = full_int_slope_zmu_condition_draws[
                          ,gl_nms]
zmu_draws_uh = full_int_slope_zmu_condition_draws[
                          ,uh_nms]
zmu_draws_ul = full_int_slope_zmu_condition_draws[
                          ,ul_nms]

zmu_mix_gh = full_mix_zmu_condition_draws[
                          ,gh_nms]
zmu_mix_gl = full_mix_zmu_condition_draws[
                          ,gl_nms]
zmu_mix_uh = full_mix_zmu_condition_draws[
                          ,uh_nms]
zmu_mix_ul = full_mix_zmu_condition_draws[
                          ,ul_nms]

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

zmu_mix_gh = apply(X = zmu_mix_gh,
                     MARGIN = 2,
                     FUN = unwrap_circular
)
zmu_mix_gl = apply(X = zmu_mix_gl,
                     MARGIN = 2,
                     FUN = unwrap_circular
)
zmu_mix_uh = apply(X = zmu_mix_uh,
                     MARGIN = 2,
                     FUN = unwrap_circular
)
zmu_mix_ul = apply(X = zmu_mix_ul,
                     MARGIN = 2,
                     FUN = unwrap_circular
)

deg_pred_gh = circular(x = deg(zmu_draws_gh[,1:dim(zmu_draws_gh)[2]] +
                                         uw_mu_circ[,1] ),
                    type = 'angles',
                    unit = 'degrees',
                    template = 'geographics',
                    modulo = '2pi',
                    zero = pi/2,
                    rotation = 'clock')
deg_pred_gl = circular(x = deg(zmu_draws_gl[,1:dim(zmu_draws_gl)[2]] +
                                         uw_mu_circ[,2]) + 
                         deg_pred_gh,
                    type = 'angles',
                    unit = 'degrees',
                    template = 'geographics',
                    modulo = '2pi',
                    zero = pi/2,
                    rotation = 'clock')
deg_pred_uh = circular(x = deg(zmu_draws_uh[,1:dim(zmu_draws_uh)[2]] +
                                         uw_mu_circ[,3]) + 
                         deg_pred_gh,
                    type = 'angles',
                    unit = 'degrees',
                    template = 'geographics',
                    modulo = '2pi',
                    zero = pi/2,
                    rotation = 'clock')
deg_pred_ul = circular(x = deg(zmu_draws_ul[,1:dim(zmu_draws_uh)[2]] +
                                         uw_mu_circ[,4]) + 
                         deg(zmu_draws_uh[,1:dim(zmu_draws_uh)[2]] +
                               uw_mu_circ[,3]) +
                         deg(zmu_draws_gl[,1:dim(zmu_draws_gl)[2]] +
                               uw_mu_circ[,2]) + 
                         deg_pred_gh,
                    type = 'angles',
                    unit = 'degrees',
                    template = 'geographics',
                    modulo = '2pi',
                    zero = pi/2,
                    rotation = 'clock')

#bimodal
mix_pred_gh = circular(x = deg(zmu_mix_gh[,1:dim(zmu_mix_gh)[2]] +
                                         uw_mu_circ_mix[,1] ),
                    type = 'angles',
                    unit = 'degrees',
                    template = 'geographics',
                    modulo = '2pi',
                    zero = pi/2,
                    rotation = 'clock')
mix_pred_gl = circular(x = deg(zmu_mix_gl[,1:dim(zmu_mix_gl)[2]] +
                                         uw_mu_circ_mix[,2]) + 
                         mix_pred_gh,
                    type = 'angles',
                    unit = 'degrees',
                    template = 'geographics',
                    modulo = '2pi',
                    zero = pi/2,
                    rotation = 'clock')
mix_pred_uh = circular(x = deg(zmu_mix_uh[,1:dim(zmu_mix_uh)[2]] +
                                         uw_mu_circ_mix[,3]) + 
                                        mix_pred_gh,
                    type = 'angles',
                    unit = 'degrees',
                    template = 'geographics',
                    modulo = '2pi',
                    zero = pi/2,
                    rotation = 'clock')
mix_pred_ul1 = circular(x = deg(zmu_mix_ul[,1:dim(zmu_mix_uh)[2]] +
                                         uw_mu_circ_mix[,4]) + 
                         deg(zmu_mix_uh[,1:dim(zmu_mix_uh)[2]] +
                               uw_mu_circ_mix[,3]) +
                         deg(zmu_mix_gl[,1:dim(zmu_mix_gl)[2]] +
                               uw_mu_circ_mix[,2]) + 
                         mix_pred_gh,
                    type = 'angles',
                    unit = 'degrees',
                    template = 'geographics',
                    modulo = '2pi',
                    zero = pi/2,
                    rotation = 'clock')
mix_pred_ul2 = circular(x = deg(zmu_mix_ul[,1:dim(zmu_mix_uh)[2]] +
                                         uw_mu_circ_mix[,4]) + 
                         deg(zmu_mix_uh[,1:dim(zmu_mix_uh)[2]] +
                               uw_mu_circ_mix[,3]) +
                         deg(zmu_mix_gl[,1:dim(zmu_mix_gl)[2]] +
                               uw_mu_circ_mix[,2]) + 
                               deg(apply(uw_fmu2_mix, MARGIN = 1, FUN = sum)) + 
                         mix_pred_gh,
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
par(pty = 's')#sometimes gets skipped? Needs to come first
par(mar =rep(0,4),
    # mfrow = rep(x = ceiling(sqrt(n_indiv)),
    #             times = 2) )
    mfrow = c(n_indiv,
                4) )
for(id in u_id)
{
  # i = which(u_id %in% '2016.07.12.14.43')
  i = grep(x = colnames(deg_pred_gh),
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
  
  # arrows.circular(x = median.circular(circular(x = deg(uw_mu_circ[,1]),
  arrows.circular(x = median.circular(circular(x = deg_pred_gh[,i],
                                               type = 'angles',
                                               unit = angle_unit,
                                               template = 'geographics',
                                               modulo = '2pi',
                                               zero = pi/2,
                                               rotation = angle_rot
  )),
  y = Softpl_to_meanvec(est_vm$kappa_Intercept +
                                 cf_vm_k_gh[id, 'Estimate']),
  length =0, 
  lwd = 1,
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
  points.circular(x = circular(x = deg_pred_gl[,i],
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
  
  # arrows.circular(x = median.circular(circular(x = deg(uw_mu_circ[,1]),
  arrows.circular(x = median.circular(circular(x = deg_pred_gl[,i],
                                               type = 'angles',
                                               unit = angle_unit,
                                               template = 'geographics',
                                               modulo = '2pi',
                                               zero = pi/2,
                                               rotation = angle_rot
  )),
  y = Softpl_to_meanvec(
                        est_vm$kappa_Intercept + est_vm$kappa_BRl +
                         cf_vm_k_gl[id, 'Estimate'] ),
  length =0, 
  lwd = 1,
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
  points.circular(x = circular(x = deg_pred_uh[,i],
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
  
  # arrows.circular(x = median.circular(circular(x = deg(uw_mu_circ[,1]),
  arrows.circular(x = median.circular(circular(x = deg_pred_uh[,i],
                                               type = 'angles',
                                               unit = angle_unit,
                                               template = 'geographics',
                                               modulo = '2pi',
                                               zero = pi/2,
                                               rotation = angle_rot
  )),
  y = Softpl_to_meanvec(est_vm$kappa_Intercept + est_vm$kappa_CLu +
                          cf_vm_k_uh[id, 'Estimate']),
  length =0, 
  lwd = 1,
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
  points.circular(x = circular(x = deg_pred_ul[,i],
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
  
  # arrows.circular(x = median.circular(circular(x = deg(uw_mu_circ[,1]),
  arrows.circular(x = median.circular(circular(x = deg_pred_ul[,i],
                                               type = 'angles',
                                               unit = angle_unit,
                                               template = 'geographics',
                                               modulo = '2pi',
                                               zero = pi/2,
                                               rotation = angle_rot
  )),
  y = Softpl_to_meanvec(
                        est_vm$kappa_Intercept + est_vm$kappa_BRl +
                         est_vm$kappa_CLu + est_vm$kappa_BRl.CLu +
                           cf_vm_k_ul[id, 'Estimate']
                        ),
  length =0, 
  lwd = 1,
  col = adjustcolor(id_cols[i], alpha.f = 1)
  )
}


# Plot coefficients -------------------------------------------------------

par(
    mfrow = c(2, 4),
    mar = c(4,4,4,4)
    )

PCestimates(angles =  deg(uw_mu_circ[,1]),
      col = 'green3',
      title = 'Green Bright (Intercept)',
      titline = 2
)
PCestimates(angles =  deg(uw_mu_circ[,2]),
      col = 'darkgreen',
      title = 'Green Dim - Green Bright',
      titline = 2
)
PCestimates(angles =  deg(uw_mu_circ[,4]),
      col = 'purple',
      title = 'UV Dim - Green Bright',
      titline = 2
)
PCestimates(angles =  deg(uw_mu_circ[,3]),
      col = 'gray40',
      title = 'UV Bright - Green Bright',
      titline = 2
)

VertHist(data = unlist(gh_kappa),
         ylab = 'kappa',
         breaks = 2e1,
         main = 'Green Bright (Intercept)',
         col = 'green3'
         )

VertHist(data = unlist(gl_kappa - gh_kappa),
         ylab = 'Δkappa',
         breaks = 2e1,
         ylim = c(-1, 1)*8,
         col = 'darkgreen',
         main = 'Green Dim - Green Bright')
abline(h = 0,
       col = 'gray')

VertHist(data = unlist(ul_kappa - gh_kappa),
         ylab = 'Δkappa',
         breaks = 2e1,
         ylim = c(-1, 1)*8,
         col = 'purple',
         main = 'UV Dim - Green Bright')
abline(h = 0,
       col = 'gray')

VertHist(data = unlist(uh_kappa - gh_kappa),
         ylab = 'Δkappa',
         breaks = 2e1,
         ylim = c(-1, 1)*8,
         col = 'gray40',
         main = 'UV Bright - Green Bright')
abline(h = 0,
       col = 'gray')

#rho version
par(
    mfrow = c(2, 4),
    mar = c(4,4,4,4)
    )

PCestimates(angles =  deg(uw_mu_circ[,1]),
      col = 'green3',
      title = 'Green Bright (Intercept)'
)
PCestimates(angles =  deg(uw_mu_circ[,2]),
      col = 'darkgreen',
      title = 'Green Dim - Green Bright'
)
PCestimates(angles =  deg(uw_mu_circ[,4]),
      col = 'purple',
      title = 'UV Dim - Green Bright'
)
PCestimates(angles =  deg(uw_mu_circ[,3]),
      col = 'gray40',
      title = 'UV Bright - Green Bright'
)

VertHist(data = A1(unlist(gh_kappa)),
         ylab = 'rho',
         ylim = c(0,1),
         breaks = 2e1,
         main = 'Green Bright (Intercept)',
         col = 'green3'
         )

VertHist(data = A1(unlist(gl_kappa)) - A1(unlist(gh_kappa)),
         ylab = 'Δrho',
         breaks = 2e1,
         ylim = c(-1, 1)*0.5,
         col = 'darkgreen',
         main = 'Green Dim - Green Bright')
abline(h = 0,
       col = 'gray')

VertHist(data = A1(unlist(ul_kappa)) - A1(unlist(gh_kappa)),
         ylab = 'Δrho',
         breaks = 2e1,
         ylim = c(-1, 1)*0.5,
         col = 'purple',
         main = 'UV Dim - Green Bright')
abline(h = 0,
       col = 'gray')

VertHist(data = A1(unlist(uh_kappa)) - A1(unlist(gh_kappa)),
         ylab = 'Δrho',
         breaks = 2e1,
         ylim = c(-1, 1)*0.5,
         col = 'gray40',
         main = 'UV Bright - Green Bright')
abline(h = 0,
       col = 'gray')




# Bimodal model -----------------------------------------------------------

par(pty = 's')#sometimes gets skipped? Needs to come first
par(mar =rep(0,4),
    # mfrow = rep(x = ceiling(sqrt(n_indiv)),
    #             times = 2) )
    mfrow = c(n_indiv,
              4) )
for(id in u_id)
{
  i = grep(x = colnames(mix_pred_gh),
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
  points.circular(x = circular(x = mix_pred_gh[,i],
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
  
  # arrows.circular(x = median.circular(circular(x = deg(uw_mu_circ[,1]),
  arrows.circular(x = median.circular(circular(x = mix_pred_gh[,i],
                                               type = 'angles',
                                               unit = angle_unit,
                                               template = 'geographics',
                                               modulo = '2pi',
                                               zero = pi/2,
                                               rotation = angle_rot
  )),
  y = Softpl_to_meanvec(est_mix$kappa_Intercept +
                          cf_mix_k_gh[id, 'Estimate']),
  length =0, 
  lwd = 1,
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
  points.circular(x = circular(x = mix_pred_gl[,i],
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
  
  # arrows.circular(x = median.circular(circular(x = deg(uw_mu_circ[,1]),
  arrows.circular(x = median.circular(circular(x = mix_pred_gl[,i],
                                               type = 'angles',
                                               unit = angle_unit,
                                               template = 'geographics',
                                               modulo = '2pi',
                                               zero = pi/2,
                                               rotation = angle_rot
  )),
  y = Softpl_to_meanvec(
    est_mix$kappa_Intercept + est_mix$kappa_BRl +
      cf_mix_k_gl[id, 'Estimate'] ),
  length =0, 
  lwd = 1,
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
  points.circular(x = circular(x = mix_pred_uh[,i],
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
  
  # arrows.circular(x = median.circular(circular(x = deg(uw_mu_circ[,1]),
  arrows.circular(x = median.circular(circular(x = mix_pred_uh[,i],
                                               type = 'angles',
                                               unit = angle_unit,
                                               template = 'geographics',
                                               modulo = '2pi',
                                               zero = pi/2,
                                               rotation = angle_rot
  )),
  y = Softpl_to_meanvec(est_mix$kappa_Intercept + est_mix$kappa_CLu +
                          cf_mix_k_uh[id, 'Estimate']),
  length =0, 
  lwd = 1,
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
  points.circular(x = circular(x = mix_pred_ul1[,i],
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
  points.circular(x = circular(x = mix_pred_ul2[,i],
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
  
  # arrows.circular(x = median.circular(circular(x = deg(uw_mu_circ[,1]),
  arrows.circular(x = median.circular(circular(x = mix_pred_ul1[,i],
                                               type = 'angles',
                                               unit = angle_unit,
                                               template = 'geographics',
                                               modulo = '2pi',
                                               zero = pi/2,
                                               rotation = angle_rot
  )),
  y = Softpl_to_meanvec(
    est_mix$kappa_Intercept + est_mix$kappa_BRl +
      est_mix$kappa_CLu + est_mix$kappa_BRl.CLu +
      cf_mix_k_ul[id, 'Estimate']
  ),
  length =0, 
  lwd = 1*plogis(median(full_mix_lambda_draws$logit_lambda)),
  col = adjustcolor(id_cols[i], alpha.f = 1)
  )
  # arrows.circular(x = median.circular(circular(x = deg(uw_mu_circ[,1]),
  arrows.circular(x = median.circular(circular(x = mix_pred_ul2[,i],
                                               type = 'angles',
                                               unit = angle_unit,
                                               template = 'geographics',
                                               modulo = '2pi',
                                               zero = pi/2,
                                               rotation = angle_rot
  )),
  y = Softpl_to_meanvec(
    est_mix$kappa_Intercept + est_mix$kappa_BRl +
      est_mix$kappa_CLu + est_mix$kappa_BRl.CLu +
      cf_mix_k_ul[id, 'Estimate']
  ),
  length =0, 
  lwd = 1*plogis(-median(full_mix_lambda_draws$logit_lambda)),
  col = adjustcolor(id_cols[i], alpha.f = 1)
  )
}
  


