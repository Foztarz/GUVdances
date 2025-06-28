# Details ---------------------------------------------------------------
#       AUTHOR:	James Foster              DATE: 2025 06 05
#     MODIFIED:	James Foster              DATE: 2025 06 11
#
#  DESCRIPTION: Functions for Foster JJ, Kühn N. & Pfeiffer K (in prep).
#               
#       INPUTS: 
#               
#      OUTPUTS: 
#
#	   CHANGES: - Stan functions removed
#             - 
#
#   REFERENCES: Edrich, W., Neumeyer, C. and von Helversen, O. (1979).
#               “Anti-sun orientation” of bees with regard to a field of ultraviolet light. 
#               J. Comp. Physiol. 134, 151–157.
#
#               Rossel, S. and Wehner, R. (1984).
#               Celestial orientation in bees: the use of spectral cues.
#               Journal of Comparative Physiology A 155, 605–613.
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
SelectCSV = function(path = file.path(ltp,'Documents', "*.csv"))
{
  # set path to file
  if(sys_win){#choose.files is only available on Windows
    message('\n\nPlease select the ".csv" file\n\n')
    Sys.sleep(0.5)#goes too fast for the user to see the message on some computers
    path_file  <- choose.files(
      default = path,#For some reason this is not possible in the "root" user
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
}

# General functions -----------------------------------------------------

#Open file with default program on any OS
# https://stackoverflow.com/a/35044209/3745353
shell.exec.OS = function(x){
  # replacement for shell.exec (doesn't exist on MAC)
  if (exists("shell.exec",where = "package:base"))
  {return(base::shell.exec(x))}else
  {comm <- paste0('open "',x,'"')
  return(system(comm))}
}

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

# Visual pigment template ------------------------------------

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



# Data summary functions --------------------------------------------------

#function for calculating mean vector params
Mvec = function(angle, units = 'radians', rotation = 'counter', type = 'angles')
{
  return(
    c(
      mu = as.numeric(mean.circular(circular(angle,units = units, rotation = rotation, type = type),
                                    control.circular = list(units = units, rotation = rotation, type = type))), 
      rho = rho.circular(circular(angle,units = units, rotation = rotation, type = type))
    )
  )
}


MLE_est = function(x)
{
  with(mle.vonmises(circular(x,
                             template = 'none'),
                    bias = TRUE),
       c(mu = as.numeric(mu), 
         kappa = kappa))
  
}

MuDiff = function(id, dt, cl, br, ref_cl = 'g', ref_br = 'h')
{
  with(
    subset(dt,
           ID %in% id),
    deg(
      mod_circular(
        rad(mu[colour %in% cl & brightn %in% br] -  
              mu[colour %in% ref_cl &
                   brightn %in% ref_br])
      )
    )
  )
}

QBox = function(quant,
                x = log10(4e14), 
                offset = -0.1,
                cols = c('green', 'darkgreen')
                )
{
  with(quant,
       {
         polygon(x = abs(offset)*c(-1,1,1,-1)+x+offset,
                 y = c(q1, q1, q3, q3),
                 col = 'white',
                 border = 'white'
         )
         polygon(x = abs(offset)*c(-1,1,1,-1)+x+offset,
                 y = c(q1, q1, q3, q3),
                 col = adjustcolor(cols[1], alpha.f = 25/256),
                 border = cols[1]
         )
         lines(x =  abs(offset)*c(-1,1)+x+offset,
               y = c(q2, q2),
               lwd = 5,
               lend = 'butt',
               col = cols[2])
         lines(x =  offset*c(1,1)+x,
               y = c(q0, q1),
               lwd = 2,
               lend = 'butt',
               col = 'black')
         lines(x =  offset*c(1,1)+x,
               y = c(q3, q4),
               lwd = 2,
               lend = 'butt',
               col = 'black')
       }
  )
}
# Circular plotting functions ---------------------------------------------


Plt_mvec = function(id)
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
}

Plt_kappa = function(id)
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
}

Plt_iskappa = function(id)
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
}

Plt_mu = function(id)
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




Plt_cmle = function(dt, col = 'black', title = '')
{
  cdt = circular(x = dt,
                 units = 'degrees',
                 rotation = 'clock',
                 zero = pi/2)
  plot_circMLE(data = cdt,
               table = circ_mle(data = cdt),
               col = c(col, col, 'gray20', 'gray20') )
  text(x = 0, y = -1.5,
       labels = title)
}

SunDiff = function(id, dt, cl, br)
{
  with(
    subset(dt,
           ID %in% id),
    deg(
      mod_circular(
        rad(mu[colour %in% cl & brightn %in% br] +  
              unique(sun_az) ) # sun azimuth in degrees
      )
    )
  )
}

IndCond = function(id, dt)
{
  with(
    subset(dt,
           ID %in% id),
    {
      length(
        unique( paste(colour, brightn) ) #find unique combinations of colour and brightness
      )
    }
  )
}

OpenCplot = function(x,
                     angle_unit = 'degrees',
                     angle_rot = 'clock',
                     n_sample = 10
)
{
  plot.circular(x = circular(x = NULL, 
                             type = 'angles',
                             unit = angle_unit,
                             # template = 'geographics',
                             modulo = '2pi',
                             zero = pi/2,
                             rotation = angle_rot
  ),
  col = NA
  )
  lines.circular(x = circular(x = 
                                seq(from = -180,
                                    to = 180,
                                    length.out = 1e3),
                              units = angle_unit, 
                              rotation = angle_rot),  
                 y = rep(x = sqrt(-log(0.05)/n_sample)-1,times = 1e3),
                 col = 'black', 
                 lty = 2,
                 lwd = 0.25,
                 lend = 'butt')
}

DiffArc = function(a1, a2, r1, r2, 
                   col1 = 'cyan4',
                   col2 = 'darkblue',
                   angle_rot = 'clock',
                   ...)
{
  ma1 = Mod360.180(a1)
  ma2 = Mod360.180(a2)
  
  if(ma2 - ma1 < 180) #fix for excess arc lengths
  {
    lines.circular(x = circular(x = seq(from = ma1, 
                                        to = ma2,
                                        length.out =1e2), 
                                type = 'angles',
                                unit = 'degrees',
                                modulo = '2pi',
                                zero = pi/2,
                                rotation = angle_rot),
                   y =  seq(from = r1, 
                            to = r2,
                            length.out = 1e2
                   )-1,
                   ...
    )
  }else
  {
    lines.circular(x = circular(x = seq(from = ma1+360, 
                                        to = ma2,
                                        length.out =1e2), 
                                type = 'angles',
                                unit = 'degrees',
                                modulo = '2pi',
                                zero = pi/2,
                                rotation = angle_rot),
                   y =  seq(from = r1, 
                            to = r2,
                            length.out = 1e2
                   )-1,
                   ...
    )
  }
  points(x = c(sin(rad(a1)) * r1,
               sin(rad(a2)) * r2
  ), 
  y = c(cos(rad(a1)) * r1,
        cos(rad(a2)) * r2
  ), 
  col = adjustcolor(col = c(col1,col2),
                    alpha.f = 0.5),
  pch= c(21, 19),
  lwd = 2
  )
}

# Set up circular formats -----------------------------------------------
pipi0 = list(units = 'radians',
             type = 'angles',
             modulo = '2pi',
             zero = 0,
             rotation = 'clock',
             template = 'none')

deg360 = list(units = 'degrees',
              type = 'angles',
              modulo = '2pi',
              zero = 0,
              rotation = 'clock',
              template = 'none')




# Generate circular parameters based on input data ------------------------
#Generate parameters for each individual
ParGenerator = function(i,
                        params,
                        nn = 1)
{
  list(mu = switch(EXPR = params$mu,#mu can be "uniform" or "vonmises" (i.e. distribution)
                   uniform = runif(
                     n = nn,
                     min = 0,
                     max = 360-.Machine$double.xmin),
                   vonmises = circular::rvonmises(
                     n = nn, 
                     mu = as.circular(x = 0,
                                      control.circular = deg360),
                     # kappa = (params$sd*pi/180)^-2),#kappa ~= 1/(sd)^2
                     kappa = A1inv(exp(((params$sd*pi/180)^2)/(-2))),
                     control.circular = deg360),#sd = sqrt(-2 * log(r))
                   runif(
                     n = nn,
                     min = 0,
                     max = 360-.Machine$double.xmin)
  ),
  kappa = exp( rep(x = log(params$kappa), times = nn) + 
                 rnorm(n = nn, mean = 0, sd = params$sd_logkappa) )
  )
}

#generate full distribution from those parameters
RVMgenerator = function(param,
                        nn = 20,
                        cc = list(units = 'degrees',
                                  type = 'angles',
                                  modulo = '2pi',
                                  zero = 0,
                                  rotation = 'clock',
                                  template = 'none'))
{
  xx= circular::rvonmises(n = nn,
                          mu = circular::as.circular(x = param$mu,
                                                     control.circular = cc), 
                          kappa = param$kappa)-pi
  return(xx)
}


# Maximum Likelihood Modelling --------------------------------------------

#Extract the parameters of the top model fit by circ_mle
MD_extract = function(md)
{
  #extract the order of models in the circm_mle results
  md_order = with(md, {rownames(results)})
  #extract the results for just the "best model" lowest AIC
  md_best = with(md, {results[md_order %in% bestmodel,]})
  #extract just the relevant parameters
  with(md_best,
       {
         return(
           list(mu1 = deg(q1), #convert to degrees
                kappa1 = k1,
                mu2 = deg(q2), #convert to degrees
                kappa2 = k2,
                weight1 = lamda,
                loglikelihood = -Likelihood #convert to true log likelihood
           )
         )
       }
  )
}

#function to calculate the likelihood of a sample, given the input ML parameters
LLcalc = function(ml,
                  angles,
                  au = 'degrees',
                  ar = 'clock',
                  ...)
{
  with(ml,
       {
         sum( # add together
           dvonmises(x = circular(x = angles,
                                  units = au,
                                  rotation = ar), # probability density for each observed angle
                     mu = mu, # ML estimated mean
                     kappa = kappa, # ML estimated concentration
                     log = TRUE) # on a log scale (i.e. add instead of multiplying)
         )
       }
  )
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
CI_vM = function(angles, #vector of angles fitted (used for sample size)
                 m1, #primary mean
                 k1, #primary concentration
                 m2 = NA, #secondary mean (ignored if NULL or NA)
                 k2 = NA, #secondary kappa
                 w1 = 1, #weighting of primary mean
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
        Mod360.180(
          quantile.circular(x = circular(x = m1_est,
                                         units = au,
                                         rotation = ar),
                            probs = probs1)
        )
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
        list(m1 = Mod360.180(
          quantile.circular(x = circular(x = m1_est,
                                         units = au,
                                         rotation = ar),
                            probs = probs1)
        ),
        m2 = Mod360.180(
          quantile.circular(x = circular(x = m2_est,
                                         units = au,
                                         rotation = ar),
                            probs = probs2)
        )
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

Table_vCI = function(ci)
{
  rowCI = function(i)
  {
    ci_lab = c('lower', 'median', 'upper')
    if(length(i) == 3) #only m1 simulated, CI extremes and median returned
    {
      nm1 = as.numeric(i)
      m2 = c(NA, NA, NA) #set m2 to NA
      # names(nm1) = names(i) #give same labels as m1 for consistency
      # names(m2) = names(i) #give same labels as m1 for consistency      
      names(nm1) = ci_lab #give same labels as m1 for consistency
      names(m2) = ci_lab #give same labels as m1 for consistency
      return(c(m1 = nm1, #remove circular formatting to collapse vector
               m2 = m2))
    }else
    { #both m1 and m2 simulated
      return(
        with(i,
             {
              nm1 = as.numeric(m1)
              nm2 = as.numeric(m2)
              # names(nm1) = names(m1)
              # names(nm2) = names(m2)
              names(nm1) = ci_lab
              names(nm2) = ci_lab
              c(m1 = nm1,#remove circular formatting to collapse vector
                m2 = nm2) #remove circular formatting to collapse vector
             }
             )
      )
    }
  }
  lst = lapply(X = ci,
         FUN = rowCI)
  do.call(what = rbind,
          lst)
}

PlotCI_vM = function(ci_vec,
                     col = 'salmon',
                     lwd = 2,
                     radius = 0.95,
                     ...)#passed to lines()
{
 angle_seq1 = 
   c(
     seq(from = ci_vec[1], #lower
         to = ci_vec[1] +
           Mod360.180(ci_vec[2]-ci_vec[1]), #median
         length.out =1e2/2),
     seq(from = ci_vec[2], #median
         to = ci_vec[2] +
           Mod360.180(ci_vec[3]-ci_vec[2]) , #upper
         length.out =1e2/2)
   )
 lines(x = radius*sin( rad(angle_seq1) ),
       y = radius*cos( rad(angle_seq1) ),
       col = col,
       lwd = 2,
       lend = 'butt',
       ...
 )
 if(!is.na(ci_vec[4]))
 {
   angle_seq2 = 
     c(
       seq(from = ci_vec[1+3],
           to = ci_vec[1+3] +
             Mod360.180(ci_vec[2+3]-ci_vec[1+3]),
           length.out =1e2/2),
       seq(from = ci_vec[2+3],
           to = ci_vec[2+3] +
             Mod360.180(ci_vec[3+3]-ci_vec[2+3]) ,
           length.out =1e2/2)
     )
 lines(x = radius*sin( rad(angle_seq2) ),
       y = radius*cos( rad(angle_seq2) ),
       col = col,
       lwd = 2,
       lend = 'butt',
       ...
 )
 }
}

PlotMV_circMLE = function(mod_par,
                          au = 'degrees',
                          ar = 'clock',
                          col1 = 'salmon',
                          col2 = col1,
                          ...) #passed to arrows.circular()
{
  m1 = with(mod_par, 
            {
              circular(mu1,
                              unit = au,
                              rotation = ar,
                              modulo = '2pi',
                              zero = pi/2
                              )
            }
            )
  with(mod_par,
       {
         arrows.circular(x = m1,
                         y = A1(kappa1),
                         col = col1,
                         lwd = 5*weight1,
                         length = 0.1,
                         ...
                         )
       }
  )
  if(!is.na(mod_par$mu2))
  {
  m2 = with(mod_par, 
            {
              circular(mu2,
                              unit = au,
                              rotation = ar,
                              modulo = '2pi',
                              zero = pi/2
                              )
            }
            )
  with(mod_par,
       {
         arrows.circular(x = m2,
                         y = A1(kappa2),
                         col = col2,
                         lwd = 5*(1-weight1),
                         length = 0.1,
                         ...
         )
       }
  )
  }
  
}

CollectDetails = function(ll,
                          dts = 'all',
                          bimod = FALSE,# whether a bimodal version is also considered
                          ...)
{
  if(!bimod)
  {
  with( subset(ll, dataset == dts),
        {
          data.frame(
            modnm = c('pairs same','pairs diff'),
            ll = c(
                   loglikelihood[model %in% 'same'],
                   loglikelihood[model %in% 'diff']
                   ),
            deviance = c(
                         -2*loglikelihood[model %in% 'same'],
                         -2*loglikelihood[model %in% 'diff']
                      ), 
            rnk = rank(
                        c(
                          -2*loglikelihood[model %in% 'same'],
                          -2*loglikelihood[model %in% 'diff'] 
                      ) ),#paired differences use fewer observations, don't include in ranking 
            df = c(1, 2)
          )
        }
  )
  }else
  {
    with( subset(ll, dataset == dts),
        {
            data.frame(
              modnm = c('pairs same','pairs diff', 'pairs multi'),
              ll = c(
                loglikelihood[model %in% 'same'],
                loglikelihood[model %in% 'diff'],
                loglikelihood[model %in% 'multi']
              ),
              deviance = c(
                -2*loglikelihood[model %in% 'same'],
                -2*loglikelihood[model %in% 'diff'],
                -2*loglikelihood[model %in% 'multi']
              ), 
              rnk = rank(
                c(
                  -2*loglikelihood[model %in% 'same'],
                  -2*loglikelihood[model %in% 'diff'], 
                  -2*loglikelihood[model %in% 'multi']
                ) ),#paired differences use fewer observations, don't include in ranking 
              df = c(1,
                     2,
                     5)
          )
        }
    )
  }
}


#function to prepare and perform the LR test
LR_calc = function(tst, mdt)
{
  #collect the deviances for h0 and h1 and calculate the difference in degrees of freedom
  lr_res = 
    with(mdt,
         {
           switch(EXPR = tst,
                  #test for uniformity (i.e. more comprehensive Rayleigh test)
                  uniformity = data.frame(
                    dev0 = deviance[modnm == 'uniform'],#null hypothesis: uniform distribution
                    dev1 = deviance[rnk == 1], #lowest rank is most likely model, any non-uniform distribution
                    d.f. = df[rnk == 1] #uniform has 0 degrees of freedom, test degrees of freedom are best model - 0
                  ),
                  #test for unpaired grand mean
                  trials_same_mean = data.frame(
                    dev0 = deviance[rnk == 2], #null hypothesis: trials don't differ
                    dev1 = deviance[modnm == 'trial mean'], #within trial obs. share a mean, expect lower deviance with more params
                    d.f. = df[modnm == 'trial mean'] -
                      df[rnk == 2] #2nd best fitting model
                  ),
                  #test for nonzero differences of pairs
                  pairs_diff_zero = data.frame(
                    dev0 = deviance[modnm == 'pairs same'], #null hypothesis: trials don't differ
                    dev1 = deviance[modnm == 'pairs diff'], #within trial obs. share a mean, expect lower deviance with more params
                    d.f. = df[modnm == 'pairs diff'] -
                      df[modnm == 'pairs same'] #grand mean has half the number of params
                  ),
                  #test for nonzero multimodal differences of pairs
                  pairs_multi_zero = data.frame(
                    dev0 = deviance[modnm == 'pairs same'], #null hypothesis: trials don't differ
                    dev1 = deviance[modnm == 'pairs multi'], #within trial obs. share a mean, expect lower deviance with more params
                    d.f. = df[modnm == 'pairs multi'] -
                      df[modnm == 'pairs same'] #grand mean has half the number of params
                  ),
                  #test for multimodal rather than unimodal differences of pairs
                  pairs_multi_diff = data.frame(
                    dev0 = deviance[modnm == 'pairs diff'], #null hypothesis: trials don't differ
                    dev1 = deviance[modnm == 'pairs multi'], #within trial obs. share a mean, expect lower deviance with more params
                    d.f. = df[modnm == 'pairs multi'] -
                      df[modnm == 'pairs diff'] #grand mean has half the number of params
                  ),
                  
           )
         }
    )
  #calculate the change in deviance (chi-squared distributed)
  lr_res = within(lr_res,
                  {
                    chi_squared = abs(unlist(dev0) - unlist(dev1))
                  }
  )
  #calculate the p value
  lr_res = within(lr_res,
                  {
                    p = pchisq(q = unlist(chi_squared),
                               df = unlist(d.f.),
                               lower.tail = FALSE)
                  }
  )
  #adjust for multiple comparisons (3 in this case)
  lr_res = within(lr_res,
                  {
                    p_adjusted = p.adjust(p = p,
                                          method = 'BH',
                                          n = length(p))#could this be flexible?
                  }
  )
  return(lr_res)
  
}
#function to prepare and perform the LR test
LR_calc_lst = function(tst, mdt, digits = 6)
{
  #collect the deviances for h0 and h1 and calculate the difference in degrees of freedom
  lr_res = 
    with(mdt,
         { data.frame(
           t(
           sapply(X = tst,
                  FUN = switch,
                  #test for uniformity (i.e. more comprehensive Rayleigh test)
                  uniformity = data.frame(
                    dev0 = deviance[modnm == 'uniform'],#null hypothesis: uniform distribution
                    dev1 = deviance[rnk == 1], #lowest rank is most likely model, any non-uniform distribution
                    d.f. = df[rnk == 1] #uniform has 0 degrees of freedom, test degrees of freedom are best model - 0
                  ),
                  #test for unpaired grand mean
                  trials_same_mean = data.frame(
                    dev0 = deviance[rnk == 2], #null hypothesis: trials don't differ
                    dev1 = deviance[modnm == 'trial mean'], #within trial obs. share a mean, expect lower deviance with more params
                    d.f. = df[modnm == 'trial mean'] -
                      df[rnk == 2] #2nd best fitting model
                  ),
                  #test for nonzero differences of pairs
                  pairs_diff_zero = data.frame(
                    dev0 = deviance[modnm == 'pairs same'], #null hypothesis: trials don't differ
                    dev1 = deviance[modnm == 'pairs diff'], #within trial obs. share a mean, expect lower deviance with more params
                    d.f. = df[modnm == 'pairs diff'] -
                      df[modnm == 'pairs same'] #grand mean has half the number of params
                  ),
                  #test for nonzero multimodal differences of pairs
                  pairs_multi_zero = data.frame(
                    dev0 = deviance[modnm == 'pairs same'], #null hypothesis: trials don't differ
                    dev1 = deviance[modnm == 'pairs multi'], #within trial obs. share a mean, expect lower deviance with more params
                    d.f. = df[modnm == 'pairs multi'] -
                      df[modnm == 'pairs same'] #grand mean has half the number of params
                  ),
                  #test for multimodal rather than unimodal differences of pairs
                  pairs_multi_diff = data.frame(
                    dev0 = deviance[modnm == 'pairs diff'], #null hypothesis: trials don't differ
                    dev1 = deviance[modnm == 'pairs multi'], #within trial obs. share a mean, expect lower deviance with more params
                    d.f. = df[modnm == 'pairs multi'] -
                      df[modnm == 'pairs diff'] #grand mean has half the number of params
                  ),
           )
           )
         )
         }
    )
  #calculate the change in deviance (chi-squared distributed)
  lr_res =  within(lr_res,
                      expr = 
                  {
                    chi_squared = abs(unlist(dev0) - unlist(dev1))
                  }
  )
  
  #calculate the p value
  lr_res = within(lr_res,
                  {
                    p = pchisq(q = unlist(chi_squared),
                               df = unlist(d.f.),
                               lower.tail = FALSE)
                  }
  )
  #adjust for multiple comparisons (3 in this case)
  lr_res = within(lr_res,
                  {
                    p_adjusted = p.adjust(p = p,
                                          method = 'BH',
                                          n = length(p))#could this be flexible?
                  }
  )
  lr_res = apply(X = lr_res,
                 MARGIN = 2,
                 FUN = unlist)#for whatever reason entries are still in list format
  return(lr_res)
  
}

#Interpret the hypothesis tests depending on size and sign of difference in deviance
H1label = function(tst, d0, d1, pa)
{
  switch(EXPR = tst,
         uniformity = if(d0 > d1 & pa <0.05) # data may be oriented, not significantly oriented, or significantly disoriented
         {'data are significantly oriented'}else
         {'data _are not_ significantly oriented'},
         trials_same_mean = if(d0 > d1 & pa <0.05) # trials significantly differ in mean, not significantly differ in mean, or share a significant mean
         {'trial means differ significantly'}else
         {'trial means _do not_ differ significantly'},
         pairs_diff_zero = if(d0 > d1 & pa <0.05) # trials significantly differ in mean, not significantly differ in mean, or share a significant mean
         {'paired trials differ significantly'}else
         {'paired trials _do not_ differ significantly'},
         pairs_multi_zero = if(d0 > d1 & pa <0.05) # trials significantly differ with multiple means, not significantly differ in mean, or share a significant mean
         {'paired trials differ significantly with multiple means'}else
         {'paired trials _do not_ differ significantly with multiple means'},
         pairs_multi_diff = if(d0 > d1 & pa <0.05) # trials significantly differ with multiple means, not significantly differ in mean, or share a significant mean
         {'paired trials differ significantly with no single mean'}else
         {'paired trials _do not_ differ significantly with multiple means'}
  )
}
