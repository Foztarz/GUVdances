# Details ---------------------------------------------------------------
#       AUTHOR:	James Foster              DATE: 2025 06 05
#     MODIFIED:	James Foster              DATE: 2025 06 05
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
                 title = '')
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
        side = 1,
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


