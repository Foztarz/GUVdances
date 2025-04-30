#FOR A 'CLEAN' RUN, PRESS ctrl+shift+F10 to RESTART Rstudio
graphics.off()
# Details ---------------------------------------------------------------
#       AUTHOR:	James Foster              DATE: 2025 03 20
#     MODIFIED:	James Foster              DATE: 2025 04 30
#
#  DESCRIPTION: Inspect and summarise individual dances.
#               
#       INPUTS: 
#               
#      OUTPUTS: Plots and test statistics
#
#	   CHANGES: - 
#
#   REFERENCES: 
#              Edrich, W., Neumeyer, C. and Von Heiversen, O. (1979).
#              “Anti-sun orientation” of bees with regard to a field of ultraviolet light. 
#              J. Comp. Physiol. 134, 151–157.
#
#              Rossel, S. and Wehner, R. (1984).
#              Celestial orientation in bees: the use of spectral cues.
#              Journal of Comparative Physiology A 155, 605–613.
#
#    EXAMPLES:  
#
# 
#TODO   ---------------------------------------------
#TODO   
#- Visualise individual dances  +
#- Summarise individual dances  +
#- Visualise condition differences  +
#- Indentify full condition individuals
#- Investigate only full conditions
#- Summarise condition differences  


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


# Inspect individual dances -----------------------------------------------

#set up a sequence for the axes
xc = seq(from = -pi, to = pi-1e-16, length.out = 1e3)
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

#set up plot
pdf_fl = file.path(dirname(path_file), 'all_dances.pdf')
pdf(file = pdf_fl)
par(mfrow = c(2, 2), mar = c(0,0,0,0))
par(pty = 's')

#add legend to 1st page
plot(x = NULL,
     xlim = c(-1,1),
     ylim = c(-1,1),
     pch = 19,
     axes = FALSE,
     xlab = '',
     ylab = '',
     main = ''
)
legend(x = 'center',
       legend = c('Green Bright',
                  'Green Dim',
                  'UV Bright',
                  'UV Dim',
                  'Sun azimuth (North = Up)'),
       col = c('green',
               'darkgreen',
               'magenta',
               'purple',
               'orange'),
       pch = c(21,21,21,21,NA),
       lty = c(NA, NA, NA, NA, 1),
       lwd = c(2,2,2,2,3)
       )
#loop through individuals
for(ii in u_id)
{
  mnv_gh = Mvec(subset(cd, ID %in% ii & colour %in% 'g' & brightn %in% 'h')$angle)
  mnv_gl = Mvec(subset(cd, ID %in% ii & colour %in% 'g' & brightn %in% 'l')$angle)
  mnv_uh = Mvec(subset(cd, ID %in% ii & colour %in% 'u' & brightn %in% 'h')$angle)
  mnv_ul = Mvec(subset(cd, ID %in% ii & colour %in% 'u' & brightn %in% 'l')$angle)
  plot(x = NULL,
       xlim = 100*c(-1,1),
       ylim = 100*c(-1,1),
       pch = 19,
       axes = FALSE,
       xlab = '',
       ylab = '',
       main = ''
       )
  abline(a = 0, b = 1, col = 'gray90')
  abline(a = 0, b = -1, col = 'gray90')
  abline(h = 0, v = 0, col = 'gray75')
  lines(x = 20*sin(xc), y = 20*cos(xc), lty = 3)
  lines(x = 40*sin(xc), y = 40*cos(xc), lty = 3)
  lines(x = 60*sin(xc), y = 60*cos(xc), lty = 3)
  lines(x = 80*sin(xc), y = 80*cos(xc), lty = 3)
  text(x = 1:4*20,
       y = c(0,0,0,0),
       labels = paste(' run', 1:4*20), 
       cex = 0.3,
       adj = c(0,1))
  with(subset(cd, ID %in% ii),
       {
  lines(x = c(0, 100*sin(sun_az_rad)),
        y = c(0, 100*cos(sun_az_rad)),
        col = adjustcolor('orange',alpha.f = 0.5),
        lwd = 2)
  points(x = run*sin(as.numeric(angle)),
         y = run*cos(as.numeric(angle)),
         bg = gray(level = 1.0,
                   alpha =  0.4),
         col = c('green', 'darkgreen', 'magenta', 'purple')
                [ifelse(colour %in% 'g',
                        yes = ifelse(brightn %in% 'h', yes = 1, no = 2),
                        no = ifelse(brightn %in% 'h', yes = 3, no = 4))],
        pch = 21,
        lwd = 2
       )
       }
  )
  lines(x = c(0,100*sin(mnv_gh['mu'])*mnv_gh['rho']), 
        y = c(0,100*cos(mnv_gh['mu'])*mnv_gh['rho']), 
        col = 'green', 
        lwd = 1)
  lines(x = c(0,100*sin(mnv_gl['mu'])*mnv_gl['rho']), 
        y = c(0,100*cos(mnv_gl['mu'])*mnv_gl['rho']), 
        col = 'darkgreen', 
        lwd = 1)
  lines(x = c(0,100*sin(mnv_uh['mu'])*mnv_uh['rho']), 
        y = c(0,100*cos(mnv_uh['mu'])*mnv_uh['rho']), 
        col = 'magenta', 
        lwd = 1)
  lines(x = c(0,100*sin(mnv_ul['mu'])*mnv_ul['rho']), 
        y = c(0,100*cos(mnv_ul['mu'])*mnv_ul['rho']), 
        col = 'purple', 
        lwd = 1)
  mtext(ii,side = 3,line = -2)
}
#save and open
dev.off()
shell.exec.OS(pdf_fl)

# Calculate mean vectors --------------------------------------------------


#calculate mean vectors
mean_vectors = aggregate(angle~ID*brightn*colour+sun_az, # N.B. Including sun azimuth drops some cases without a time stamp
                         data = cd,
                         FUN = rho.circular
)
#correct names
mean_vectors = within(mean_vectors,
                      {mean_vector = angle; rm(angle)} # anlge now indicates a mean vector, not an angle
)
#add kappa
MLE_est = function(x)
{
          with(mle.vonmises(circular(x,
                                      template = 'none'),
                             bias = TRUE),
                c(mu = as.numeric(mu), 
                  kappa = kappa))
  
}

mle_estimates = aggregate(angle~ID*brightn*colour+sun_az,
                          data = cd,
                          FUN = MLE_est
)
#add to the summary table and 
#calculate inverse softplus kappa
mean_vectors = within(mean_vectors,
                      {
                        mu = deg(mle_estimates$angle[,'mu'])
                        kappa = mle_estimates$angle[,'kappa']
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
  
invisible(
  {
    lapply(X = u_id,
           FUN = Plt_mvec
           )
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

invisible(
  {
    lapply(X = u_id,
           FUN = Plt_kappa
           )
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

invisible(
  {
    lapply(X = u_id,
           FUN = Plt_iskappa
           )
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
  

invisible(
  {
    lapply(X = u_id,
           FUN = Plt_mu
           )
  })

abline(h = c(-180, -90, 0, 90,  180),
       lty = c(1,3,1,3,1))
#plot angles

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

par(mfrow = c(2,2), mar = c(0,0,0,0))
par(pty = 's')
PCfun(angles = subset(mean_vectors, colour %in% 'g' & brightn %in% 'h')$mu,
              col = 'green')
PCfun(angles = subset(mean_vectors, colour %in% 'g' & brightn %in% 'l')$mu,
              col = 'darkgreen')
PCfun(angles = subset(mean_vectors, colour %in% 'u' & brightn %in% 'h')$mu,
              col = 'magenta')
PCfun(angles = subset(mean_vectors, colour %in% 'u' & brightn %in% 'l')$mu,
              col = 'purple')


#average orientation changes
#relative to green bright
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

mu_diff_gh = sapply(X = u_id,
                    FUN = MuDiff,
                    dt = mean_vectors,
                    cl = 'g',
                    br = 'h')
unique(mu_diff_gh)
# [[1]]
# [1] 0

mu_diff_gl = sapply(X = u_id,
                    FUN = MuDiff,
                    dt = mean_vectors,
                    cl = 'g',
                    br = 'l')
mu_diff_uh = sapply(X = u_id,
                    FUN = MuDiff,
                    dt = mean_vectors,
                    cl = 'u',
                    br = 'h')
mu_diff_ul = sapply(X = u_id,
                    FUN = MuDiff,
                    dt = mean_vectors,
                    cl = 'u',
                    br = 'l')

par(mfrow = c(2,2), mar = c(0,0,0,0))
par(pty = 's')
# PCfun(angles = subset(mean_vectors, colour %in% 'g' & brightn %in% 'h')$mu,
#       col = 'green')
PCfun(angles = unlist(mu_diff_gl),
      col = 'darkgreen',
      shrink = 2,
      title = 'Green Dim - Green Bright')
PCfun(angles = unlist(mu_diff_uh),
      col = 'magenta',
      shrink = 2,
      title = 'UV Bright - Green Bright')
PCfun(angles = unlist(mu_diff_ul),
      col = 'purple',
      shrink = 2,
      title = 'UV Dim - Green Bright')

circ_mu_diff_ul = circular(x = unlist(mu_diff_ul),
                           units = 'degrees',
                           rotation = 'clock',
                           zero = pi/2)
cml_ul = circ_mle(circ_mu_diff_ul)
plot_circMLE(data = circ_mu_diff_ul, 
             table = cml_ul)
text(labels = 'Bimodal Model of UV Dim - Green Bright',
     x = 0,
     y = -1.5)

#Other contrasts
mu_diff_uhl = sapply(X = u_id,
                    FUN = MuDiff,
                    dt = mean_vectors,
                    cl = 'u',
                    br = 'l',
                    ref_cl = 'u',
                    ref_br = 'h')
mu_diff_gul = sapply(X = u_id,
                    FUN = MuDiff,
                    dt = mean_vectors,
                    cl = 'u',
                    br = 'l',
                    ref_cl = 'u',
                    ref_br = 'h')

par(mfrow = c(2,2), mar = c(0,0,0,0))
par(pty = 's')
PCfun(angles = unlist(mu_diff_gl),
      col = 'darkgreen',
      shrink = 2,
      title = 'Green Dim - Green Bright')
PCfun(angles = unlist(mu_diff_uhl),
      col = 'purple',
      shrink = 2,
      title = 'UV Dim - UV Bright')
PCfun(angles = unlist(mu_diff_uh),
      col = 'gray40',
      shrink = 2,
      title = 'UV Bright - Green Bright')
PCfun(angles = unlist(mu_diff_gul),
      col = 'gray25',
      shrink = 2,
      title = 'UV Dim - Green Dim')


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

par(mfrow = c(2,2), mar = c(0,0,0,0))
par(pty = 's')
Plt_cmle(dt = unlist(mu_diff_gl),
      col = 'darkgreen',
      title = 'Green Dim - Green Bright')
Plt_cmle(dt = unlist(mu_diff_uhl),
      col = 'purple',
      title = 'UV Dim - UV Bright')
Plt_cmle(dt = unlist(mu_diff_uh),
      col = 'gray40',
      title = 'UV Bright - Green Bright')
Plt_cmle(dt = unlist(mu_diff_gul),
      col = 'gray25',
      title = 'UV Dim - Green Dim')




# Dances relative to sun azimuth ------------------------------------------
#relative to green bright
SunDiff = function(id, dt, cl, br)
{
  with(
    subset(dt,
           ID %in% id),
    deg(
      mod_circular(
        rad(mu[colour %in% cl & brightn %in% br] -  
              unique(sun_az) ) # sun azimuth in degrees
      )
    )
  )
}

sun_diff_gh = sapply(X = u_id,
                    FUN = SunDiff,
                    dt = mean_vectors,
                    cl = 'g',
                    br = 'h')
sun_diff_gl = sapply(X = u_id,
                    FUN = SunDiff,
                    dt = mean_vectors,
                    cl = 'g',
                    br = 'l')
sun_diff_uh = sapply(X = u_id,
                    FUN = SunDiff,
                    dt = mean_vectors,
                    cl = 'u',
                    br = 'h')
sun_diff_ul = sapply(X = u_id,
                    FUN = SunDiff,
                    dt = mean_vectors,
                    cl = 'u',
                    br = 'l')


par(mfrow = c(2,2), mar = c(0,0,0,0))
par(pty = 's')
PCfun(angles = unlist(sun_diff_gh),
      col = 'green',
      shrink = 2.5,
      title = 'Green Bright Geographic Bearing')
PCfun(angles = unlist(sun_diff_gl),
      col = 'darkgreen',
      shrink = 2.5,
      title = 'Green Dim Geographic Bearing')
PCfun(angles = unlist(sun_diff_uh),
      col = 'magenta',
      shrink = 2.5,
      title = 'UV Bright Geographic Bearing')
PCfun(angles = unlist(sun_diff_ul),
      col = 'purple',
      shrink = 2.5,
      title = 'UV Dim Geographic Bearing')


# Isolate full condition indviduals ---------------------------------------
length(mu_diff_gl) - length(unlist(mu_diff_gl))
#around 144 individuals did not complete gl and gh, resulting in nulls

#Most individuals that made it to bright green 
full_ids = unique(subset(mean_vectors,
                         brightn %in% 'h' &
                           colour %in% 'g')$ID)
length(full_ids)#52 individuals
#extract just those individuals
# full_mean_vectors = subset(mean_vectors, ID %in% full_ids)
full_mu_diff_gl = mu_diff_gl[u_id %in% full_ids]
full_mu_diff_uh = mu_diff_uh[u_id %in% full_ids]
full_mu_diff_ul = mu_diff_ul[u_id %in% full_ids]

#still some missing, actually more haphazard than that
length(full_mu_diff_gl) - length(unlist(full_mu_diff_gl))
length(full_mu_diff_uh) - length(unlist(full_mu_diff_uh))
length(full_mu_diff_ul) - length(unlist(full_mu_diff_ul))


# Estimate Coefficients for Model -----------------------------------------

