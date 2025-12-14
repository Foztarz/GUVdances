#FOR A 'CLEAN' RUN, PRESS ctrl+shift+F10 to RESTART Rstudio
graphics.off()
# Details ---------------------------------------------------------------
#       AUTHOR:	James Foster              DATE: 2025 11 02
#     MODIFIED:	James Foster              DATE: 2025 12 13
#
#  DESCRIPTION: Sort data by weather based on weather reports near Marbug
#               
#       INPUTS: 
#               
#      OUTPUTS: Plots and test statistics
#
#	   CHANGES: - 
#
#   REFERENCES: 
# 
#
#    EXAMPLES:  
#
# 
#TODO   ---------------------------------------------
#TODO   
#- Read in weather data +
#- Organise weather data  +
#- Plot mean vectors by weather category  +
#- Import Neu-Ulrichstein data
#- Model comparison cloudy non-cloudy


# Set up workspace --------------------------------------------------------
angle_unit = 'degrees'
angle_rot = 'clock'

#Check the operating system and assign a logical flag (T or F)
sys_win = Sys.info()[['sysname']] == 'Windows'
#find file location
here_path = tryCatch(expr = #look in the folder containing this file: sys.frame(1)$ofile
                       {file.path(dirname(sys.frame(1)$ofile))},
                     error = function(e)
                     {#if that fails, try to find the "Documents" folder
                       file.path(Sys.getenv('HOME'))
                     }
)
path_functions = file.path(here_path,
                           'GUV_functions.R',
                           fsep = if(sys_win){'\\'}else{'/'})
path_functions = tryCatch(expr = #try to load functions from the folder containing this file
           {
             source(file = path_functions)
             return(file.path(here_path,
                              'GUV_functions.R',
                              fsep = if(sys_win){'\\'}else{'/'}
             )
             )
           },
         error = function(e)
         {#if that fails, ask the user
           
           path_functions =  
             if(sys_win){#choose.files is only available on Windows
               message('\n\nPlease select the "GUV_functions.R" file\n\n')
               Sys.sleep(0.5)#goes too fast for the user to see the message on some computers
               choose.files(
                 default = file.path(here_path, "*.R"),#For some reason this is not possible in the "root" user
                 caption = 'Please select the "GUV_functions.R" file'
               )
             }else{
               message('\n\nPlease select the "GUV_functions.R" file\n\n')
               Sys.sleep(0.5)#goes too fast for the user to see the message on some computers
               file.choose(new=F)
             }
           #show the user the path they have selected
           if(is.null(path_functions) | !length(path_functions))
           {stop('No file selected.')}else
           {print(path_functions)}
           source(path_functions)
           return(path_functions)#overwrite with user selected
         }
)

###  Find Gießen data ----------------------------------------------------------
path_file = file.path(dirname(path_functions), "1Data/colour_dance_reorg.csv")
if(file.exists(path_file))
{
  print(path_file)
}else
{
  # set path to file
  if(sys_win){#choose.files is only available on Windows
    message('\n\nPlease select the ".csv" file\n\n')
    Sys.sleep(0.5)#goes too fast for the user to see the message on some computers
    path_file = choose.files(
      default = file.path(here_path, "*.csv"),#For some reason this is not possible in the "root" user
      caption = 'Please select the "colour_dance_reorg.csv" file'
    )
  }else{
    message('\n\nPlease select the "colour_dance_reorg.csv" file\n\n')
    Sys.sleep(0.5)#goes too fast for the user to see the message on some computers
    path_file = file.choose(new=F)
  }
  #show the user the path they have selected
  if(is.null(path_file) | !length(path_file))
  {stop('No file selected.')}else
  {print(path_file)}
}

weather_file = file.path(dirname(path_file), "Gießen_Cloud_19490101_20241231_01639.txt")
if(file.exists(weather_file))
{
  print(weather_file)
}else
{
  # set path to file
  if(sys_win){#choose.files is only available on Windows
    message('\n\nPlease select the ".txt" file\n\n')
    Sys.sleep(0.5)#goes too fast for the user to see the message on some computers
    weather_file = choose.files(
      default = file.path(here_path, "*.txt"),#For some reason this is not possible in the "root" user
      caption = 'Please select the "Gießen_Cloud_19490101_20241231_01639.txt" file'
    )
  }else{
    message('\n\nPlease select the "Gießen_Cloud_19490101_20241231_01639.txt" file\n\n')
    Sys.sleep(0.5)#goes too fast for the user to see the message on some computers
    weather_file = file.choose(new=F)
  }
  #show the user the path they have selected
  if(is.null(weather_file) | !length(weather_file))
  {stop('No file selected.')}else
  {print(weather_file)}
}

###  Find Gießen data ----------------------------------------------------------

weather_nu = file.path(dirname(path_file), "Weather-Neu-Ulrichstein_wetter-online.de.csv")
if(file.exists(weather_nu))
{
  print(weather_nu)
}else
{
  # set path to file
  if(sys_win){#choose.files is only available on Windows
    message('\n\nPlease select the ".csv" file\n\n')
    Sys.sleep(0.5)#goes too fast for the user to see the message on some computers
    weather_nu = choose.files(
      default = file.path(here_path, "*.csv"),#For some reason this is not possible in the "root" user
      caption = 'Please select the "Weather-Neu-Ulrichstein_wetter-online.de.csv" file'
    )
  }else{
    message('\n\nPlease select the "Weather-Neu-Ulrichstein_wetter-online.de.txt" file\n\n')
    Sys.sleep(0.5)#goes too fast for the user to see the message on some computers
    weather_nu = file.choose(new=F)
  }
  #show the user the path they have selected
  if(is.null(weather_nu) | !length(weather_nu))
  {stop('No file selected.')}else
  {print(weather_nu)}
}

# Read in the data and format ---------------------------------------------


## Dance data ------------------------------------------------------------

#select the reorganised data
cd = read.table(file = path_file, 
                header = T, 
                sep  = ',')
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

#full condition
luc = sapply(u_id,
             FUN = IndCond,
             dt = cd)

#Most individuals that made it to bright green 
full_ids = u_id[luc == 4]
length(full_ids)#19 individuals
#extract just those individuals

## Weather data ----------------------------------------------------------

#select the reorganised data
wd = read.table(file = weather_file, 
                header = T, 
                sep  = ';')#select the reorganised data

nd = read.table(file = weather_nu, 
                header = T, 
                sep  = ',')

## Add time --------------------------------------------------------------


Time2Hour = function(x)
{
  tt = strsplit(x, split = '\\-|\\ |\\:')[[1]] #split at literal '-', ' ' and ':'
  tt = paste(tt[-(length(tt) + #remove minutes and seconds
                    c(-1,0))],
                    collapse = '')
  tt = as.numeric(tt) #convert to numeric
  return(tt)
}

daterange = with(cd,
                 {
                 range(cet_time,
                       na.rm = TRUE)  
                 }
                 )

daterange = sapply(X = daterange,
                   FUN = Time2Hour)

#choose just the relevant time period

wd = subset(wd,
            MESS_DATUM > (daterange[1] -1) &
            MESS_DATUM < (daterange[2] +1)
            )


## Add ID ----------------------------------------------------------------

nd = within(nd,
            {
              Time = sub(Time,
                         pattern = ':',
                         replacement = '.')
              ID = ifelse(is.na(Time),
                          yes = paste0(Date, '.__.__'),
                          no = paste0(Date, '.', Time)
                          )
            }
            )


# Find corresponding data -------------------------------------------------
cd = within(cd,
            {
            hour = sapply(cet_time,
                          FUN = Time2Hour)
            }
            )

wd_subs = subset(wd,
                 MESS_DATUM %in% cd$hour
                 )

GetCloud = function(hr,dt)
{
  if(hr %in% dt$MESS_DATUM)
  {
  cc = subset(dt, MESS_DATUM %in% hr)
  return(as.numeric(cc$V_N))
  }else
  {
    return(NA) #some videos are missing a date
  }
}

cd = within(cd,
            {
            cloud8 = do.call(what = rbind,
                             args = lapply(hour,
                            GetCloud,
                            dt = wd_subs
                            ))
            cloud_prop = cloud8/8
            }
)

cd = within(cd,
            {
              nu_obs = ifelse(ID %in% nd$ID,
                              yes = nd$Wetter,
                              no = NA)
            }
            )

cd = within(cd,
            {
              nu_lev = sapply(nu_obs,
                              FUN = switch,
                              `blauer Himmel` = 1,
                              `leicht wolkig` = 2,
                              `bewölkt` = 3,
                              `NA` = NA
                              )
            }
            )


#find data for full condition individuals
full_cd = subset(cd, ID %in% full_ids)

# Inspect data ------------------------------------------------------------
#The two sites report very different weather!
with(cd,
     stripchart(nu_lev~cloud_prop,
                vertical = TRUE,
                method = 'stack',
                offset = 0.015,
          xlab = 'Gießen cloud data (0-8/8)',
          ylab = 'Neu-Ulrichstein cloud data (1-3)',
          col = gray(0, alpha = 0.02),
          pch = 1)
     )


## Gießen data -----------------------------------------------------------

with(cd,
     {
hist(x = cloud_prop,
     breaks = 8,
     xlab = 'Cloud cover (proportion)',
     main = 'Cloud cover for all dances',
     probability = TRUE,
     border = NA,
     ylim = c(-0.5,5))
boxplot(x = cloud_prop,
        at = -0.25,
        horizontal = TRUE,
        col = adjustcolor('blue', alpha.f = 0.2),
        border = gray(0.25, 1.0),
        add = TRUE)
     }
)
#about half of data is below 2/8, might be a good cut

with(full_cd,
     {
hist(x = cloud_prop,
     breaks = 8,
     xlab = 'Cloud cover (proportion)',
     main = 'Cloud cover for full-condition dances',
     probability = TRUE,
     border = NA,
     ylim = c(-0.5,5))
boxplot(x = cloud_prop,
        at = -0.25,
        horizontal = TRUE,
        col = adjustcolor('blue', alpha.f = 0.2),
        border = gray(0.25, 1.0),
        add = TRUE)
     }
)
#about nearly all of data is above 2/8, half are completely overcast
with(full_cd, quantile(cloud8, na.rm = TRUE))

# cloud_threshold = 0.26
cloud_threshold = 6.5/8


## Neu-Ulrichstein data --------------------------------------------------


with(cd,
     {
       hist(x = nu_lev,
            breaks = 8,
            xlab = 'Cloud cover (1-3)',
            main = 'Cloud cover for all dances',
            probability = TRUE,
            border = NA,
            ylim = c(-0.5,5))
       boxplot(x = nu_lev,
               at = -0.25,
               horizontal = TRUE,
               col = adjustcolor('blue', alpha.f = 0.2),
               border = gray(0.25, 1.0),
               add = TRUE)
     }
)
#about half of data is below 3, might be a good cut

with(full_cd,
     {
       hist(x = nu_lev,
            breaks = 8,
            xlab = 'Cloud cover (1-3)',
            main = 'Cloud cover for full-condition dances',
            probability = TRUE,
            border = NA,
            ylim = c(-0.5,5))
       boxplot(x = nu_lev,
               at = -0.25,
               horizontal = TRUE,
               col = adjustcolor('blue', alpha.f = 0.2),
               border = gray(0.25, 1.0),
               add = TRUE)
     }
)
#about equal propotions in each
with(full_cd, quantile(nu_lev, na.rm = TRUE))


# Calculate mean vectors --------------------------------------------------


mean_vectors = aggregate(angle~ID*brightn*colour+sun_az+sun_el_rad+cloud_prop + nu_lev, # N.B. Including sun azimuth drops some cases without a time stamp
                         data = cd,
                         FUN = rho.circular
)
#correct names
mean_vectors = within(mean_vectors,
                      {mean_vector = angle; rm(angle)} # anlge now indicates a mean vector, not an angle
)
#add kappa
mle_estimates = aggregate(angle~ID*brightn*colour+sun_az+sun_el+cloud_prop + nu_lev,
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
#find the full condition individuals
mean_vectors_full = subset(mean_vectors,
                           subset = ID %in% full_ids)


# Inspect mean vectors ----------------------------------------------------


## By DWD Gießen data ----------------------------------------------------
par(pty = 's')
par(mfrow = c(2,4), mar = c(0,0,0,0))
PCfun(angles = subset(mean_vectors_full, 
                      subset = colour %in% 'g' & brightn %in% 'h'&
                        cloud_prop <cloud_threshold)$mu,
      col = 'green',
      shrink = 1.5,
      title = 'Green Bright\nClear',
      titleline = -0.5)
PCfun(angles = subset(mean_vectors_full, 
                      subset = colour %in% 'g' & brightn %in% 'l'&
                        cloud_prop <cloud_threshold)$mu,
      col = 'darkgreen',
      shrink = 1.5,
      title = 'Green Dim')
PCfun(angles = subset(mean_vectors_full, 
                      subset = colour %in% 'u' & brightn %in% 'h'&
                        cloud_prop <cloud_threshold)$mu,
      col = 'magenta',
      shrink = 1.5,
      title = 'UV Bright')
PCfun(angles = subset(mean_vectors_full, 
                      subset = colour %in% 'u' & brightn %in% 'l'&
                        cloud_prop <cloud_threshold)$mu,
      col = 'purple',
      shrink = 1.5,
      title = 'UV Dim')

PCfun(angles = subset(mean_vectors_full, 
                      subset = colour %in% 'g' & brightn %in% 'h'&
                        cloud_prop >cloud_threshold)$mu,
      col = 'green',
      shrink = 1.5,
      title = '\n\nGreen Bright\nCloudy',
      titleline = -0.5)
PCfun(angles = subset(mean_vectors_full, 
                      subset = colour %in% 'g' & brightn %in% 'l'&
                        cloud_prop >cloud_threshold)$mu,
      col = 'darkgreen',
      shrink = 1.5,
      title = 'Green Dim')
PCfun(angles = subset(mean_vectors_full, 
                      subset = colour %in% 'u' & brightn %in% 'h'&
                        cloud_prop >cloud_threshold)$mu,
      col = 'magenta',
      shrink = 1.5,
      title = 'UV Bright')
PCfun(angles = subset(mean_vectors_full, 
                      subset = colour %in% 'u' & brightn %in% 'l'&
                        cloud_prop >cloud_threshold)$mu,
      col = 'purple',
      shrink = 1.5,
      title = 'UV Dim')

## By Wetter-online Neu-Ulrichstein data ---------------------------------------
par(pty = 's')
par(mfrow = c(3,4), mar = c(0,0,0,0))
PCfun(angles = subset(mean_vectors_full, 
                      subset = colour %in% 'g' & brightn %in% 'h'&
                        nu_lev ==1)$mu,
      col = 'green',
      shrink = 1.5,
      title = 'Green Bright\nClear',
      titleline = -0.5)
PCfun(angles = subset(mean_vectors_full, 
                      subset = colour %in% 'g' & brightn %in% 'l'&
                        nu_lev ==1)$mu,
      col = 'darkgreen',
      shrink = 1.5,
      title = 'Green Dim')
PCfun(angles = subset(mean_vectors_full, 
                      subset = colour %in% 'u' & brightn %in% 'h'&
                        nu_lev ==1)$mu,
      col = 'magenta',
      shrink = 1.5,
      title = 'UV Bright')
PCfun(angles = subset(mean_vectors_full, 
                      subset = colour %in% 'u' & brightn %in% 'l'&
                        nu_lev ==1)$mu,
      col = 'purple',
      shrink = 1.5,
      title = 'UV Dim')

PCfun(angles = subset(mean_vectors_full, 
                      subset = colour %in% 'g' & brightn %in% 'h'&
                        nu_lev ==2)$mu,
      col = 'green',
      shrink = 1.5,
      title = 'Green Bright\nLight cloud',
      titleline = -0.5)
PCfun(angles = subset(mean_vectors_full, 
                      subset = colour %in% 'g' & brightn %in% 'l'&
                        nu_lev ==2)$mu,
      col = 'darkgreen',
      shrink = 1.5,
      title = 'Green Dim')
PCfun(angles = subset(mean_vectors_full, 
                      subset = colour %in% 'u' & brightn %in% 'h'&
                        nu_lev ==2)$mu,
      col = 'magenta',
      shrink = 1.5,
      title = 'UV Bright')
PCfun(angles = subset(mean_vectors_full, 
                      subset = colour %in% 'u' & brightn %in% 'l'&
                        nu_lev ==2)$mu,
      col = 'purple',
      shrink = 1.5,
      title = 'UV Dim')

PCfun(angles = subset(mean_vectors_full, 
                      subset = colour %in% 'g' & brightn %in% 'h'&
                        nu_lev ==3)$mu,
      col = 'green',
      shrink = 1.5,
      title = 'Green Bright\nOvercast',
      titleline = -0.5)
PCfun(angles = subset(mean_vectors_full, 
                      subset = colour %in% 'g' & brightn %in% 'l'&
                        nu_lev ==3)$mu,
      col = 'darkgreen',
      shrink = 1.5,
      title = 'Green Dim')
PCfun(angles = subset(mean_vectors_full, 
                      subset = colour %in% 'u' & brightn %in% 'h'&
                        nu_lev ==3)$mu,
      col = 'magenta',
      shrink = 1.5,
      title = 'UV Bright')
PCfun(angles = subset(mean_vectors_full, 
                      subset = colour %in% 'u' & brightn %in% 'l'&
                        nu_lev ==3)$mu,
      col = 'purple',
      shrink = 1.5,
      title = 'UV Dim')


#  Differences in mean angle ----------------------------------------------
#DWD Gießen
mean_vectors_full_clear = subset(mean_vectors_full,
                                 cloud_prop < cloud_threshold)
mean_vectors_full_cloudy = subset(mean_vectors_full,
                                 cloud_prop > cloud_threshold)
#Wetter-online Neu-Ulrichstein
mean_vectors_full_bh = subset(mean_vectors_full, #blauer Himmel
                                 nu_lev == 1)
mean_vectors_full_lw = subset(mean_vectors_full, #leicht wolkig
                                 nu_lev == 2)
mean_vectors_full_bw = subset(mean_vectors_full, #bewölkt
                                 nu_lev == 3)


## DWD Gießen ------------------------------------------------------------
#Clear
#Plot changes in mean dance angle
mu_diff_gl_clear = sapply(X = full_ids,
                    FUN = MuDiff,
                    dt = mean_vectors_full_clear,
                    cl = 'g',
                    br = 'l')
mu_diff_uh_clear = sapply(X = full_ids,
                    FUN = MuDiff,
                    dt = mean_vectors_full_clear,
                    cl = 'u',
                    br = 'h')
mu_diff_ul_clear = sapply(X = full_ids,
                    FUN = MuDiff,
                    dt = mean_vectors_full_clear,
                    cl = 'u',
                    br = 'l')


#Other contrasts
mu_diff_uhl_clear = sapply(X = full_ids,
                     FUN = MuDiff,
                     dt = mean_vectors_full_clear,
                     cl = 'u',
                     br = 'l',
                     ref_cl = 'u',
                     ref_br = 'h')
mu_diff_gul_clear = sapply(X = full_ids,
                     FUN = MuDiff,
                     dt = mean_vectors_full_clear,
                     cl = 'u',
                     br = 'l',
                     ref_cl = 'g',
                     ref_br = 'h')

#Cloudy
#Plot changes in mean dance angle
mu_diff_gl_cloudy = sapply(X = full_ids,
                    FUN = MuDiff,
                    dt = mean_vectors_full_cloudy,
                    cl = 'g',
                    br = 'l')
mu_diff_uh_cloudy = sapply(X = full_ids,
                    FUN = MuDiff,
                    dt = mean_vectors_full_cloudy,
                    cl = 'u',
                    br = 'h')
mu_diff_ul_cloudy = sapply(X = full_ids,
                    FUN = MuDiff,
                    dt = mean_vectors_full_cloudy,
                    cl = 'u',
                    br = 'l')


#Other contrasts
mu_diff_uhl_cloudy = sapply(X = full_ids,
                     FUN = MuDiff,
                     dt = mean_vectors_full_cloudy,
                     cl = 'u',
                     br = 'l',
                     ref_cl = 'u',
                     ref_br = 'h')
mu_diff_gul_cloudy = sapply(X = full_ids,
                     FUN = MuDiff,
                     dt = mean_vectors_full_cloudy,
                     cl = 'u',
                     br = 'l',
                     ref_cl = 'g',
                     ref_br = 'h')

par(pty = 's')
par(mfrow = c(2,4), mar = c(0,0,0,0))
PCfun(angles = unlist(mu_diff_gl_clear),
      col = 'darkgreen',
      shrink = 1.5,
      title = 'Green Dim - Green Bright\nClear',
      titleline = -0.5)
PCfun(angles = unlist(mu_diff_uhl_clear),
      col = 'purple',
      shrink = 1.5,
      title = 'UV Dim - UV Bright')
PCfun(angles = unlist(mu_diff_uh_clear),
      col = 'gray40',
      shrink = 1.5,
      title = 'UV Bright - Green Bright')
PCfun(angles = unlist(mu_diff_gul_clear),
      col = 'gray25',
      shrink = 1.5,
      title = 'UV Dim - Green Dim')
PCfun(angles = unlist(mu_diff_gl_cloudy),
      col = 'darkgreen',
      shrink = 1.5,
      title = 'Green Dim - Green Bright\nCloudy',
      titleline = -0.5)
PCfun(angles = unlist(mu_diff_uhl_cloudy),
      col = 'purple',
      shrink = 1.5,
      title = 'UV Dim - UV Bright')
PCfun(angles = unlist(mu_diff_uh_cloudy),
      col = 'gray40',
      shrink = 1.5,
      title = 'UV Bright - Green Bright')
PCfun(angles = unlist(mu_diff_gul_cloudy),
      col = 'gray25',
      shrink = 1.5,
      title = 'UV Dim - Green Dim')

## Wetter-online Neu-Ulrichstein -----------------------------------------------
#Clear
#Plot changes in mean dance angle
mu_diff_gl_bh = sapply(X = full_ids,
                    FUN = MuDiff,
                    dt = mean_vectors_full_bh,
                    cl = 'g',
                    br = 'l')
mu_diff_uh_bh = sapply(X = full_ids,
                    FUN = MuDiff,
                    dt = mean_vectors_full_bh,
                    cl = 'u',
                    br = 'h')
mu_diff_ul_bh = sapply(X = full_ids,
                    FUN = MuDiff,
                    dt = mean_vectors_full_bh,
                    cl = 'u',
                    br = 'l')


#Other contrasts
mu_diff_uhl_bh = sapply(X = full_ids,
                     FUN = MuDiff,
                     dt = mean_vectors_full_bh,
                     cl = 'u',
                     br = 'l',
                     ref_cl = 'u',
                     ref_br = 'h')
mu_diff_gul_bh = sapply(X = full_ids,
                     FUN = MuDiff,
                     dt = mean_vectors_full_bh,
                     cl = 'u',
                     br = 'l',
                     ref_cl = 'g',
                     ref_br = 'h')

#Partly Cloudy
#Plot changes in mean dance angle
mu_diff_gl_lw = sapply(X = full_ids,
                    FUN = MuDiff,
                    dt = mean_vectors_full_lw,
                    cl = 'g',
                    br = 'l')
mu_diff_uh_lw = sapply(X = full_ids,
                    FUN = MuDiff,
                    dt = mean_vectors_full_lw,
                    cl = 'u',
                    br = 'h')
mu_diff_ul_lw = sapply(X = full_ids,
                    FUN = MuDiff,
                    dt = mean_vectors_full_lw,
                    cl = 'u',
                    br = 'l')


#Other contrasts
mu_diff_uhl_lw = sapply(X = full_ids,
                     FUN = MuDiff,
                     dt = mean_vectors_full_lw,
                     cl = 'u',
                     br = 'l',
                     ref_cl = 'u',
                     ref_br = 'h')
mu_diff_gul_lw = sapply(X = full_ids,
                     FUN = MuDiff,
                     dt = mean_vectors_full_lw,
                     cl = 'u',
                     br = 'l',
                     ref_cl = 'g',
                     ref_br = 'h')

#Overcast
#Plot changes in mean dance angle
mu_diff_gl_bw = sapply(X = full_ids,
                    FUN = MuDiff,
                    dt = mean_vectors_full_bw,
                    cl = 'g',
                    br = 'l')
mu_diff_uh_bw = sapply(X = full_ids,
                    FUN = MuDiff,
                    dt = mean_vectors_full_bw,
                    cl = 'u',
                    br = 'h')
mu_diff_ul_bw = sapply(X = full_ids,
                    FUN = MuDiff,
                    dt = mean_vectors_full_bw,
                    cl = 'u',
                    br = 'l')


#Other contrasts
mu_diff_uhl_bw = sapply(X = full_ids,
                     FUN = MuDiff,
                     dt = mean_vectors_full_bw,
                     cl = 'u',
                     br = 'l',
                     ref_cl = 'u',
                     ref_br = 'h')
mu_diff_gul_bw = sapply(X = full_ids,
                     FUN = MuDiff,
                     dt = mean_vectors_full_bw,
                     cl = 'u',
                     br = 'l',
                     ref_cl = 'g',
                     ref_br = 'h')

par(pty = 's')
par(mfrow = c(3,4), mar = c(0,0,0,0))
PCfun(angles = unlist(mu_diff_gl_bh),
      col = 'darkgreen',
      shrink = 1.5,
      title = 'Green Dim - Green Bright\nClear',
      titleline = -0.5)
PCfun(angles = unlist(mu_diff_uhl_bh),
      col = 'purple',
      shrink = 1.5,
      title = 'UV Dim - UV Bright')
PCfun(angles = unlist(mu_diff_uh_bh),
      col = 'gray40',
      shrink = 1.5,
      title = 'UV Bright - Green Bright')
PCfun(angles = unlist(mu_diff_gul_bh),
      col = 'gray25',
      shrink = 1.5,
      title = 'UV Dim - Green Dim')
PCfun(angles = unlist(mu_diff_gl_lw),
      col = 'darkgreen',
      shrink = 1.5,
      title = 'Green Dim - Green Bright\nPartly Cloudy',
      titleline = -0.5)
PCfun(angles = unlist(mu_diff_uhl_lw),
      col = 'purple',
      shrink = 1.5,
      title = 'UV Dim - UV Bright')
PCfun(angles = unlist(mu_diff_uh_lw),
      col = 'gray40',
      shrink = 1.5,
      title = 'UV Bright - Green Bright')
PCfun(angles = unlist(mu_diff_gul_lw),
      col = 'gray25',
      shrink = 1.5,
      title = 'UV Dim - Green Dim')
PCfun(angles = unlist(mu_diff_gl_bw),
      col = 'darkgreen',
      shrink = 1.5,
      title = 'Green Dim - Green Bright\nOvercast',
      titleline = -0.5)
PCfun(angles = unlist(mu_diff_uhl_bw),
      col = 'purple',
      shrink = 1.5,
      title = 'UV Dim - UV Bright')
PCfun(angles = unlist(mu_diff_uh_bw),
      col = 'gray40',
      shrink = 1.5,
      title = 'UV Bright - Green Bright')
PCfun(angles = unlist(mu_diff_gul_bw),
      col = 'gray25',
      shrink = 1.5,
      title = 'UV Dim - Green Dim')


# Fit MLE models ----------------------------------------------------------


## DWD Gießen data -------------------------------------------------------

### Organise heading changes -----
#collect all heading differences
pair_diffs_lst_cc = lapply(X = list(
  #bright vs dim green
  g_hl_clear = mu_diff_gl_clear, 
  #bright vs dim UV
  uv_hl_clear = mu_diff_uhl_clear,
  #bright green vs bright uv 
  uvg_h_clear = mu_diff_uh_clear,
  #dim green vs bright UV
  uvg_l_clear = mu_diff_gul_clear,
  #bright vs dim green
  g_hl_cloud = mu_diff_gl_cloudy, 
  #bright vs dim UV
  uv_hl_cloud = mu_diff_uhl_cloudy,
  #bright green vs bright uv 
  uvg_h_cloud = mu_diff_uh_cloudy,
  #dim green vs bright UV
  uvg_l_cloud = mu_diff_gul_cloudy), 
  FUN = unlist) 
#Provide appropriate data format
angle_unit = 'degrees' 
angle_rot = 'clock'
#make sure they are in circular format
pair_diffs_lst_cc = lapply(X = pair_diffs_lst_cc,
                        FUN = circular,
                        units = angle_unit,
                        rotation = angle_rot)
### Fit maximum likelihood von Mises distributions ------
#Distribution for same means
ml_same_lst_cc = lapply(X = pair_diffs_lst_cc,
                     mle.vonmises,
                     bias = TRUE,
                     mu = circular(x = 0,
                                   units = angle_unit,
                                   rotation = angle_rot)
)
#Distribution for different means
ml_diff_lst_cc = lapply(X = pair_diffs_lst_cc,
                     mle.vonmises,
                     bias = TRUE)
#Distribution for multiple different means
mod_diff_lst_cc = lapply(X = pair_diffs_lst_cc,
                      FUN = circ_mle,#converts circular format data to radians
                      niter = 1e4)#some stochastic variability for green stim
#Forced multiple different means
bim_diff_lst_cc = lapply(X = pair_diffs_lst_cc,
                      FUN = circ_mle,#converts circular format data to radians
                      exclude = c('M1','M2A', 'M2B', 'M2C'),#no uniform or unimodal
                      niter = 1e4)#some stochastic variability for green stim


#Extract the parameters of the model with the lowest AIC
mod_diff_par_cc = t( sapply(X = mod_diff_lst_cc,
                         FUN = MD_extract) )
bim_diff_par_cc = t( sapply(X = bim_diff_lst_cc,
                         FUN = MD_extract) )

# Calculate likelihood 
#calculate all in list
all_ll_cc = mapply(FUN = LLcalc,
                ml = c(ml_same_lst_cc, ml_diff_lst_cc), #apply to both sets of ML estimates
                angles = pair_diffs_lst_cc, #repeat datasets
                au = angle_unit, #important, use the same frame of reference as the original
                ar = angle_rot)

#inspect resulting parameters
mod_diff_results_cc = do.call(rbind, c(ml_same_lst_cc, ml_diff_lst_cc))
#add easy to parse names
rownames(mod_diff_results_cc) = paste( sort( #distinguish same and different mean
  rep(x = c( 'same', 'diff'), 
      times =  length(pair_diffs_lst_cc) ),
  decreasing = TRUE
),
rep(x = names(pair_diffs_lst_cc),
    times = 2) 
)
#display version
rspar_cc = mod_diff_results_cc[,c('mu', 'kappa', 'se.mu', 'se.kappa')]
rspar_cc = cbind(rspar_cc, loglikelihood = all_ll_cc)
rspar_cc = apply(X = rspar_cc, MARGIN = 2, FUN = unlist)
#reorder for easy comparison
rspar_cc = rspar_cc[c(rbind(1:8, 1:8+8)),]#put the different directly after the same
rspar_cc[,3:4] = NA
rspar_cc = cbind(rspar_cc, weight1 = 1.0)
colnames(rspar_cc) = colnames(bim_diff_par_cc)[c(1:4, 6:5)]
all_res_cc = rbind(data.frame(rspar_cc),
                data.frame(bim_diff_par_cc) )
all_res_cc = apply(all_res_cc, MARGIN = 2, FUN = unlist)
rn_res_cc = rownames(all_res_cc)#NB in R4.5 need to make a copy of rownames
rn_res_cc[17:24] = paste('multi', rownames(bim_diff_par_cc))
rownames(all_res_cc) = rn_res_cc #
#not sure why some angles are bigger than 360°?
all_res_cc[,c(1,3)] = apply(all_res_cc[,c(1,3)], 
                            FUN = Mod360.180,
                            MARGIN = 2)
all_res_cc = all_res_cc[c(1:2, 17, 3:4, 18, 5:6, 19, 7:8, 20,
                          9:10, 21, 11:12, 22, 13:14, 23, 15:16, 24
                          ),]

all_res_cc = within(data = data.frame(all_res_cc),
                 {
                   model = c('same',
                             'diff',
                             'multi')
                   dataset = c(sapply(X = names(pair_diffs_lst_cc),
                                      FUN = rep,
                                      times = 3) )
                 }
)
# print(rspar, digits = 3)
# print(mod_diff_par, digits = 3)
# print(bim_diff_par, digits = 3)
print(all_res_cc, digits = 3)


### Likelihood ratio tests for each condition-weather set -------------------
# Perform LR tests
mod_details_cc= lapply(X = names(pair_diffs_lst_cc),
                    FUN = CollectDetails,
                    ll = all_res_cc,
                    bimod = TRUE)
names(mod_details_cc) = names(pair_diffs_lst_cc)

#set up tests for hypotheses
lr_tests_cc = c('pairs_diff_zero', #is the turn different from zero?
             'pairs_multi_diff', #is the turn in one or two directions?
             'pairs_multi_zero' #are the two directions different from zero?
)#multiple tests


#Perform collect likelihood ratio tests for each hypothesis
all_lr_results_lst_cc = lapply(X =mod_details_cc,
                               FUN = LR_calc_lst,
                               tst = lr_tests_cc)
all_lr_results_cc = do.call(what = rbind,
                            args = all_lr_results_lst_cc
)
all_lr_results_cc = within(data.frame(all_lr_results_cc),
                           {
                             comparison = c( sapply(X = names(mod_details_cc), 
                                                    FUN = rep,
                                                    times = length(lr_tests_cc)
                             ) )
                             test = rep(lr_tests_cc, 
                                        times = length(names(mod_details_cc))
                             ) 
                           }
)

#add labels to interpret each test according to p value and sign of comparison
all_lr_results_cc = within(all_lr_results_cc,
                           {
                             result = mapply(tst = lr_tests_cc,
                                             d0 = unlist(dev0),
                                             d1 = unlist(dev1),
                                             pa = p_adjusted,
                                             FUN = H1label
                             )
                           }
)
print_lr_results_cc = all_lr_results_cc
print_lr_results_cc = within(print_lr_results_cc,
                             {
                               rm(dev0,dev1, p)
                               test = sub(sub(sub(test, 
                                                  pattern = '_', 
                                                  replacement = ' '),
                                              pattern = '_', 
                                              replacement = ' '),
                                          pattern = 'pairs ',
                                          replacement = '')
                               mu1 = round(all_res_cc$mu1, 1)
                             }
)
print_lr_results_cc = within(print_lr_results_cc,
                             {mu2 = round(all_res_cc$mu2, 1)}
)
rownames(print_lr_results_cc) = with(print_lr_results_cc, paste(comparison, test))#better to keep them for the html
#Round numeric columns
# print_lr_results = apply(print_lr_results, MARGIN = 2, FUN = unlist)
num_col_cc = sapply(X = print_lr_results_cc, 
                    FUN = is.numeric)
print_lr_results_cc[,num_col_cc] = apply(X = print_lr_results_cc[,num_col_cc],
                                         MARGIN = 2, 
                                         FUN = round, 
                                         digits = 3)
#Print the results
print(print_lr_results_cc[,-6],
      digits = 3)
# View(print_lr_results_cc)

#save results
write.table(x = print_lr_results_cc, 
            file = '2Results/LRT_paired_diffs_Gweather.csv',
            sep = ',',
            row.names = TRUE)
##Conclusions:
## Under neither set of weather conditions do the green or bright stimuli
## differ significantly.
## The UV stimuli differ significantly, with a single mean for clear skies
## and two means for cloudy skies.
## The dim stimuli differ significantly for the cloudy but not the clear conditions.

## Wetter-online Neu-Ulrichstein data ------------------------------------------

### Organise heading changes -----
#collect all heading differences
pair_diffs_lst_cc_nu = lapply(X = list(
  #bright vs dim green
  g_hl_clear = mu_diff_gl_bh, 
  #bright vs dim UV
  uv_hl_clear = mu_diff_uhl_bh,
  #bright green vs bright uv 
  uvg_h_clear = mu_diff_uh_bh,
  #dim green vs bright UV
  uvg_l_clear = mu_diff_gul_bh,
  #bright vs dim green
  g_hl_part = mu_diff_gl_lw, 
  #bright vs dim UV
  uv_hl_part = mu_diff_uhl_lw,
  #bright green vs bright uv 
  uvg_h_part = mu_diff_uh_lw,
  #dim green vs bright UV
  uvg_l_part = mu_diff_gul_lw,
  #bright vs dim green
  g_hl_cloud = mu_diff_gl_bw, 
  #bright vs dim UV
  uv_hl_cloud = mu_diff_uhl_bw,
  #bright green vs bright uv 
  uvg_h_cloud = mu_diff_uh_bw,
  #dim green vs bright UV
  uvg_l_cloud = mu_diff_gul_bw
  ), 
  FUN = unlist) 
#Provide appropriate data format
angle_unit = 'degrees' 
angle_rot = 'clock'
#make sure they are in circular format
pair_diffs_lst_cc_nu = lapply(X = pair_diffs_lst_cc_nu,
                        FUN = circular,
                        units = angle_unit,
                        rotation = angle_rot)
### Fit maximum likelihood von Mises distributions ------
#Distribution for same means
ml_same_lst_cc_nu = lapply(X = pair_diffs_lst_cc_nu,
                     mle.vonmises,
                     bias = TRUE,
                     mu = circular(x = 0,
                                   units = angle_unit,
                                   rotation = angle_rot)
)
#Distribution for different means
ml_diff_lst_cc_nu = lapply(X = pair_diffs_lst_cc_nu,
                     mle.vonmises,
                     bias = TRUE)
#Distribution for multiple different means
mod_diff_lst_cc_nu = lapply(X = pair_diffs_lst_cc_nu,
                      FUN = circ_mle,#converts circular format data to radians
                      niter = 1e4)#some stochastic variability for green stim
#Forced multiple different means
bim_diff_lst_cc_nu = lapply(X = pair_diffs_lst_cc_nu,
                      FUN = circ_mle,#converts circular format data to radians
                      exclude = c('M1','M2A', 'M2B', 'M2C'),#no uniform or unimodal
                      niter = 1e4)#some stochastic variability for green stim


#Extract the parameters of the model with the lowest AIC
mod_diff_par_cc_nu = t( sapply(X = mod_diff_lst_cc_nu,
                         FUN = MD_extract) )
bim_diff_par_cc_nu = t( sapply(X = bim_diff_lst_cc_nu,
                         FUN = MD_extract) )

# Calculate likelihood 
#calculate all in list
all_ll_cc_nu = mapply(FUN = LLcalc,
                ml = c(ml_same_lst_cc_nu, ml_diff_lst_cc_nu), #apply to both sets of ML estimates
                angles = pair_diffs_lst_cc_nu, #repeat datasets
                au = angle_unit, #important, use the same frame of reference as the original
                ar = angle_rot)

#inspect resulting parameters
mod_diff_results_cc_nu = do.call(rbind, c(ml_same_lst_cc_nu, ml_diff_lst_cc_nu))
#add easy to parse names
rownames(mod_diff_results_cc_nu) = paste( sort( #distinguish same and different mean
  rep(x = c( 'same', 'diff'), 
      times =  length(pair_diffs_lst_cc_nu) ),
  decreasing = TRUE
),
rep(x = names(pair_diffs_lst_cc_nu),
    times = 2) 
)
#display version
rspar_cc_nu = mod_diff_results_cc_nu[,c('mu', 'kappa', 'se.mu', 'se.kappa')]
rspar_cc_nu = cbind(rspar_cc_nu, loglikelihood = all_ll_cc_nu)
rspar_cc_nu = apply(X = rspar_cc_nu, MARGIN = 2, FUN = unlist)
#reorder for easy comparison
rspar_cc_nu = rspar_cc_nu[c(rbind(1:12, 1:12+12)),]#put the different directly after the same
rspar_cc_nu[,3:4] = NA
rspar_cc_nu = cbind(rspar_cc_nu, weight1 = 1.0)
colnames(rspar_cc_nu) = colnames(bim_diff_par_cc_nu)[c(1:4, 6:5)]
all_res_cc_nu = rbind(data.frame(rspar_cc_nu),
                data.frame(bim_diff_par_cc_nu) )
all_res_cc_nu = apply(all_res_cc_nu, MARGIN = 2, FUN = unlist)
rn_res_cc_nu = rownames(all_res_cc_nu)#NB in R4.5 need to make a copy of rownames
rn_res_cc_nu[25:36] = paste('multi', rownames(bim_diff_par_cc_nu))
rownames(all_res_cc_nu) = rn_res_cc_nu #
#Easier range to interpret
all_res_cc_nu[,c(1,3)] = apply(all_res_cc_nu[,c(1,3)], 
                            FUN = Mod360.180,
                            MARGIN = 2)
all_res_cc_nu = all_res_cc_nu[c(1:2, 25, 3:4, 26, 5:6, 27, 7:8, 28,
                          9:10, 29, 11:12, 30, 13:14, 31, 15:16, 32,
                          17:18, 33, 19:20, 34, 21:22, 35, 23:24, 36
                          ),]

all_res_cc_nu = within(data = data.frame(all_res_cc_nu),
                 {
                   model = c('same',
                             'diff',
                             'multi')
                   dataset = c(sapply(X = names(pair_diffs_lst_cc_nu),
                                      FUN = rep,
                                      times = 3) )
                 }
)
# print(rspar, digits = 3)
# print(mod_diff_par, digits = 3)
# print(bim_diff_par, digits = 3)
print(all_res_cc_nu, digits = 3)


### Likelihood ratio tests for each condition-weather set -------------------
# Perform LR tests
mod_details_cc_nu= lapply(X = names(pair_diffs_lst_cc_nu),
                    FUN = CollectDetails,
                    ll = all_res_cc_nu,
                    bimod = TRUE)
names(mod_details_cc_nu) = names(pair_diffs_lst_cc_nu)

#set up tests for hypotheses
lr_tests_cc_nu = c('pairs_diff_zero', #is the turn different from zero?
             'pairs_multi_diff', #is the turn in one or two directions?
             'pairs_multi_zero' #are the two directions different from zero?
)#multiple tests

#Perform collect likelihood ratio tests for each hypothesis
all_lr_results_lst_cc_nu = lapply(X =mod_details_cc_nu,
                            FUN = LR_calc_lst,
                            tst = lr_tests_cc_nu)
all_lr_results_cc_nu = do.call(what = rbind,
                         args = all_lr_results_lst_cc_nu
)
all_lr_results_cc_nu = within(data.frame(all_lr_results_cc_nu),
                        {
                          comparison = c( sapply(X = names(mod_details_cc_nu), 
                                                 FUN = rep,
                                                 times = length(lr_tests_cc_nu)
                          ) )
                          test = rep(lr_tests_cc_nu, 
                                     times = length(names(mod_details_cc_nu))
                          ) 
                        }
)

#add labels to interpret each test according to p value and sign of comparison
all_lr_results_cc_nu = within(all_lr_results_cc_nu,
                        {
                          result = mapply(tst = lr_tests_cc_nu,
                                          d0 = unlist(dev0),
                                          d1 = unlist(dev1),
                                          pa = p_adjusted,
                                          FUN = H1label
                          )
                        }
)
print_lr_results_cc_nu = all_lr_results_cc_nu
print_lr_results_cc_nu = within(print_lr_results_cc_nu,
                          {
                            rm(dev0,dev1, p)
                            test = sub(sub(sub(test, 
                                               pattern = '_', 
                                               replacement = ' '),
                                           pattern = '_', 
                                           replacement = ' '),
                                       pattern = 'pairs ',
                                       replacement = '')
                            mu1 = round(all_res_cc_nu$mu1, 1)
                          }
)
print_lr_results_cc_nu = within(print_lr_results_cc_nu,
                          {mu2 = round(all_res_cc_nu$mu2, 1)}
)
rownames(print_lr_results_cc_nu) = with(print_lr_results_cc_nu, paste(comparison, test))#better to keep them for the html
#Round numeric columns
# print_lr_results = apply(print_lr_results, MARGIN = 2, FUN = unlist)
num_col_cc_nu = sapply(X = print_lr_results_cc_nu, 
                 FUN = is.numeric)
print_lr_results_cc_nu[,num_col_cc_nu] = apply(X = print_lr_results_cc_nu[,num_col_cc_nu],
                                   MARGIN = 2, 
                                   FUN = round, 
                                   digits = 3)
#Print the results
print(print_lr_results_cc_nu[,-6],
      digits = 3)
# View(print_lr_results_cc)

#save results
write.table(x = print_lr_results_cc_nu, 
            file = '2Results/LRT_paired_diffs_NUweather.csv',
            sep = ',',
            row.names = TRUE)
##Conclusions:
## Under neither set of weather conditions do the green or bright stimuli
## differ significantly.
## The UV stimuli differ significantly, with a single mean for clear skies
## and two means for cloudy skies.
## The dim stimuli differ significantly for the cloudy but not the clear conditions.

# Compare with null models excluding weather ------------------------------

#Null 
all_res_null = read.table(file = file.path(dirname(path_functions), '/2Results',
                                    'MLE_paired_diffs.csv'),
                          sep = ',',
                          header = TRUE)

mod_details_null= lapply(X = with(all_res_null, unique(dataset)),
                          FUN = CollectDetails,
                          ll = all_res_null,
                          bimod = TRUE)
names(mod_details_null) = with(all_res_null, unique(dataset))

print(all_res_null, digits = 3)
# View(print_lr_results_null)
#Select likelihoods for models:
# Ghl same
# Uhl multi
# GUh same
# GUl multi


rn_res_null = rownames(all_res_null)
rn_res_cc = rownames(all_res_cc)
rn_res_cc_nu = rownames(all_res_cc_nu)

best_names_null = c('same green_hl',
                    'multi uv_hl',
                    'same uvg_h',
                    'multi uvg_l')

best_names_cc = c("same g_hl_clear", 
                  "same g_hl_cloud",
                  "diff uv_hl_clear",
                  "multi uv_hl_cloud",
                  "same uvg_h_clear",
                  "same uvg_h_cloud",
                  "multi uvg_l_clear",
                  "multi uvg_l_cloud"
)

best_names_cc_nu = c("same g_hl_clear", 
                  "same g_hl_part",
                  "multi g_hl_cloud",
                  "multi uv_hl_clear",
                  "same uv_hl_part",
                  "multi uv_hl_cloud",
                  "diff uvg_h_clear",
                  "same uvg_h_part",
                  "same uvg_h_cloud",
                  "same uvg_l_clear",
                  "diff uvg_l_part",
                  "multi uvg_l_cloud"
)

mod_details_dt_null= do.call(rbind,mod_details_null)
rownames(mod_details_dt_null) = rn_res_null

mod_details_dt_cc= do.call(rbind,mod_details_cc)
rownames(mod_details_dt_cc) = rn_res_cc

mod_details_dt_cc_nu= do.call(rbind,mod_details_cc_nu)
rownames(mod_details_dt_cc_nu) = rn_res_cc_nu

best_null = mod_details_dt_null[best_names_null,]
best_cc = mod_details_dt_cc[best_names_cc,]
best_cc_nu = mod_details_dt_cc_nu[best_names_cc_nu,]


# Calculate joint likelihood for null models
# Calculate joint likelihood for null models
joint_null = apply(best_null[,c('ll', 'deviance', 'df')],
                   FUN = sum,
                   MARGIN = 2)

## DWD Gießen ------------------------------------------------------------

joint_cc = apply(best_cc[,c('ll', 'deviance', 'df')],
                   FUN = sum,
                   MARGIN = 2)

delta_weather = joint_cc - joint_null
lrt_joint_weather = with(data.frame(t(delta_weather)), 
                 data.frame(
                   dev_null = joint_null['deviance'],
                   dev_cloud = joint_cc['deviance'], #within trial obs. share a mean, expect lower deviance with more params
                   d.f. = df, #grand mean has half the number of params
                   chi_squared = deviance,
                   p = pchisq(q= abs(deviance), df = df, lower.tail = FALSE)
                 )
)


#Calculate joint likelihood and parameters (probabilities multiplied = log-likelihoods summed)
ll_g_hl_cc = apply(best_cc[c("same g_hl_clear", 
                              "same g_hl_cloud"),
                            c('ll', 'deviance', 'df')],
                   FUN = sum,
                   MARGIN = 2)
ll_u_hl_cc = apply(best_cc[c("diff uv_hl_clear",
                              "multi uv_hl_cloud"),
                           c('ll', 'deviance','df' )],
                 FUN = sum,
                 MARGIN = 2)
ll_uvg_h_cc = apply(best_cc[c("same uvg_h_clear",
                               "same uvg_h_cloud"),
                            c('ll','deviance', 'df' )],
                  FUN = sum,
                  MARGIN = 2)
ll_uvg_l_cc = apply(best_cc[c("multi uvg_l_clear",
                                 "multi uvg_l_cloud"),
                            c('ll','deviance', 'df' )],
                  FUN = sum,
                  MARGIN = 2)

#Likelihood ratio tests
delta_g_hl_cc = ll_g_hl_cc - 
                  best_null[best_names_null[1],c('ll','deviance', 'df' )]
delta_u_hl_cc = ll_u_hl_cc - 
                  best_null[best_names_null[2],c('ll','deviance', 'df' )]
delta_uvg_h_cc = ll_uvg_h_cc - 
                  best_null[best_names_null[3],c('ll','deviance', 'df' )]
delta_uvg_l_cc = ll_uvg_l_cc - 
                  best_null[best_names_null[4],c('ll','deviance', 'df' )]

lrt_g_hl = with(delta_g_hl_cc, 
                data.frame(
                  dev_null = best_null[best_names_null[1],'deviance'],
                  dev_cloud = ll_g_hl_cc['deviance'], #within trial obs. share a mean, expect lower deviance with more params
                  d.f. = df, #grand mean has half the number of params
                  chi_squared = deviance,
                  p = pchisq(q= abs(deviance), df = df, lower.tail = FALSE)
                )
              )
    
lrt_u_hl = with(delta_u_hl_cc, 
                data.frame(
                  dev_null = best_null[best_names_null[2],'deviance'],
                  dev_cloud = ll_u_hl_cc['deviance'], #within trial obs. share a mean, expect lower deviance with more params
                  d.f. = df, #grand mean has half the number of params
                  chi_squared = deviance,
                  p = pchisq(q= abs(deviance), df = df, lower.tail = FALSE)
                )
              )    
lrt_uvg_h = with(delta_uvg_h_cc, 
                data.frame(
                  dev_null = best_null[best_names_null[3],'deviance'],
                  dev_cloud = ll_uvg_h_cc['deviance'], #within trial obs. share a mean, expect lower deviance with more params
                  d.f. = df, #grand mean has half the number of params
                  chi_squared = deviance,
                  p = pchisq(q= abs(deviance), df = df, lower.tail = FALSE)
                )
              )    
lrt_uvg_l = with(delta_uvg_l_cc, 
                data.frame(
                  dev_null = best_null[best_names_null[4],'deviance'],
                  dev_cloud = ll_uvg_l_cc['deviance'], #within trial obs. share a mean, expect lower deviance with more params
                  d.f. = df, #grand mean has half the number of params
                  chi_squared = deviance,
                  p = pchisq(q= abs(deviance), df = df, lower.tail = FALSE)
                )
              )
lr_weather = rbind(lrt_g_hl, lrt_u_hl, lrt_uvg_h, lrt_uvg_l)
rownames(lr_weather) = best_names_null
lr_weather = within(lr_weather,
                    {
                      p_adjusted = p.adjust(p = p,
                                            method = 'BH',
                                            n = length(p))
                    })    
lr_weather = rbind(within(lrt_joint_weather, {p_adjusted = NA}), lr_weather)
rn_lr_weather = rownames(lr_weather)
rn_lr_weather[1] = 'joint models'
rownames(lr_weather) = rn_lr_weather
print_lr_weather = apply(X = lr_weather,
                         MARGIN = 2, 
                         FUN = round, 
                         digits = 3)
#Print the results
print(print_lr_weather,
      digits = 3)
#In the case of comparing the two bright stimuli
#the split by weather helps explain the small leftward turn.
#For the other condition, the models that do not account for weather
#better explain the data.
#save results
write.table(x = print_lr_results_cc, 
            file = file.path(dirname(path_functions),
                             '2Results/LRT_joint_Gweather.csv'),
            sep = ',',
            row.names = TRUE)

## Wetter-online NU ------------------------------------------------------------


# Calculate joint likelihood for weather models
joint_cc_nu = apply(best_cc_nu[,c('ll', 'deviance', 'df')],
                   FUN = sum,
                   MARGIN = 2)

delta_weather_nu = joint_cc_nu - joint_null
lrt_joint_weather_nu = with(data.frame(t(delta_weather_nu)), 
                 data.frame(
                   dev_null = joint_null['deviance'],
                   dev_cloud = joint_cc_nu['deviance'], #within trial obs. share a mean, expect lower deviance with more params
                   d.f. = df, #grand mean has half the number of params
                   chi_squared = deviance,
                   p = pchisq(q= abs(deviance), df = df, lower.tail = FALSE)
                 )
)




#Calculate joint likelihood and parameters (probabilities multiplied = log-likelihoods summed)
ll_g_hl_cc_nu = apply(best_cc_nu[c("same g_hl_clear", 
                                "same g_hl_part", 
                              "multi g_hl_cloud"),
                            c('ll', 'deviance', 'df')],
                   FUN = sum,
                   MARGIN = 2)
ll_u_hl_cc_nu = apply(best_cc_nu[c("multi uv_hl_clear",
                              "same uv_hl_part",
                              "multi uv_hl_cloud"),
                           c('ll', 'deviance','df' )],
                 FUN = sum,
                 MARGIN = 2)
ll_uvg_h_cc_nu = apply(best_cc_nu[c("diff uvg_h_clear",
                               "same uvg_h_part",
                               "same uvg_h_cloud"),
                            c('ll','deviance', 'df' )],
                  FUN = sum,
                  MARGIN = 2)
ll_uvg_l_cc_nu = apply(best_cc_nu[c("same uvg_l_clear",
                                 "diff uvg_l_cloud",
                                 "multi uvg_l_cloud"),
                            c('ll','deviance', 'df' )],
                  FUN = sum,
                  MARGIN = 2)

#Likelihood ratio tests
delta_g_hl_cc_nu = ll_g_hl_cc_nu - 
                  best_null[best_names_null[1],c('ll','deviance', 'df' )]
delta_u_hl_cc_nu = ll_u_hl_cc_nu - 
                  best_null[best_names_null[2],c('ll','deviance', 'df' )]
delta_uvg_h_cc_nu = ll_uvg_h_cc_nu - 
                  best_null[best_names_null[3],c('ll','deviance', 'df' )]
delta_uvg_l_cc_nu = ll_uvg_l_cc_nu - 
                  best_null[best_names_null[4],c('ll','deviance', 'df' )]

lrt_g_hl_nu = with(delta_g_hl_cc_nu, 
                data.frame(
                  dev_null = best_null[best_names_null[1],'deviance'],
                  dev_cloud = ll_g_hl_cc_nu['deviance'], #within trial obs. share a mean, expect lower deviance with more params
                  d.f. = df, #grand mean has half the number of params
                  chi_squared = deviance,
                  p = pchisq(q= abs(deviance), df = df, lower.tail = FALSE)
                )
              )
    
lrt_u_hl_nu = with(delta_u_hl_cc_nu, 
                data.frame(
                  dev_null = best_null[best_names_null[2],'deviance'],
                  dev_cloud = ll_u_hl_cc_nu['deviance'], #within trial obs. share a mean, expect lower deviance with more params
                  d.f. = df, #grand mean has half the number of params
                  chi_squared = deviance,
                  p = pchisq(q= abs(deviance), df = df, lower.tail = FALSE)
                )
              )    
lrt_uvg_h_nu = with(delta_uvg_h_cc_nu, 
                data.frame(
                  dev_null = best_null[best_names_null[3],'deviance'],
                  dev_cloud = ll_uvg_h_cc_nu['deviance'], #within trial obs. share a mean, expect lower deviance with more params
                  d.f. = df, #grand mean has half the number of params
                  chi_squared = deviance,
                  p = pchisq(q= abs(deviance), df = df, lower.tail = FALSE)
                )
              )    
lrt_uvg_l_nu = with(delta_uvg_l_cc_nu, 
                data.frame(
                  dev_null = best_null[best_names_null[4],'deviance'],
                  dev_cloud = ll_uvg_l_cc_nu['deviance'], #within trial obs. share a mean, expect lower deviance with more params
                  d.f. = df, #grand mean has half the number of params
                  chi_squared = deviance,
                  p = pchisq(q= abs(deviance), df = df, lower.tail = FALSE)
                )
              )
lr_weather_nu = rbind(lrt_g_hl_nu, lrt_u_hl_nu, lrt_uvg_h_nu, lrt_uvg_l_nu)
rownames(lr_weather_nu) = best_names_null
lr_weather_nu = within(lr_weather_nu,
                    {
                      p_adjusted = p.adjust(p = p,
                                            method = 'BH',
                                            n = length(p))
                    })    
lr_weather_nu = rbind(within(lrt_joint_weather_nu, {p_adjusted = NA}), lr_weather_nu)
rn_lr_weather_nu = rownames(lr_weather_nu)
rn_lr_weather_nu[1] = 'joint models'
rownames(lr_weather_nu) = rn_lr_weather_nu
print_lr_weather_nu = apply(X = lr_weather_nu,
                         MARGIN = 2, 
                         FUN = round, 
                         digits = 3)
#Print the results
print(print_lr_weather_nu,
      digits = 3)
#In the case of comparing the two bright stimuli
#the split by weather helps explain the small leftward turn.
#For the other condition, the models that do not account for weather
#better explain the data.
#save results
write.table(x = print_lr_results_cc_nu, 
            file = file.path(dirname(path_functions),
                             '2Results/LRT_joint_NUweather.csv'),
            sep = ',',
            row.names = TRUE)
