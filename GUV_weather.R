#FOR A 'CLEAN' RUN, PRESS ctrl+shift+F10 to RESTART Rstudio
graphics.off()
# Details ---------------------------------------------------------------
#       AUTHOR:	James Foster              DATE: 2025 11 02
#     MODIFIED:	James Foster              DATE: 2025 11 02
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
#- Read in weather data
#- Organise weather data
#- Plot mean vectors by weather category


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

tryCatch(expr = #try to load functions from the folder containing this file
           {
             source(file = file.path(here_path,
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
         }
)

path_file = file.path(here_path, "1Data/colour_dance_reorg.csv")
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

weather_file = file.path(here_path, "1Data/Gießen_Cloud_19490101_20241231_01639.txt")
if(file.exists(weather_file))
{
  print(weather_file)
}else
{
  # set path to file
  if(sys_win){#choose.files is only available on Windows
    message('\n\nPlease select the ".csv" file\n\n')
    Sys.sleep(0.5)#goes too fast for the user to see the message on some computers
    path_file = choose.files(
      default = file.path(here_path, "*.csv"),#For some reason this is not possible in the "root" user
      caption = 'Please select the "Gießen_Cloud_19490101_20241231_01639.txt" file'
    )
  }else{
    message('\n\nPlease select the "Gießen_Cloud_19490101_20241231_01639.txt" file\n\n')
    Sys.sleep(0.5)#goes too fast for the user to see the message on some computers
    path_file = file.choose(new=F)
  }
  #show the user the path they have selected
  if(is.null(weather_file) | !length(weather_file))
  {stop('No file selected.')}else
  {print(weather_file)}
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


## Weather data ----------------------------------------------------------

#select the reorganised data
wd = read.table(file = weather_file, 
                header = T, 
                sep  = ';')

## Add time --------------------------------------------------------------

# IDtocode = function(x)
# {
#   tp = strsplit(x = as.character(x),
#                 split = '\\.')[[1]]
#   tcode = paste0(tp[1], '-', tp[2], '-', tp[3], ' ',
#                  tp[4], ':', tp[5], ':00'
#   )
#   return(tcode)
# }
# 
# cd = within(cd,
#             {
#               time = sub(pattern = '^...........',
#                          x = ID,
#                          replacement = '')
#               time_code = sapply(X = ID, FUN = IDtocode)
#               cet_time =  as.POSIXct(time_code, tz = "Europe/Berlin", format = "%Y-%m-%d %H:%M:%OS")#CET time
#               utc_time =  as.POSIXct(cet_time, tz = "UTC")#UTC time
#               rm(time_code); rm(time)
#             }
# )

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


# Inspect data ------------------------------------------------------------

hist(x = cd$cloud_prop,
     breaks = 8,
     xlab = 'Cloud cover (proportion)',
     probability = TRUE,
     border = NA)
boxplot(x = cd$cloud_prop,
        horizontal = TRUE,
        col = adjustcolor('blue', alpha.f = 0.2),
        border = gray(0.25, 1.0),
        add = TRUE)
#about half of data is below 2/8, might be a good cut


# Calculate mean vectors --------------------------------------------------
#full condition
luc = sapply(u_id,
             FUN = IndCond,
             dt = cd)

#Most individuals that made it to bright green 
full_ids = u_id[luc == 4]
length(full_ids)#19 individuals
#extract just those individuals

#find data for full condition individuals
full_cd = subset(cd, ID %in% full_ids)

mean_vectors = aggregate(angle~ID*brightn*colour+sun_az+sun_el_rad+cloud_prop, # N.B. Including sun azimuth drops some cases without a time stamp
                         data = cd,
                         FUN = rho.circular
)
#correct names
mean_vectors = within(mean_vectors,
                      {mean_vector = angle; rm(angle)} # anlge now indicates a mean vector, not an angle
)
#add kappa
mle_estimates = aggregate(angle~ID*brightn*colour+sun_az+sun_el+cloud_prop,
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

par(mfrow = c(2,4), mar = c(0,0,0,0))
PCfun(angles = subset(mean_vectors_full, 
                      subset = colour %in% 'g' & brightn %in% 'h'&
                        cloud_prop <0.26)$mu,
      col = 'green',
      shrink = 1.5,
      title = 'Green Bright\nClear')
PCfun(angles = subset(mean_vectors_full, 
                      subset = colour %in% 'g' & brightn %in% 'l'&
                        cloud_prop <0.26)$mu,
      col = 'darkgreen',
      shrink = 1.5,
      title = 'Green Dim')
PCfun(angles = subset(mean_vectors_full, 
                      subset = colour %in% 'u' & brightn %in% 'h'&
                        cloud_prop <0.26)$mu,
      col = 'magenta',
      shrink = 1.5,
      title = 'UV Bright')
PCfun(angles = subset(mean_vectors_full, 
                      subset = colour %in% 'u' & brightn %in% 'l'&
                        cloud_prop <0.26)$mu,
      col = 'purple',
      shrink = 1.5,
      title = 'UV Dim')

PCfun(angles = subset(mean_vectors_full, 
                      subset = colour %in% 'g' & brightn %in% 'h'&
                        cloud_prop >0.26)$mu,
      col = 'green',
      shrink = 1.5,
      title = 'Green Bright\nCloudy')
PCfun(angles = subset(mean_vectors_full, 
                      subset = colour %in% 'g' & brightn %in% 'l'&
                        cloud_prop >0.26)$mu,
      col = 'darkgreen',
      shrink = 1.5,
      title = 'Green Dim')
PCfun(angles = subset(mean_vectors_full, 
                      subset = colour %in% 'u' & brightn %in% 'h'&
                        cloud_prop >0.26)$mu,
      col = 'magenta',
      shrink = 1.5,
      title = 'UV Bright')
PCfun(angles = subset(mean_vectors_full, 
                      subset = colour %in% 'u' & brightn %in% 'l'&
                        cloud_prop >0.26)$mu,
      col = 'purple',
      shrink = 1.5,
      title = 'UV Dim')
