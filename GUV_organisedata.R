#FOR A 'CLEAN' RUN, PRESS ctrl+shift+F10 to RESTART Rstudio
graphics.off()
# Details ---------------------------------------------------------------
#       AUTHOR:	James Foster              DATE: 2024 08 13
#     MODIFIED:	James Foster              DATE: 2026 01 03
#
#  DESCRIPTION: Reorganise data and print out to csv
#               
#       INPUTS: 
#               
#      OUTPUTS: csv
#
#	   CHANGES: - added sun elevation
#	            - convert anticlockwise to clockwise
#	            - correct time codes with underscore
#
#   REFERENCES: Natalie Cooper & Pen-Yuan Hsing, 2017 
#               A Guide to Reproducible Code in Ecology and Evolution
#               British Ecological Society 
#               https://britishecologicalsociety.org/wp-content/uploads/2017/12/guide-to-reproducible-code.pdf
#
#    EXAMPLES:  
#
# 
#TODO   ---------------------------------------------
#TODO   
#- Load data  +
#- Reorganise +
#- Add solar azimuth +

# Useful functions --------------------------------------------------------

## Load packages ----------------------------------------------------------
#needs installing before first use (in Rstudio, see automatic message)
suppressMessages(#these are disturbing users unnecessarily
  {
    require(circular)#package for handling circular data
    require(oce)#package for calculating sun azimuth
  }
)


## General functions -----------------------------------------------------

Mod360.180 = function(x)
{#use atan2 to convert any angle to the range (-180,180)
  deg(
    atan2(y = sin(rad(x)),
          x = cos(rad(x))
    )
  )
}



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
    caption = 'Please select the ".csv" file'
  )
}else{
  message('\n\nPlease select the ".csv" file\n\n')
  Sys.sleep(0.5)#goes too fast for the user to see the message on some computers
  path_file <- file.choose(new=F)
}
#show the user the path they have selected
if(is.null(path_file) | !length(path_file))
{stop('No file selected.')}else
{print(path_file)}


# Read in the data and format ---------------------------------------------

cd = read.table(file = path_file, 
                header = T, 
                sep  = ',')
View(cd)

#measured_angle is anticlockwise, convert to clockwise for comparison with sun azimuth
cd = within(cd,
            {bearing = (-measured_angle %% 360)})
#check relationship
par(pty = 's')
with(cd, plot(measured_angle, bearing))

#format other variables
cd = within(cd,
            {
              ID = as.factor(bee) # beedance identifier as a factor
              date = as.factor(date) # date as a factor
              signed_angle = Mod360.180(bearing)  # bearing between -180 and 180
              angle = circular(rad(signed_angle),# bearing between -pi and pi
                               rotation = 'clock') # circular format suppresses later warnings
              rm(bee); rm(signed_angle);
            }
)

#remove ul2
cd = subset(cd, 
            !(light_type %in% 'ul2') )

u_id = with(cd, unique(ID)) # unique beedances
length(u_id)#169 beedances

#check which individuals experienced both conditions (not necessary w/ brms?)
idx_both  = with(cd, 
                 {
                   unique(subset(ID, light_type == 'ul')) %in% 
                     unique(subset(ID, light_type == 'uh'))
                 }
)
#extract the names
u_id_both = u_id[idx_both]


# Reorganise data ---------------------------------------------------------

## Split condition into component parts ----------------------------------
cd = within(cd,
            {
              col_bright = do.call(what = rbind,
                                   args = strsplit(x = light_type,
                                                   split = '')
              )
              colour = col_bright[,1]
              brightn = col_bright[,2]
              rm(col_bright); rm(light_type); 
            }
)


## Add time --------------------------------------------------------------

IDtocode = function(x)
{
  tp = strsplit(x = as.character(x),
                    split = '\\.')[[1]]
  tcode = paste0(tp[1], '-', tp[2], '-', tp[3], ' ',
                 tp[4], ':', 
                 sub(pattern = '\\_.*', replacement = '',
                     x = tp[5]), #some of these have an extra _1 or _2
                 ':00'
                 )
  return(tcode)
}

cd = within(cd,
            {
              time = sub(pattern = '^...........',
                         x = ID,
                         replacement = '')
              time_code = sapply(X = ID, FUN = IDtocode)
              cet_time =  as.POSIXct(time_code, tz = "Europe/Berlin", format = "%Y-%m-%d %H:%M:%OS")#CET time
              utc_time =  as.POSIXct(cet_time, tz = "UTC")#UTC time
              rm(time_code); rm(time)
            }
)

# Add sun azimuth ---------------------------------------------------------
# Breitengrad: 50.806299; Längengrad: 8.811341

#handling function (sunAngle is bad with NAs and vectors)
GetSaz = function(tm,
                  lon,
                  lat)
{
  if(!is.na(tm))
  {
    oce::sunAngle(
    t = tm,
    longitude = lon,
    latitude = lat,
    useRefraction = TRUE)$azimuth
  }else
  {
    tm
  }
}
#handling function (sunAngle is bad with NAs and vectors)
GetSel = function(tm,
                  lon,
                  lat)
{
  if(!is.na(tm))
  {
    oce::sunAngle(
    t = tm,
    longitude = lon,
    latitude = lat,
    useRefraction = TRUE)$altitude
  }else
  {
    tm
  }
}

#add sun azimuth
cd = within(cd, 
            {
            sun_az = mapply(
                   FUN = GetSaz,
                   tm = utc_time,
                   lon = 8.811341,
                   lat = 50.806299)
            }
            )

#add sun azimuth
cd = within(cd, 
            {
            sun_el = mapply(
                   FUN = GetSel,
                   tm = utc_time,
                   lon = 8.811341,
                   lat = 50.806299)
            }
            )

#convert to signed radians for modelling
cd = within(cd, 
            {
            sun_az_rad = circular(rad(Mod360.180(sun_az)),# bearing between -pi and pi
                             rotation = 'clock') # circular format suppresses later warnings
            sun_el_rad = circular(rad(Mod360.180(sun_el)),# for trignometric predictions
                             rotation = 'clock') # circular format suppresses later warnings
            rm(utc_time)
            }
            )


# Write out table ---------------------------------------------------------

write.table(x = cd,
            file = file.path(dirname(path_file), "colour_dance_reorg.csv"),
            sep = ',',
            row.names = FALSE
            )
