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
#- Read in weather data +
#- Organise weather data  +
#- Plot mean vectors by weather category  +
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
with(cd,
     {
hist(x = cloud_prop,
     breaks = 8,
     xlab = 'Cloud cover (proportion)',
     probability = TRUE,
     border = NA,
     ylim = c(-0.5,5))
boxplot(x = ccloud_prop,
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


#  Differences in mean angle ----------------------------------------------
mean_vectors_full_clear = subset(mean_vectors_full,
                                 cloud_prop < 0.26)
mean_vectors_full_cloudy = subset(mean_vectors_full,
                                 cloud_prop > 0.26)

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

par(mfrow = c(2,4), mar = c(0,0,0,0))
PCfun(angles = unlist(mu_diff_gl_clear),
      col = 'darkgreen',
      shrink = 1.5,
      title = 'Green Dim - Green Bright\nClear')
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
      title = 'Green Dim - Green Bright\nCloudy')
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


# Fit MLE models ----------------------------------------------------------
## Organise heading changes -----
#collect all heading differences
pair_diffs_lst = lapply(X = list(
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
pair_diffs_lst = lapply(X = pair_diffs_lst,
                        FUN = circular,
                        units = angle_unit,
                        rotation = angle_rot)
# Fit maximum likelihood von Mises distributions ------
#Distribution for same means
ml_same_lst = lapply(X = pair_diffs_lst,
                     mle.vonmises,
                     bias = TRUE,
                     mu = circular(x = 0,
                                   units = angle_unit,
                                   rotation = angle_rot)
)
#Distribution for different means
ml_diff_lst = lapply(X = pair_diffs_lst,
                     mle.vonmises,
                     bias = TRUE)
#Distribution for multiple different means
mod_diff_lst = lapply(X = pair_diffs_lst,
                      FUN = circ_mle,#converts circular format data to radians
                      niter = 1e4)#some stochastic variability for green stim
#Forced multiple different means
bim_diff_lst = lapply(X = pair_diffs_lst,
                      FUN = circ_mle,#converts circular format data to radians
                      exclude = c('M1','M2A', 'M2B', 'M2C'),#no uniform or unimodal
                      niter = 1e4)#some stochastic variability for green stim


#Extract the parameters of the model with the lowest AIC
mod_diff_par = t( sapply(X = mod_diff_lst,
                         FUN = MD_extract) )
bim_diff_par = t( sapply(X = bim_diff_lst,
                         FUN = MD_extract) )

# Calculate likelihood 
#calculate all in list
all_ll = mapply(FUN = LLcalc,
                ml = c(ml_same_lst, ml_diff_lst), #apply to both sets of ML estimates
                angles = pair_diffs_lst, #repeat datasets
                au = angle_unit, #important, use the same frame of reference as the original
                ar = angle_rot)

#inspect resulting parameters
mod_diff_results = do.call(rbind, c(ml_same_lst, ml_diff_lst))
#add easy to parse names
rownames(mod_diff_results) = paste( sort( #distinguish same and different mean
  rep(x = c( 'same', 'diff'), 
      times =  length(pair_diffs_lst) ),
  decreasing = TRUE
),
rep(x = names(pair_diffs_lst),
    times = 2) 
)
#display version
rspar = mod_diff_results[,c('mu', 'kappa', 'se.mu', 'se.kappa')]
rspar = cbind(rspar, loglikelihood = all_ll)
rspar = apply(X = rspar, MARGIN = 2, FUN = unlist)
#reorder for easy comparison
rspar = rspar[c(rbind(1:8, 1:8+8)),]#put the different directly after the same
rspar[,3:4] = NA
rspar = cbind(rspar, weight1 = 1.0)
colnames(rspar) = colnames(bim_diff_par)[c(1:4, 6:5)]
all_res = rbind(data.frame(rspar),
                data.frame(bim_diff_par) )
all_res = apply(all_res, MARGIN = 2, FUN = unlist)
rownames(all_res[17:24,]) = paste('multi', rownames(bim_diff_par))
#TODO reorganise this
all_res = all_res[c(1:2, 9, 3:4, 10, 5:6, 11, 7:8, 12),]

all_res = within(data = data.frame(all_res),
                 {
                   model = c('same',
                             'diff',
                             'multi')
                   dataset = c(sapply(X = names(pair_diffs_lst),
                                      FUN = rep,
                                      times = 3) )
                 }
)
# print(rspar, digits = 3)
# print(mod_diff_par, digits = 3)
# print(bim_diff_par, digits = 3)
print(all_res, digits = 3)

