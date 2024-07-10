#FOR A 'CLEAN' RUN, PRESS ctrl+shift+F10 to RESTART Rstudio
graphics.off()
# Details ---------------------------------------------------------------
#       AUTHOR:	James Foster              DATE: 2024 07 10
#     MODIFIED:	James Foster              DATE: 2024 07 10
#
#  DESCRIPTION: Load dance angles, calculate mean vectors and fit beta regression
#               to the mean vectors via brms.
#               
#       INPUTS: 
#               
#      OUTPUTS: Plots and test statistics
#
#	   CHANGES: - 
#   REFERENCES: Bürkner, P.-C., 2018. 
#               Advanced Bayesian Multilevel Modeling with the R Package brms.
#               The R Journal 10, 395–411. 
#               https://doi.org/10.32614/RJ-2018-017
# 
#               
#
#    EXAMPLES:  
#
# 
#TODO   ---------------------------------------------
#TODO   
#- Load data  +
#- Random intercepts model
#- Run model
#- Extract results  
#- Get modelling consistent 
#- Eliminate divergent transitions
#- Handwrite hypothesis tests for all contrasts
# Useful functions --------------------------------------------------------

# . Load packages ----------------------------------------------------------
#needs installing before first use (in Rstudio, see automatic message)
suppressMessages(#these are disturbing users unnecessarily
  {
    require(circular)#package for handling circular data
    require(cmdstanr)#package for Bayesian modelling via Stan
    require(brms)#package for preparing Stan models
  }
)


# . General functions -----------------------------------------------------

Mod360.180 = function(x)
{#use atan2 to convert any angle to the range (-180,180)
  deg(
    atan2(y = sin(rad(x)),
          x = cos(rad(x))
    )
  )
}




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

cd = within(cd,
            {
              ID = as.factor(bee) # beedance identifier as a factor
              date = as.factor(date) # date as a factor
              signed_angle = Mod360.180(bearing)  # bearing between -180 and 180
              angle = circular(rad(signed_angle),# bearing between -pi and pi
                               rotation = 'clock') # circular format suppresses later warnings
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
cd = within(cd,
            {
              col_bright = do.call(what = rbind,
                                   args = strsplit(x = light_type,
                                                   split = '')
              )
              colour = col_bright[,1]
              brightn = col_bright[,2]
            }
)

# Make into standata ------------------------------------------------------
form_v1 = bf( angle ~ light_type + (light_type |bee),
              kappa ~ light_type + (light_type |bee),
              family = von_mises)			

stan_data  = make_standata(
                formula = form_v1,
                data = cd
              )

br_code = make_stancode(
  formula = form_v1,
  data = cd
)	
#brms can't handle random effects in von_mises distributions (& others)
#save this stancode and rewrite it in the style of "handwritten.uvMEmises6-1.stan"
write.table(br_code, 
            file = file.path(dirname(path_file), 
                             'BRMS_vonmises_v1.stan'),
            col.names = F, row.names = F, quote = F)					


# Fit model ---------------------------------------------------------------
mod_v1 = cmdstanr::cmdstan_model(stan_file = 
                                   file.path( dirname(path_file), 
                                             'ME_vonmises_v1.stan')#I wrote this myself
                                 )
fit_v1 = mod_v1$sample(data = stan_data,
                       chains = 4,
                       iter_warmup = 2000,
                       iter_sampling = 2000,
                       adapt_delta = 0.95
                      )

