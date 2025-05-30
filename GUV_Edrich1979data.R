#FOR A 'CLEAN' RUN, PRESS ctrl+shift+F10 to RESTART Rstudio
graphics.off()
# Details ---------------------------------------------------------------
#       AUTHOR:	James Foster              DATE: 2025 05 30
#     MODIFIED:	James Foster              DATE: 2025 05 30
#
#  DESCRIPTION: Load data and plot
#               
#       INPUTS: 
#               
#      OUTPUTS: csv
#
#	   CHANGES: - 
#
#   REFERENCES: Edrich, W., Neumeyer, C. and von Helversen, O. (1979).
#               “Anti-sun orientation” of bees with regard to a field of ultraviolet light. 
#               J. Comp. Physiol. 134, 151–157.
#
#
#    EXAMPLES:  
#
# 
#TODO   ---------------------------------------------
#TODO   
#- Load data  
#- Plot
#- Fit curves




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
  message('\n\nPlease select the "Edrich1979_brightness-accuracy.csv" file\n\n')
  Sys.sleep(0.5)#goes too fast for the user to see the message on some computers
  path_file  <- choose.files(
    default = file.path(ltp,'Documents', "*.csv"),#For some reason this is not possible in the "root" user
    caption = 'Please select the "Edrich1979_brightness-accuracy.csv" file'
  )
}else{
  message('\n\nPlease select the ".Edrich1979_brightness-accuracy.csv" file\n\n')
  Sys.sleep(0.5)#goes too fast for the user to see the message on some computers
  path_file <- file.choose(new=F)
}
#show the user the path they have selected
if(is.null(path_file) | !length(path_file))
{stop('No file selected.')}else
{print(path_file)}


# Read in the data and format ---------------------------------------------

ed = read.table(file = path_file, 
                header = T, 
                sep  = ',')
View(ed)

with(ed,
     {
     plot(x = log10(intensity),
          y = accuracy,
          xlim = c(9,14),
          ylim = c(0,1),
          pch = 21,
          col = 'black',
          bg = sapply(paste(wavelength),
                       FUN = switch,
                       `354` = 'purple',
                       `429` = 'blue',
                       `535` = 'darkgreen',
                       'red')
          )
     abline(h = c(0,1))
     }
     )

