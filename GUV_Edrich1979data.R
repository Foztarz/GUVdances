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
          cex = 1.5,
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


# Fit models --------------------------------------------------------------
#prepare for model fitting
ed = within(ed,
            {
              logit_accuracy = qlogis(accuracy)
              log10_intensity = log10(intensity)
            }
            )


# Subset data by wavelength -----------------------------------------------

ed_uv = subset(ed, wavelength <400)
ed_bl = subset(ed, wavelength == 429)
ed_gr = subset(ed, wavelength >500)



#fit simple GLM
lm_uv = with(ed_uv,
             lm(logit_accuracy ~ log10_intensity) 
             )

lm_bl = with(ed_bl,
             lm(logit_accuracy ~ log10_intensity) 
             )

lm_gr = with(ed_gr,
             lm(logit_accuracy ~ log10_intensity) 
             )

xl = seq(from = 9,
         to = 14,
         length.out = 1e3)

nd = data.frame(log10_intensity = xl)
nd = merge(x = ed, y = nd, all.y = TRUE)

with(ed,
     {
       plot(x = log10_intensity,
            y = logit_accuracy,
            xlim = c(9,14),
            ylim = c(-5,5),
            pch = 21,
            cex = 1.5,
            col = 'black',
            bg = sapply(paste(wavelength),
                        FUN = switch,
                        `354` = 'purple',
                        `429` = 'blue',
                        `535` = 'darkgreen',
                        'red')
       )
       abline(h = c(0,1))
       abline(lm_uv,
              col = 'magenta',
              lwd = 3)
       abline(lm_bl,
              col = 'cyan',
              lwd = 3)
       abline(lm_gr,
              col = 'green',
              lwd = 3)
     }
)

with(ed,
     {
       plot(x = log10_intensity,
            y = accuracy,
            xlim = c(9,14),
            ylim = c(0,1),
            pch = 21,
            cex = 1.5,
            col = 'black',
            bg = sapply(paste(wavelength),
                        FUN = switch,
                        `354` = 'purple',
                        `429` = 'blue',
                        `535` = 'darkgreen',
                        'red')
       )
       abline(h = c(0,1))
       lines(x = xl,
             y = plogis(xl*coef(lm_uv)[2]+
                          coef(lm_uv)[1]),
              col = 'magenta',
              lwd = 3)
       lines(x = xl,
             y = plogis(xl*coef(lm_bl)[2]+
                          coef(lm_bl)[1]),
              col = 'cyan',
              lwd = 3)
       lines(x = xl,
             y = plogis(xl*coef(lm_gr)[2]+
                          coef(lm_gr)[1]),
              col = 'green',
              lwd = 3)
     }
)


#rough estimate psychometric
ed_uv = within(ed_uv,
               {
               rel_logit_acc = qlogis(
                               (accuracy - min(accuracy))/
                                max(accuracy - min(accuracy))
                                     )
               }
              )

ed_bl = within(ed_bl,
               {
               rel_logit_acc = qlogis(
                               (accuracy - min(accuracy))/
                                max(accuracy - min(accuracy))
                                     )
               }
              )

ed_gr = within(ed_gr,
               {
               rel_logit_acc = qlogis(
                               (accuracy - min(accuracy))/
                                max(accuracy - min(accuracy))
                                     )
               }
              )

relm_uv = with(subset(ed_uv, is.finite(rel_logit_acc)),
             lm(rel_logit_acc ~ log10_intensity) 
)

relm_bl = with(subset(ed_bl, is.finite(rel_logit_acc)),
             lm(rel_logit_acc ~ log10_intensity) 
)

relm_gr = with(subset(ed_gr, is.finite(rel_logit_acc)),
             lm(rel_logit_acc ~ log10_intensity) 
)



with(ed,
     {
       plot(x = log10_intensity,
            y = accuracy,
            xlim = c(9,14),
            ylim = c(0,1),
            pch = 21,
            cex = 1.5,
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
with(ed_uv,
     {
      lines(x = xl,
           y = diff(range(accuracy))*
             plogis(xl*coef(relm_uv)[2]+
                        coef(relm_uv)[1]) +
             min(accuracy),
           col = 'magenta',
           lwd = 3)
     }
)
with(ed_bl,
     {
      lines(x = xl,
           y = diff(range(accuracy))*
             plogis(xl*coef(relm_bl)[2]+
                        coef(relm_bl)[1]) +
             min(accuracy),
           col = 'cyan',
           lwd = 3)
     }
)
with(ed_gr,
     {
      lines(x = xl,
           y = diff(range(accuracy))*
             plogis(xl*coef(relm_gr)[2]+
                        coef(relm_gr)[1]) +
             min(accuracy),
           col = 'green',
           lwd = 3)
     }
)



# Psychometric version ----------------------------------------------------


