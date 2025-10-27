

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


# Load data -----
#Dance angle data
cd = read.table(file = '1Data/colour_dance_reorg.csv', 
                header = T, 
                sep  = ',')

# Data formatting ------

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

# Plot all dances in one figure -------------------------------------------
cx = 0.3



#set up a sequence for the axes
xc = seq(from = -pi, to = pi-1e-16, length.out = 1e3)
par(mfrow = c(22, 8),
# par(mfrow = c(8, 22), # landscape version
    mar = c(0,0,0,0))
par(pty = 's')

#add legend to 1st page
plot(x = NULL,
     xlim = c(-1,1),
     ylim = c(-1,1),
     pch = 19,
     axes = FALSE,
     xlab = '',
     ylab = '',
     main = '',
     cex = cx
)
legend(x = 'center',
       legend = c('Green Bright',
                  'Green Dim',
                  'UV Bright',
                  'UV Dim',
                  'Implied North\n(-sun azimuth)'),
       col = c('green',
               'darkgreen',
               'magenta',
               'purple',
               'darkred'),
       pch = c(21,21,21,21,NA),
       lty = c(NA, NA, NA, NA, 1),
       lwd = c(2,2,2,2,3),
       bty = 'n',
       cex = cx
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
       main = '',
       cex = cx
  )
  abline(a = 0, b = 1, col = gray(0.9, alpha = cx), lwd = cx)
  abline(a = 0, b = -1, , col = gray(0.9, alpha = cx), lwd = cx)
  abline(h = 0, v = 0, , col = gray(0.75, alpha = cx), lwd = cx)
  lines(x = 20*sin(xc), y = 20*cos(xc), lty = 3, lwd = cx, col = gray(0, alpha = cx))
  lines(x = 40*sin(xc), y = 40*cos(xc), lty = 3, lwd = cx, col = gray(0, alpha = cx))
  lines(x = 60*sin(xc), y = 60*cos(xc), lty = 3, lwd = cx, col = gray(0, alpha = cx))
  lines(x = 80*sin(xc), y = 80*cos(xc), lty = 3, lwd = cx, col = gray(0, alpha = cx))
  text(x = 1:4*20,
       y = c(0,0,0,0),
       labels = paste(' run', 1:4*20), 
       cex = 0.35*cx,
       adj = c(0,1))
  with(subset(cd, ID %in% ii),
       {
         lines(x = c(0, 100*sin(-sun_az_rad)),
               y = c(0, 100*cos(-sun_az_rad)),
               col = adjustcolor('darkred',alpha.f = 0.5),
               lwd = 2*cx)
         points(x = run*sin(as.numeric(angle)),
                y = run*cos(as.numeric(angle)),
                bg = gray(level = 1.0,
                          alpha =  0.4),
                col = c('green', 'darkgreen', 'magenta', 'purple')
                [ifelse(colour %in% 'g',
                        yes = ifelse(brightn %in% 'h', yes = 1, no = 2),
                        no = ifelse(brightn %in% 'h', yes = 3, no = 4))],
                pch = 21,
                lwd = cx,
                # lwd = 2*cex,
                cex = cx*0.5
         )
       }
  )
  lines(x = c(0,100*sin(mnv_gh['mu'])*mnv_gh['rho']), 
        y = c(0,100*cos(mnv_gh['mu'])*mnv_gh['rho']), 
        col = 'green', 
        lwd = 1*cx)
  lines(x = c(0,100*sin(mnv_gl['mu'])*mnv_gl['rho']), 
        y = c(0,100*cos(mnv_gl['mu'])*mnv_gl['rho']), 
        col = 'darkgreen', 
        lwd = 1*cx)
  lines(x = c(0,100*sin(mnv_uh['mu'])*mnv_uh['rho']), 
        y = c(0,100*cos(mnv_uh['mu'])*mnv_uh['rho']), 
        col = 'magenta', 
        lwd = 1*cx)
  lines(x = c(0,100*sin(mnv_ul['mu'])*mnv_ul['rho']), 
        y = c(0,100*cos(mnv_ul['mu'])*mnv_ul['rho']), 
        col = 'purple', 
        lwd = 1*cx)
  mtext(ii,side = 1,line = -1, cex = cx)
}
