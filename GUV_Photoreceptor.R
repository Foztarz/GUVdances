require(sfsmisc)


# Define a function for photoreceptor sensitivity -------------------------


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


# Define a function for photoreceptor response ----------------------------
#params based on 250327_006 from Kolyfetis et al., 2025
#https://docs.google.com/spreadsheets/d/18r4MJe30J2Hgh7MQHX0us8KqTKhufxmgN_C4IpbRGmc/edit?gid=2023063010#gid=2023063010&range=A1:F1
PRresponse = function(
            quanta = 10^seq(from = 10, to  = 14, length.out = 21), #photons/cm2/s
            vmax = 50,#mV
            k50 = 0.03*1e14,#inflection point
            n = 1.07)#slope
  {
mV = (vmax*quanta^n)/ 
      (k50^n + quanta^n)
}
#this sets sensitivity in the following way
## rel_quanta = quanta/max(quanta)
## yy = mV/vmax
## sensitivity = k50*(yy/(1-yy))^(n^-1)


# Set up wavelengths of interest ------------------------------------------
# Response peak for honeybee photoreceptors Menzel et al., 1986
lmax_u = 343
lmax_g = 543
# Emission peak for each LED
led_u = 365
led_g = 528
# A 20nm (uniform) band for each 
wlband = seq(from = -10, to  = 10, by = 0.1)


# Predict sensitivity curves -------------------------------------------------
# across whole spectrum
u_photopig = predict(StavengaSpline(lambda_max = lmax_u))
g_photopig = predict(StavengaSpline(lambda_max = lmax_g))
# specifically for the LEDs
#green photoreceptor
gPR_uLED = predict(StavengaSpline(lambda_max = lmax_g), x = led_u + wlband)
gPR_gLED = predict(StavengaSpline(lambda_max = lmax_g), x = led_g + wlband)
#UV photoreceptor
uPR_gLED = predict(StavengaSpline(lambda_max = lmax_u), x = led_g + wlband)
uPR_uLED = predict(StavengaSpline(lambda_max = lmax_u), x = led_u + wlband)

# integrated sensitivity
sens_gPR_uLED = integrate.xy(gPR_uLED)
sens_gPR_gLED = integrate.xy(gPR_gLED)
g_rel_u = sens_gPR_uLED/sens_gPR_gLED

# Plot response curves ----------------------------------------------------
xx = 10^seq(from = 11, to  = 15, length.out = 1e3) #photons/cm2/s
example_g_equal = PRresponse(k50 = 3e12/g_rel_u, quanta = xx)
example_g_lows = PRresponse(k50 = 3e13/g_rel_u, quanta = xx)
example_u = PRresponse(k50 = 3e12, quanta = xx)

#set margins all around
par(mar = c(4.5,4.5,4.5,4.5))


## Assume equal absolute sensitivity -------------------------------------
#open plot with green photoreceptor
plot(x = xx,
     y = example_g_equal,
     type = 'l',
     col = 'green',
     lwd = 5,
     log = 'x',
     ylim = c(0, 50),
     main = 'equal absolute sensitivity',
     xlab = 'UV intensity (photons/cm2/s)',
     ylab = 'photoreceptor response (mV)')
#add UV photoreceptor
lines(x = xx,
     y = example_u,
     col = 'magenta',
     lwd = 5)
#add reference lines
abline(h = c(0,0.5,1.0)*50,
       v = 10^c(11:15),
       col = 'gray')
#add LEDs
abline(v = c(1.2e12, 4e12, 4e14),
       col = c('magenta4', 'purple', 'magenta'))
#label curves
text(x = 10^c(12, 13.5, 13.7),
     y = c(25, 25, 30),
     labels = c('UV\nphotor.',
                'green\nphotor.',
                'log ratio'),
     col = c('magenta', 'green', 'purple4')
     )
#add log ratio
#open new plot
par(new = TRUE)
plot(x = xx,
     y = log(example_g_equal/example_u),
     type = 'l',
     ylim = c(-1.55, 0),
     log = 'x',
     lwd = 5,
     col = 'purple4',
     axes = FALSE,
     xlab = '',
     ylab = '')
#add axis
axis(side = 4,
     col = 'purple4')
#add axis label
mtext(side = 4,
      line = 2,
      text = 'log(green/UV)',
      col = 'purple4')


## Account for sensitivity difference -----------------------------------


plot(x = xx,
     y = example_g_lows,
     type = 'l',
     col = 'green',
     lwd = 5,
     log = 'x',
     ylim = c(0, 50),
     main = 'green lower absolute sensitivity',
     xlab = 'UV intensity (photons/cm2/s)',
     ylab = 'photoreceptor response (mV)')
lines(x = xx,
     y = example_u,
     col = 'magenta',
     lwd = 5)
#add reference lines
abline(h = c(0,0.5,1.0)*50,
       v = 10^c(11:15),
       col = 'gray')
#add LEDs
abline(v = c(1.2e12, 4e12, 4e14),
       col = c('magenta4', 'purple', 'magenta'))
#label curves
text(x = 10^c(12.2, 14.4, 13.7),
     y = c(25, 25, 25),
     labels = c('UV\nphotor.',
                'green\nphotor.',
                'log ratio'),
     col = c('magenta', 'green', 'purple4')
)
#add log ratio
#open new plot
par(new = TRUE)
plot(x = xx,
     y = log(example_g_lows/example_u),
     type = 'l',
     ylim = c(-4.0, 0),
     log = 'x',
     lwd = 5,
     col = 'purple4',
     axes = FALSE,
     xlab = '',
     ylab = '')

axis(side = 4,
     col = 'purple4')
mtext(side = 4,
      line = 2,
      text = 'log(green/UV)',
      col = 'purple4')


# Plot for Fig. 4 ----------------------------------------------------------
par(mfrow = c(1,2),
    mar = c(5,5,0,3))
plot(u_photopig,
     type = 'l',
     lwd = 5,
     col = 'magenta',
     xlim = c(300, 700),
     ylim = c(0,1),
     xlab = 'wavelength (nm)',
     ylab = 'photoreceptor sensitivity\n(relative)',
     )
lines(g_photopig,
      lwd = 5,
      col = 'green')

with(uPR_uLED,
  polygon(x = c(x, rev(x)),
          y = c(y, rep(0, length(y))),
          col = adjustcolor('magenta',
                            alpha.f = 0.2),
          border = NA)
)
with(uPR_gLED,
  polygon(x = c(x, rev(x)),
          y = c(y, rep(0, length(y))),
          col = adjustcolor('magenta',
                            alpha.f = 0.2),
          border = NA)
)

with(gPR_uLED,
  polygon(x = c(x, rev(x)),
          y = c(y, rep(0, length(y))),
          col = adjustcolor('green',
                            alpha.f = 0.2),
          border = NA)
)
with(gPR_gLED,
  polygon(x = c(x, rev(x)),
          y = c(y, rep(0, length(y))),
          col = adjustcolor('green',
                            alpha.f = 0.2),
          border = NA)
)
abline(h = c(0,1.0),
       v = 100*c(3:7),
       col = 'gray')


plot(x = xx,
     y = example_g_equal,
     type = 'l',
     col = 'green',
     lwd = 5,
     log = 'x',
     ylim = c(0, 50),
     xlab = 'UV intensity (photons/cm2/s)',
     ylab = 'photoreceptor response (mV)')
lines(x = xx,
      y = example_u,
      col = 'magenta',
      lwd = 5)
#add reference lines
abline(h = c(0,0.5,1.0)*50,
       v = 10^c(11:15),
       col = 'gray')
#add LEDs
abline(v = c(1.2e12, 4e12, 4e14),
       col = c('magenta4', 'purple', 'magenta'))
#label curves
text(x = 10^c(12.2, 14.4, 13.7),
     y = c(25, 25, 25),
     labels = c('UV\nphotor.',
                'green\nphotor.',
                'ratio'),
     col = c('magenta', 'green', 'purple4')
)
#add ratio
#open new plot
par(new = TRUE)
plot(x = xx,
     y = example_g_equal/example_u,
     type = 'l',
     ylim = c(0, 1),
     log = 'x',
     lwd = 5,
     col = 'purple4',
     axes = FALSE,
     xlab = '',
     ylab = '')

axis(side = 4,
     col = 'purple4')
mtext(side = 4,
      line = 2,
      text = 'green/UV',
      col = 'purple4')
