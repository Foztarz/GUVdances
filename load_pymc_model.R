
# Load the required packages ----------------------------------------------
# require(reticulate)
# reticulate::py_install(c(
#   "numpy==1.26.4",
#   "pandas==2.1.4",
#   "scipy==1.11.4",
#   "arviz==0.17.1",
#   "packaging"
# ), method = "pip")
require(reticulate)
require(posterior)
require(bayesplot)
# require(ncdf4)
# require(cmdstanr)


# Import PyMC model -------------------------------------------------------

# Load the NetCDF file
mod_path = file.path("C:/Users/James Foster/Documents/GitHub/GUVdances/",
                     "GUV_unwrap_model.nc")

az = reticulate::import("arviz")
idata = az$from_netcdf(mod_path)
posterior_df = reticulate::py_to_r(idata$posterior$to_dataframe())

# Initial Metadata -----------------------------------------------------

n_ind = 169
n_chains = 4
obs_per_chain = nrow(posterior_df) /
                    (n_ind * n_chains) # Calculated automatically

# Individual changes every row (1,2...169, 1,2...169)
posterior_df = within(posterior_df,
                      {
                Individual_ID = rep(x = 1:n_ind, 
                         times = n_chains * obs_per_chain)
    
    # Iteration changes every 169 rows (1,1... (169 times), 2,2... (169 times))
                .iteration    = rep(x = rep(x = 1:obs_per_chain, 
                                 each = n_ind), 
                         times = n_chains)

# Chain changes every (169 * obs_per_chain) rows
                .chain        = rep(x = 1:n_chains, 
                                    each = obs_per_chain * n_ind)
                      }
)
# Define Parameter Groupings
global_params = c(
  "Intercept", "B", "C", "B:C", 
  "kappa_Intercept", "kappa_B", "kappa_C", "kappa_B:C",
  "1|Individual_kappa", "B|Individual_kappa", "C|Individual_kappa", "B:C|Individual_kappa",
  "kappa_1|Individual_sigma", "kappa_B|Individual_sigma", 
  "kappa_C|Individual_sigma", "kappa_B:C|Individual_sigma",
  "kappa_1|Individual_offset", "kappa_B|Individual_offset", 
  "kappa_C|Individual_offset", "kappa_B:C|Individual_offset"
)

ind_params = c(
  "1|Individual", "B|Individual", "C|Individual", "B:C|Individual",
  "kappa_1|Individual", "kappa_B|Individual", "kappa_C|Individual", "kappa_B:C|Individual"
)


# Reshape posterior draws -------------------------------------------------

# Reshape to Wide Format
list_of_inds = split(x = posterior_df, 
                     f = posterior_df$Individual_ID)

# Start wide df with metadata and population-level parameters
existing_globals = intersect(x = global_params,
                             y = colnames(posterior_df))
# take the first individual's data as a template for the wide data frame
draws_wide  =  list_of_inds[[1]][, c(".chain", ".iteration", existing_globals)]

# Pivot individual-level parameters into columns
for(param in ind_params)
  {
  if(param %in% colnames(posterior_df))
    {
    clean_base = gsub("\\|", "_", param)
    for(i in 1:n_ind) {
      draws_wide[[paste0(clean_base, "_Ind_", i)]] = list_of_inds[[i]][[param]]
    }
  }
}

# Export to draws dataframe
draws_final = posterior::as_draws_df(draws_wide)


# Plot draws --------------------------------------------------------------
require(circular)
unwrap_circular_deg = function(x)
{
  mux = mean.circular(x = circular(x = x, template = 'none'))
  centx = atan2(y = sin(x - mux),
                x = cos(x  - mux))
  unwrx = centx + mux
  return(deg(unwrx))
}


mcmc_trace(draws_final,
           pars = c('Intercept', 'B', 'C', 'B:C'),
           transformations = unwrap_circular_deg)

mcmc_trace(draws_final,
           pars = paste0("1_Individual_Ind_", sample(1:169, 6)),
           transformations = unwrap_circular_deg)

mcmc_trace(draws_final,
           pars = paste0("kappa_",
                         c('Intercept', 'B', 'C', 'B:C')) )
mcmc_trace(draws_final,
           pars = paste0( c('1', 'B', 'C', 'B:C'),
                          '|Individual_kappa') )
mcmc_trace(draws_final,
           pars = paste0('kappa_',
                         c('1', 'B', 'C', 'B:C'),
                         '|Individual_sigma') )



# Extract predictions -----------------------------------------------------

#A function to find the circular median in (-pi, pi)
MedCDraws = function(x,
                     templ = 'none', #assumes radians
                     medn = 1) #take the 1st median
{
  median.circular( circular(x = mod_circular(x),
                            template = templ)
  )[medn]
}
QCDraws = function(x,
                   probs = c(0.025, 0.5, 0.975),
                   templ = 'none',
                   medn = 1)
{
  quantile.circular( circular(x = mod_circular(x),
                              template = templ),
                     probs = probs
  )
}

# Setup containers
u_id <- 1:n_ind
mu_id_gh <- kappa_id_gh <- numeric(length(u_id))
mu_id_uh <- kappa_id_uh <- numeric(length(u_id))
mu_id_gl <- kappa_id_gl <- numeric(length(u_id))
mu_id_ul <- kappa_id_ul <- numeric(length(u_id))
for(i in 1:length(u_id)) {
  
  # 1. Construct the exact column names
  # Note: These must match colnames(draws_final) exactly
  mu_i   <- paste0("1_Individual_Ind_", i)
  mu_bi  <- paste0("B_Individual_Ind_", i)
  mu_ci  <- paste0("C_Individual_Ind_", i)
  mu_bci <- paste0("B:C_Individual_Ind_", i)
  
  ka_i   <- paste0("kappa_1_Individual_Ind_", i)
  ka_bi  <- paste0("kappa_B_Individual_Ind_", i)
  ka_ci  <- paste0("kappa_C_Individual_Ind_", i)
  ka_bci <- paste0("kappa_B:C_Individual_Ind_", i)
  
  # 2. Extract and Calculate
  # Using draws_final[[name]] extracts the column as a numeric vector
  
  # Green High (Reference)
  mu_id_gh[i] <- MedCDraws(draws_final[["Intercept"]] + draws_final[[mu_i]])
  kappa_id_gh[i] <- median(draws_final[["kappa_Intercept"]] + draws_final[[ka_i]])
  
  # UV High (C effect)
  mu_id_uh[i] <- MedCDraws(draws_final[["Intercept"]] + draws_final[[mu_i]] + 
                             draws_final[["C"]] + draws_final[[mu_ci]])
  
  kappa_id_uh[i] <- median(draws_final[["kappa_Intercept"]] + draws_final[[ka_i]] + 
                             draws_final[["kappa_C"]] + draws_final[[ka_ci]])
  
  # Green Low (B effect)
  mu_id_gl[i] <- MedCDraws(draws_final[["Intercept"]] + draws_final[[mu_i]] + 
                             draws_final[["B"]] + draws_final[[mu_bi]])
  
  kappa_id_gl[i] <- median(draws_final[["kappa_Intercept"]] + draws_final[[ka_i]] + 
                             draws_final[["kappa_B"]] + draws_final[[ka_bi]])
  
  # UV Low (All effects + Interactions)
  mu_id_ul[i] <- -MedCDraws(draws_final[["Intercept"]] + draws_final[[mu_i]] + 
                             draws_final[["B"]] + draws_final[["C"]] + draws_final[["B:C"]] + 
                             draws_final[[mu_bi]] + draws_final[[mu_ci]] + draws_final[[mu_bci]])
  
  kappa_id_ul[i] <- median(draws_final[["kappa_Intercept"]] + draws_final[[ka_i]] + 
                             draws_final[["kappa_B"]] + draws_final[["kappa_C"]] + draws_final[["kappa_B:C"]] + 
                             draws_final[[ka_bi]] + draws_final[[ka_ci]] + draws_final[[ka_bci]])
}
# Plot predictions --------------------------------------------------------


#the radian equivalent
mod_circular = function(x)
{
  atan2(y = sin(x),
        x = cos(x))
}

OpenCplot = function(x,
                     angle_unit = 'degrees',
                     angle_rot = 'clock',
                     n_sample = 10
)
{
  plot.circular(x = circular(x = NULL, 
                             type = 'angles',
                             unit = angle_unit,
                             # template = 'geographics',
                             modulo = '2pi',
                             zero = pi/2,
                             rotation = angle_rot
  ),
  col = NA
  )
  lines.circular(x = circular(x = 
                                seq(from = -180,
                                    to = 180,
                                    length.out = 1e3),
                              units = angle_unit, 
                              rotation = angle_rot),  
                 y = rep(x = sqrt(-log(0.05)/n_sample)-1,times = 1e3),
                 col = 'black', 
                 lty = 2,
                 lwd = 0.25,
                 lend = 'butt')
}

#add contours
Draws2Cont = function(draws,
                      palette = 'Heat 2',
                      nlevels = 20,
                      x_string = 'sin(Intercept)*A1(softplus(Intercept_kappa))',
                      y_string = 'cos(Intercept)*A1(softplus(Intercept_kappa))',
                      alpha = 200/255,
                      ngrid = 25, # defaults to a 25x25 grid
                      cropc = FALSE, #crop region outside circle
                      denstype = 'relative' # 'normalised' or 'relative' (normalised fails to plot low densities)
)
{
  kdc = with(draws,
             {
               MASS::kde2d(x = eval(str2lang(x_string)),
                           y = eval(str2lang(y_string)),
                           n = ngrid)
             }
  )
  if(cropc)
  {
    xy = with(kdc, expand.grid(x = x,y= y))#find coordinates of z variable
    idc = with(xy, x^2+y^2 < 1.0)
    #crop edge of circle
    kdc = within(kdc,
                 {z[!idc] = 0})
  }
  with(kdc,
       {
         .filled.contour(x = x,
                         y = y,
                         z = z,
                         levels = (1:nlevels)*
                           switch(EXPR = denstype,
                                  relative = max(z/nlevels),
                                  normalised = sum(z/nlevels),#warning, fails to plot low densities
                                  max(z/nlevels)
                           ),
                         col = hcl.colors(n = nlevels,
                                          palette = palette,
                                          rev = TRUE,
                                          alpha = alpha)
         )
       }
  )
}


# Population example ------------------------------------------------------


OpenCplot()

# Inside your plotting loop:
Draws2Cont(draws_final,
           x_string = paste0('sin( (Intercept)) * A1(softplus(kappa_Intercept))'),
           y_string = paste0('cos( (Intercept)) * A1(softplus(kappa_Intercept))')
             )

with(draws_final,
     {
  arrows.circular(
    x = circular(MedCDraws(Intercept),
                 units = 'radians',
                 rotation = 'clock',
                 zero = pi/2),
    y = A1(softplus(median(kappa_Intercept))), # Ensure softplus is applied here
    lwd = 5,
    length = 0.1/1.25,
    col = adjustcolor('green3', alpha.f = 0.8)
  )
}
)

# Individual example ------------------------------------------------------


OpenCplot()

# Inside your plotting loop:
Draws2Cont(draws_final,
           x_string = paste0('sin( -1*(Intercept + `',
                             mu_i,
                             '`)) * A1(softplus(kappa_Intercept + ',
                             ka_i,
                             '))'),
           y_string = paste0('cos(-1*(Intercept + `', 
                             mu_i, 
                             '`)) * A1(softplus(kappa_Intercept + ', 
                             ka_i, 
                             '))')
)

arrows.circular(
  x = circular(mu_id_gh[i],
               units = 'radians',
               rotation = 'clock',
               zero = pi/2),
  y = A1(softplus(kappa_id_gh[i])), # Ensure softplus is applied here
  lwd = 5,
  length = 0.1/1.25,
  col = adjustcolor('green3', alpha.f = 0.8)
)
