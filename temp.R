Find confidence intervals for turn distributions

For each change distribution, find the maximum likelihood von Mises distribution. For the UV high to UV low change, this will be a bimodal distribution.

{r ML heading change}
#collect all heading differences

# Test for differences between conditions ---------------------------------
#collect all heading differences
pair_diffs_lst = lapply(list(full_mu_diff_gl,
                             full_mu_diff_uh,
                             full_mu_diff_ul,
                             full_mu_diff_gul),
                        FUN = unlist 
)
names(pair_diffs_lst) = c('green_hl',
                          'uvg_h',
                          'uvg_hl',
                          'uvg_l'
)
angle_unit = 'degrees'
angle_rot = 'clock'
#make sure they are in circular format
pair_diffs_lst = lapply(X = pair_diffs_lst,
                        FUN = circular,
                        units = angle_unit,
                        rotation = angle_rot)

#fit the maximum likelihood distribution for differences between pairs
# if angles shift in one direction between trials,
# this distribution should have significantly higher likelihood
ml_diff_lst = lapply(X = pair_diffs_lst,
                     mle.vonmises,
                     bias = TRUE)

#fit the maximum likelihood distribution for differences centred on zero
# if there is no consistent shift between trials,
# this distribution should have similar likelihood,
# with one less free parameter (expected mean of zero) 
ml_same_lst = lapply(X = pair_diffs_lst,
                     mle.vonmises,
                     bias = TRUE,
                     mu = circular(x = 0,
                                   units = angle_unit,
                                   rotation = angle_rot)
)

#inspect resulting parameters
rspar = do.call(rbind, c(ml_diff_lst, ml_same_lst))
rownames(rspar) = paste( sort(rep(c('diff', 'same'), length(pair_diffs_lst)) ),
                         rep(names(pair_diffs_lst), times = 2) )

print(rspar)

#also fit the maximum likelihood model
mod_diff_lst = lapply(X = pair_diffs_lst,
                      circ_mle)

# Extract model properties ------------------------------------------------
MD_extract = function(md)
{
  md_order = with(md, {rownames(results)})
  md_best = with(md, {results[md_best %in% bestmodel]})
  with(md_best,
       {
         return(
           list(q1,
                k1,
                q2,
                k2,
                w1 = lamda
                )
         )
       }
       )
}


#For two different populations with different concentrations around the mean
md_order = rownames(cmle$results)
m5b = cmle$results[md_order %in% 'M5B',]


## Calculate distance from mixture model means ----------------------------

#generic mean angle calculator
MeanRvm = function(n, mu = circular(0), kappa, au = 'degrees', ar = 'clock')
{
  mean.circular(rvonmises(n = n, 
                          mu = circular(mu, units = au, rotation = ar), 
                          kappa = kappa,
                          control.circular = list(units = au, rotation = ar)))
}


#simulate von Mises draws from sample
SimVM_CI = function(angles, 
                    mu = NULL,
                    n = 1e4,
                    au = 'degrees', 
                    ar = 'clock',
                    calc_q = TRUE,
                    speedup_parallel = TRUE
)
{
  ml = mle.vonmises(x = circular(x = angles,
                                 units = au,
                                 rotation = ar),
                    bias = TRUE, 
                    mu = if(!is.null(mu))
                    {
                      circular(x = mu,
                               units = au,
                               rotation = ar
                      )
                    }else{mu}
  )
  if(speedup_parallel) #3x faster
  {
    cl = parallel::makePSOCKcluster(parallel::detectCores()-1)
    parallel::clusterExport(cl = cl, 
                            varlist = c('mean.circular',
                                        'circular',
                                        'rvonmises'),
                            envir = .GlobalEnv
    )
    parallel::clusterExport(cl = cl, 
                            varlist = c('MeanRvm',
                                        'angles',
                                        'mu',
                                        'n',
                                        'au',
                                        'ar'),
                            envir = environment()
    )
    smp = with(ml, 
               parallel::parSapply(cl = cl,
                                   X = 1:n,
                                   FUN = function(i)
                                   {
                                     eval.parent(
                                       {
                                         MeanRvm(n = length(angles), 
                                                 mu = mu, 
                                                 kappa = kappa,
                                                 au = au,
                                                 ar = ar)
                                       }
                                     )
                                   },
                                   simplify = 'array'
               )
    )
    parallel::stopCluster(cl)
  }else
  {
    smp = with(ml, 
               replicate(n = n, 
                         MeanRvm(n = length(angles), 
                                 mu = mu, 
                                 kappa = kappa,
                                 au = au,
                                 ar = ar)
               )
    )
  }
  return(
    if(calc_q)
    {
      Mod360.180(
        quantile.circular(x = circular(x = smp,
                                       units = au,
                                       rotation = ar),
                          probs = sort(c(c(0,1)+c(1,-1)*0.05, 0.5)))
      )
    }else
    {
      sapply(X = smp, FUN = Mod360.180)
    }
  )
}#simulate von Mises draws from sample

#TODO make version for bimodal data
SimMVM_CI = function(angles, 
                    mu = NULL,
                    n = 1e4,
                    au = 'degrees', 
                    ar = 'clock',
                    calc_q = TRUE,
                    speedup_parallel = TRUE,
                    ... #input to circ_mle
)
{
  mod = circ_mle(x = circular(x = angles,
                                 units = au,
                                 rotation = ar),
                 ...
  )
  
  mod_ex = MD_extract(mod)
  
  if(speedup_parallel) #3x faster
  {
    cl = parallel::makePSOCKcluster(parallel::detectCores()-1)
    parallel::clusterExport(cl = cl, 
                            varlist = c('mean.circular',
                                        'circular',
                                        'rvonmises'),
                            envir = .GlobalEnv
    )
    parallel::clusterExport(cl = cl, 
                            varlist = c('MeanRvm',
                                        'angles',
                                        'mu',
                                        'n',
                                        'au',
                                        'ar'),
                            envir = environment()
    )
    smp = with(mod_ex, 
               parallel::parSapply(cl = cl,
                                   X = 1:n,
                                   FUN = function(i)
                                   {
                                     eval.parent(
                                       {
                                         c(
                                         m1 = MeanRvm(n = round(
                                                           length(angles)*w1 ), 
                                                 mu = q1, 
                                                 kappa = k1,
                                                 au = au,
                                                 ar = ar),
                                         m2 = MeanRvm(n = round(
                                                           length(angles)*
                                                                   (1-w1) ), 
                                                 mu = q2, 
                                                 kappa = k1,
                                                 au = au,
                                                 ar = ar)
                                         )
                                       }
                                     )
                                   },
                                   simplify = 'array'
               )
    )
    parallel::stopCluster(cl)
  }else
  {
    smp = with(ml, 
               replicate(n = n, 
                         MeanRvm(n = length(angles), 
                                 mu = mu, 
                                 kappa = kappa,
                                 au = au,
                                 ar = ar)
               )
    )
  }
  return(
    if(calc_q)
    {
      Mod360.180(
        quantile.circular(x = circular(x = smp,
                                       units = au,
                                       rotation = ar),
                          probs = sort(c(c(0,1)+c(1,-1)*0.05, 0.5)))
      )
    }else
    {
      sapply(X = smp, FUN = Mod360.180)
    }
  )
}


# Plot the ML distributions -----------------------------------------------



## Calculate likelihood ------------------------------------------------------------

#function to calculate the likelihood of a sample, given the input ML parameters
LLcalc = function(ml,
                  angles,
                  au = 'degrees',
                  ar = 'clock',
                  ...)
{
  with(ml,
       {
         sum( # add together
           dvonmises(x = circular(x = angles,
                                  units = au,
                                  rotation = ar), # probability density for each observed angle
                     mu = mu, # ML estimated mean
                     kappa = kappa, # ML estimated concentration
                     log = TRUE) # on a log scale (i.e. add instead of multiplying)
         )
       }
  )
}
#calculate all in list
all_ll = mapply(FUN = LLcalc,
                ml = c(ml_diff_lst, ml_same_lst), #apply to both sets of ML estimates
                angles = pair_diffs_lst, #repeat datasets
                au = angle_unit, #important, use the same frame of reference as the original
                ar = angle_rot)
#add statistics
all_ll = data.frame(ll = all_ll,
                    deviance = -all_ll*2,#calculate deviance
                    model = c(rep('diff',4), #label the models
                              rep('same',4) ),
                    dataset = names(all_ll) #the dataset names were saved in the output
)

# Perform likelihood ratio test -------------------------------------------

#The likelihood ratio test compares models by estimating the likelihood gained
#for each additional parameter.
#https://en.wikipedia.org/wiki/Likelihood-ratio_test
# In our case, the ML von Mises requires two parameters (mean and concentration), 
# while the von Mises around the expected mean requires one (concentration).
# How much does the extra parameter increase the likelihood?

CollectDetails = function(ll, dts = 'all', ...)
{
  with( subset(ll, dataset == dts),
        {
          data.frame(
            modnm = c('pairs diff', 'pairs same'),
            ll = c(ll[model %in% 'diff'],
                   ll[model %in% 'same']),
            deviance = c(deviance[model %in% 'diff'], 
                         deviance[model %in% 'same']),
            rnk = rank(c(deviance[model %in% 'diff'], 
                         deviance[model %in% 'same']) ),#paired differences use fewer observations, don't include in ranking 
            df = c(2, 1)
          )
        }
  )
}

#compile model details
mod_details = lapply(X = names(pair_diffs_lst), #search for the dataset names
                     FUN = CollectDetails,
                     ll  = all_ll) #likelihood data comes from our statistics
names(mod_details) = names(pair_diffs_lst) #names are not read out this time

#set up tests for hypotheses
lr_tests = 'pairs_diff_zero' #just this one test (could also test against uniform)

#Perform collect likelihood ratio tests for each hypothesis
all_results = data.frame(datasets = names(pair_diffs_lst),
                         t(
                           sapply(X =mod_details,
                                  FUN = LR_calc,
                                  tst = lr_tests)
                         )
)

#add labels to interpret each test according to p value and sign of comparison
all_results = within(all_results,
                     {
                       result = mapply(tst = lr_tests,
                                       d0 = unlist(dev0),
                                       d1 = unlist(dev1),
                                       pa = p_adjusted,
                                       FUN = H1label
                       )
                     }
)
#add calculated means
Exmu = function(ml){round(ml$mu, 3)}
all_results = within(all_results,
                     {
                       angular_diff = 
                         sapply(X = ml_diff_lst,
                                FUN = Exmu
                         )
                     }
)
#print the results for the user                     
with(all_results,
     {
       print(
         data.frame(dataset = datasets,
                    chi_squared = round(unlist(chi_squared), 3),
                    d.f = unlist(d.f.),
                    p = round(unlist(p_adjusted), 4),
                    result = result)
       )
     }
)

# Save result -------------------------------------------------------------
#save result
res_table = apply(X = all_results,
                      MARGIN = 1:2,
                      FUN = unlist)