Find confidence intervals for turn distributions

For each change distribution, find the maximum likelihood von Mises distribution. For the UV high to UV low change, this will be a bimodal distribution.

{r ML heading change}

pair_diffs_lst = with(rdata,
                      {
                        list(green_hl = mu_diff_gl,
                             uvg_h = mu_diff_ugh,
                             uv_hl = mu_diff_ul,
                              uvg_l = mu_diff_ugl 
                        )
                      }
)

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

