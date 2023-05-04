#' OLS Imputation of Potential Outcomes in Paired Experiments
#'
#' This function is used to impute potential outcomes for the \code{p_loop} function. 
#' Leaves out each pair and imputes its potential outcomes using ordinary least squares
#' regression on the remaining pairs. Pairs are treated as units when making predictions.
#' @param ordered A matrix of pair experimental data that his been processed by the \code{pair} function, with the treatment pair first.
#' @param assigned A matrix of pair experimental data that his been processed by the \code{pair} function.
#' @param n_assigned A matrix of pair experimental data cluster sizes that has been processed by the \code{pair} function.
#' @export

p_ols_v12 = function(ordered, assigned, n_assigned){
  
  #set up observed data for estimating v1 and v2
  obs_v1 = assigned
  
  #rename columns to switch the label between unit 1 and 2 in each pair
  obs_v2 = assigned %>%
    rename_with(~ str_replace(.x,"\\d$",as.character(3-as.numeric(str_extract(.x,"\\d$")))))
  
  n1 = n_assigned$n1
  n2 = n_assigned$n2
  
  #model data with treated unit first, so sum is Y1+Y2 and difference is Y1-Y2
  mod_dat <- ordered %>%
    mutate(sum = Y1 + Y2,
           dif = Y1 - Y2) %>%
    select(-P, -Y1, -Y2)
  
  #setup for leave-one-pair-out estimation
  M = nrow(ordered)
  v1 = v2 = rep(0,M)
  dif_var = which(colnames(mod_dat) == "dif")
  sum_var = which(colnames(mod_dat) == "sum")
  
  for(i in 1:M){
    ols_sum = lm(sum ~ . , mod_dat[-i,-dif_var])
    ols_dif = lm(dif ~ . , mod_dat[-i,-sum_var])
    
    v1[i] = (n1[i] - n2[i])/2*predict(ols_sum,obs_v1[i,]) + (n1[i] + n2[i])/2*predict(ols_dif,obs_v1[i,])
    v2[i] = (n2[i] - n1[i])/2*predict(ols_sum,obs_v2[i,]) + (n2[i] + n1[i])/2*predict(ols_dif,obs_v2[i,])
  }
  return(data.frame(v1,v2))
}
