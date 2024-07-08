#' Random Forest Imputation of Potential Outcomes in Paired Experiments (Pairs as Units)
#'
#' This function is used to impute potential outcomes for the \code{p_loop} function. 
#' Leaves out each pair and imputes its potential outcomes using random forests
#' on the remaining pairs. Pairs are treated as units when making predictions.
#' @param ordered A matrix of pair experimental data that his been processed by the \code{pair} function, with the treatment pair first.
#' @param assigned A matrix of pair experimental data that his been processed by the \code{pair} function.
#' @param n_assigned A matrix of pair experimental data cluster sizes that has been processed by the \code{pair} function.
#' @param reparm If set to \code{TRUE}, covariates will be parameterized as the pairwise differences and means for each covariate.
#' @param weighted_imp If set to \code{TRUE}, cluster sizes will be used as weights in imputation models.
#' @export

p_rf_v12 = function(ordered, assigned, n_assigned, reparm = TRUE, weighted_imp = F){
  
  if(reparm == TRUE){
    obs_v1 = reparam(assigned)
    #only differences between the covariates will now be the opposite sign
    obs_v2 = obs_v1 %>%
      mutate(across(ends_with("dif"), ~-.x))
    ordered = reparam(ordered)
  }else{
    obs_v1 = assigned
    #rename columns to switch the label between unit 1 and 2 in each pair
    obs_v2 = assigned %>%
      rename_with(~ str_replace(.x,"\\d$",as.character(3-as.numeric(str_extract(.x,"\\d$")))))
  }
  
  #model data with treated unit first, so sum is Y1+Y2 and difference is Y1-Y2
  mod_dat <- ordered %>%
    mutate(sum = Y1 + Y2,
           dif = Y1 - Y2) %>%
    select(-P, -Y1, -Y2)
  
  n1 = n_assigned$n1
  n2 = n_assigned$n2
  
  M = nrow(ordered)
  v1 = v2 = rep(0,M)
  dif_var = which(colnames(mod_dat) == "dif")
  sum_var = which(colnames(mod_dat) == "sum")
  
  for(i in 1:M){
    if(weighted_imp){
      w = mod_dat[-i,]$n1 + mod_dat[-i,]$n2
      
      rf_sum = randomForest::randomForest(sum ~ . -n1 -n2, mod_dat[-i,-dif_var], weights = w)
      rf_dif = randomForest::randomForest(dif ~ . -n1 -n2, mod_dat[-i,-sum_var], weights = w)
    }else{
      rf_sum = randomForest::randomForest(sum ~ . -n1 -n2, mod_dat[-i,-dif_var])
      rf_dif = randomForest::randomForest(dif ~ . -n1 -n2, mod_dat[-i,-sum_var])
    }
   
    v1[i] = (n1[i] - n2[i])/2*predict(rf_sum,obs_v1[i,]) + (n1[i] + n2[i])/2*predict(rf_dif,obs_v1[i,])
    v2[i] = (n2[i] - n1[i])/2*predict(rf_sum,obs_v2[i,]) + (n2[i] + n1[i])/2*predict(rf_dif,obs_v2[i,])
  }
  
  return(data.frame(v1,v2))
  
}
