#' OLS Imputation of Potential Outcomes in Paired Experiments
#'
#' This function is used to impute potential outcomes for the \code{p_loop} function.
#' Leaves out each pair and imputes its potential outcomes using ordinary least squares
#' regression on the remaining pairs. Pairs are treated as units when making predictions.
#' @param ordered A matrix of pair experimental data that his been processed by the \code{pair} function, with the treatment pair first.
#' @param assigned A matrix of pair experimental data that his been processed by the \code{pair} function.
#' @param n_assigned A matrix of pair experimental data cluster sizes that has been processed by the \code{pair} function.
#' @param weighted_imp A Boolean indicating whether to weight imputation models
#' @export

p_ols_v12 = function(ordered, assigned, n_assigned, weighted_imp = F){

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

  coefs_sum <- coefs_dif <- matrix(nrow=M,ncol=ncol(mod_dat)-3)

  for(i in 1:M){
    
    if(weighted_imp){
      ols_sum = lm(sum ~ . -n1 -n2, mod_dat[-i,-dif_var], weights = n1+n2)
      ols_dif = lm(dif ~ . -n1 -n2, mod_dat[-i,-sum_var], weights = n1+n2)
    }else{
      ols_sum = lm(sum ~ . -n1 -n2, mod_dat[-i,-dif_var])
      ols_dif = lm(dif ~ . -n1 -n2, mod_dat[-i,-sum_var])
    }
    
    coefs_sum[i,] <- coef(ols_sum)
    coefs_dif[i,] <- coef(ols_dif)

    v1[i] = (n1[i] - n2[i])/2*predict(ols_sum,obs_v1[i,]) + (n1[i] + n2[i])/2*predict(ols_dif,obs_v1[i,])
    v2[i] = (n2[i] - n1[i])/2*predict(ols_sum,obs_v2[i,]) + (n2[i] + n1[i])/2*predict(ols_dif,obs_v2[i,])
  }
  colnames(coefs_sum) <- colnames(coefs_dif) <- names(coef(ols_sum))


  out <- data.frame(v1,v2)

  attr(out,"coefs") <- list(coefs_sum=coefs_sum,coefs_dif=coefs_dif)

  return(out)
}
