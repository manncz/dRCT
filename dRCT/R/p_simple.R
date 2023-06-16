#' Simple Unbiased Estimator for Paired Clustered Randomized Trial
#'
#' This function is used to impute potential outcomes for the \code{p_loop} function. Pairs are treated
#' as the unit, and leave one out mean imputation is used, with no covariate adjustment.
#' @param ordered A matrix of pair experimental data that his been processed by the \code{pair} function, with the treatment pair first.
#' @param assigned A matrix of pair experimental data that his been processed by the \code{pair} function.
#' @param n_assigned A matrix of pair experimental data cluster sizes that has been processed by the \code{pair} function.
#' @export

p_simple <- function(ordered, assigned, n_assigned, true_sum = NULL){
  
  M <- nrow(assigned)
  v1 <- v2 <- rep(NA, M)
  
  n1 = n_assigned$n1
  n2 = n_assigned$n2
  
  for(i in 1:M){
    #if no differing weights, this will get us back to the simple difference estimator
    #this is not really necessary since d=0 in this case, but flows through to the variance estimator
    loo_mean_sum <- mean((ordered$Y1[-i] + ordered$Y2[-i])/2)
    loo_mean_dif <- mean((ordered$Y1[-i] - ordered$Y2[-i])/2)
    
    if(!isnull(true_sum)) loo_mean_sum = true_sum[i]
      
    v1[i] <- (n1[i] - n2[i])*loo_mean_sum + (n1[i] + n2[i])*loo_mean_dif
    v2[i] <- (n2[i] - n1[i])*loo_mean_sum + (n1[i] + n2[i])*loo_mean_dif
    
  }
  
  return(data.frame(v1,v2))
}