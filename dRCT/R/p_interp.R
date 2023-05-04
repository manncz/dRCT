#' Interpolates between Sets of Potential Differences in Paired Experiments
#'
#' This function is used within the \code{p_rf_interp} and \code{p_ols_interp} imputation methods for the \code{p_loop} function. 
#' Given two sets of imputed potential differences, this function returns a single set of predicted differences by picking weights
#' that minimizes the mean squared error between the observed differences and the weighted average imputed differences. This minimization
#' is done using leave-one-out cross validation.
#' 
#' @param Tr The treatment assignment vector for the pairs, where the value is 1 if the first unit in the pair is assigned to treatment and 0 otherwise.
#' @param v1 A vector of the treatment minus control differences for the pairs where the first unit is assigned to treatment.
#' @param v2 A vector of the treatment minus control differences for the pairs where the second unit is assigned to treatment.
#' @param v12a A matrix with 2 columns. The first column contains a set of imputed differences had the first unit in each pair been treated.  The first column contains a set of imputed differences had the second unit in each pair been treated.
#' @param v12b A second set of imputed potential differences.
#' @export

p_interp = function(Tr,v1,v2,v12a,v12b){
  
  #two imputation methods (a and b) for v1
  v1a = v12a[,1][Tr == 1]
  v1b = v12b[,1][Tr == 1]
  
  #two imputation methods (a and b) for v2
  v2a = v12a[,2][Tr == 0]
  v2b = v12b[,2][Tr == 0]
  
  
  #find the value for the numerator and denominator for each pair
  #there is no observed v2 for Tr == 0, so set to 0 and vs. versa
  num_v1 = num_v2 = denom_v1 = denom_v2 = rep(0,length(Tr))
  
  num_v1[Tr == 1] = (v1-v1b)*(v1a-v1b)
  denom_v1[Tr == 1] = (v1a-v1b)^2
  num_v2[Tr == 0] = (v2-v2b)*(v2a-v2b)
  denom_v2[Tr == 0] = (v2a-v2b)^2
  
  #sum over all units
  nsum_v1 = sum(num_v1)
  dsum_v1 = sum(denom_v1)
  nsum_v2 = sum(num_v2)
  dsum_v2 = sum(denom_v2)
  
  #calculate alpha for v1 and v2
  alpha_v1 = (nsum_v1 - num_v1)/(dsum_v1 - denom_v1)
  alpha_v2 = (nsum_v2 - num_v2)/(dsum_v2 - denom_v2)
  alpha_v12 = cbind(alpha_v1,alpha_v2)
  
  #restrinct to the interval [0,1]
  alpha_v12 = ifelse(alpha_v12 < 0,0,ifelse(alpha_v12 > 1, 1, alpha_v12))
  
  #take weighted average of two interpolations
  v12 = v12a*alpha_v12 + v12b*(1-alpha_v12)
  out = list(v12=data.frame(v12), alpha =alpha_v12)
  
  return(out)
}
