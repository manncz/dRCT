#' Interpolates between Sets of Potential Differences in Paired Experiments
#'
#' This function is used within the \code{p_rf_interp} and \code{p_ols_interp} imputation methods for the \code{p_loop} function. 
#' Given two sets of imputed potential differences, this function returns a single set of predicted differences by picking weights
#' that minimizes the mean squared error between the observed differences and the weighted average imputed differences. This minimization
#' is done using leave-one-out cross validation.
#' 
#' @param Tr The treatment assignment vector for the pairs, where the value is 1 if the first unit in the pair is assigned to treatment and 0 otherwise.
#' @param vobs A vector of the treatment minus control differences for all pairs.
#' @param v12a A matrix with 2 columns. The first column contains a set of imputed differences had the first unit in each pair been treated.  The second column contains a set of imputed differences had the second unit in each pair been treated.
#' @param v12b A second set of imputed potential differences.
#' @export

p_interp = function(Tr,vobs,v12a,v12b){
  
  va = vb = rep(0, length(Tr))
  
  #two imputation methods (a and b) for v1 when the first unit is treated
  va[Tr == 1] = v12a[,1][Tr == 1]
  vb[Tr == 1] = v12b[,1][Tr == 1]
  
  #two imputation methods (a and b) for v2 when the second unit was treated
  va[Tr == 0] = v12a[,2][Tr == 0]
  vb[Tr == 0] = v12b[,2][Tr == 0]
  
  #find the value for the numerator and denominator for each pair
  num = (vobs-vb)*(va-vb)
  denom = (va-vb)^2
  
  #sum over all units
  nsum = sum(num)
  dsum = sum(denom)
  
  #calculate alpha
  alpha = (nsum - num)/(dsum - denom)
  
  #restrinct to the interval [0,1]
  alpha = ifelse(alpha < 0,0,ifelse(alpha > 1, 1, alpha))
  alpha_i = cbind(alpha, alpha)
  
  #take weighted average of two interpolations
  v12 = v12a*alpha_i + v12b*(1-alpha_i)
  out = list(v12=data.frame(v12), alpha =alpha)
  
  return(out)
}
