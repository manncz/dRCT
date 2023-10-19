#' Variance Estimate for the P-LOOP Estimator
#'
#' This function implements the variance estimate for the P-LOOP estimator and is used within \code{p_loop}.
#' @param assigned A matrix of pair experimental data that has been processed by the \code{pair} function
#' @param v1 A vector of the predicted treatment minus control differences for the pairs where the first unit is assigned to treatment.
#' @param v2 A vector of the predicted treatment minus control differences for the pairs where the second unit is assigned to treatment.
#' @param n_assigned A matrix of cluster sizes that has been processed by the \code{pair} function
#' @export

p_loop_var = function(assigned,v1,v2,n_assigned){
  n1 = n_assigned$n1
  n2 = n_assigned$n2
  
  V1 = 2*((assigned$Y1*n1)[assigned$Tr == 1] - (assigned$Y2*n2)[assigned$Tr == 1])
  V2 = 2*((assigned$Y2*n2)[assigned$Tr == 0] - (assigned$Y1*n1)[assigned$Tr == 0])
  
  S1 = sum((2*v1[assigned$Tr == 1] - V1)^2)
  S2 = sum((2*v2[assigned$Tr == 0] - V2)^2)
  
  N = sum(n1,n2)
  varhat = (S1 + S2)/N^2
  
  return(varhat)
}

