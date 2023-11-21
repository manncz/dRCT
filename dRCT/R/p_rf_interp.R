#' Random Forest Imputation of Potential Outcomes in Paired Experiments
#'
#' This function is used to impute potential outcomes for the \code{p_loop} function. 
#' Leaves out each pair and imputes its potential outcomes using random forests
#' on the remaining pairs. This function interpolates between the \code{p_rf_v12} and \code{p_rf_po} methods.
#' @param ordered A matrix of pair experimental data that his been processed by the \code{pair} function, with the treatment pair first.
#' @param assigned A matrix of pair experimental data that his been processed by the \code{pair} function.
#' @param n_assigned A matrix of pair experimental data cluster sizes that has been processed by the \code{pair} function.
#' @param reparm If set to \code{TRUE}, covariates will be parameterized as the pairwise differences and means for each covariate.
#' @export

p_rf_interp = function(ordered, assigned, n_assigned, reparm = TRUE){
  # two different predictions of v1 and v2
  pred1 = p_rf_v12(ordered, assigned, n_assigned, reparm)
  pred2 = p_rf_po(ordered, assigned, n_assigned)
  
  # observed v1 and v2
  Tr = assigned$Tr
  v = rep(0, length(Tr))
  v[Tr == 1] = (assigned$Y1*n_assigned$n1)[Tr == 1]-(assigned$Y2*n_assigned$n2)[Tr == 1]
  v[Tr == 0] = (assigned$Y2*n_assigned$n2)[Tr == 0]-(assigned$Y1*n_assigned$n1)[Tr == 0]
  
  # interpolation
  interpolation = p_interp(Tr,vobs=v,pred1,pred2)
  v12 = interpolation$v12
  
  attr(v12,"fitInfo") <- list(v12_info=attributes(pred1),
                              po_info=attributes(pred2),
                              alpha=interpolation$alpha)
  
  return(v12)
}