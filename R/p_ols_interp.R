#' OLS Imputation of Potential Outcomes in Paired Experiments (Interpolated)
#'
#' This function is used to impute potential outcomes for the \code{p_loop} function.
#' Leaves out each pair and imputes its potential outcomes using ordinary least squares
#' regression on the remaining pairs. This function interpolates between the \code{p_ols_v12} and \code{p_ols_po} methods.
#' @param ordered A matrix of pair experimental data that his been processed by the \code{pair} function, with the treatment pair first.
#' @param assigned A matrix of pair experimental data that his been processed by the \code{pair} function.
#' @param n_assigned A matrix of pair experimental data cluster sizes that has been processed by the \code{pair} function.
#' @param weighted_imp If set to \code{TRUE}, cluster sizes will be used as weights in imputation models.
#' @export

p_ols_interp = function(ordered, assigned, n_assigned, weighted_imp = F){
  # two different predictions of v1 and v2
  pred1 = p_ols_v12(ordered, assigned, n_assigned, weighted_imp)
  pred2 = p_ols_po(ordered, assigned, n_assigned, weighted_imp)

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
