#' P-LOOP Estimator
#'
#' Covariate-adjusted estimator of the average treatment effect in pair randomized experiments. Given a set of experimental data (including observed outcome, treatment assignment, pair assignments, and covariates), uses a method \code{pred} to impute the potential outcomes for each observation, which are then used to estimate the average treatment effect. 
#' @param Y A vector of mean experimental outcomes for each unit.
#' @param Tr The treatment assignment vector.
#' @param Z A matrix or vector of covariates.
#' @param P A vector encoding the pair assignments. Observations with the same value of \code{P} will be assumed to be in the same pair.
#' @param n A vector encoding the cluster sizes.
#' @param pred The prediction algorithm used to impute potential outcomes. By default, this is \code{p_rf_interp}, which uses random forests and interpolates between methods \code{p_rf_ab} and \code{p_rf_s} that treat the pairs as units or treats the individuals as units when making predictions. Another option is \code{p_ols_interp}, which uses linear regression. User written imputation methods may be used as well.
#' @param ... Arguments to be passed to the imputation method.
#' @import dplyr 
#' @import tidyr 
#' @import stringr
#' @export

p_loop = function(Y,Tr,Z=NULL,P,n=NULL,pred = p_simple,...){
  #if n isn't provided assume units have size 1
  if(is.null(n)){n = rep(1, length(Y))}
  
  #if covariates aren't provided, set as 1 as a placeholder
  if(is.null(Z)){
    Z = rep(1, length(Y))
    message("Note: unadjusted estimator is automatically used (pred = p_simple) because no covariates were provided")
    pred = p_simple
  }
  
  Z <- as.matrix(Z)
  colnames(Z) <- paste0("Z", 1:(ncol(Z)))
  
  # reshape data
  pair_out = pair(Y,Tr,Z,P,n)
  
  # data with treatment unit first
  ordered = pair_out$ordered
  
  # data as assigned with Tr indicator
  assigned = pair_out$agg
  n_assigned = pair_out$n_assigned
  
  # interpolate v1 and v2
  v12 = pred(ordered, assigned=assigned, n_assigned=n_assigned, ...)
  v1= v12$v1
  v2= v12$v2
  
  # treatment effect estimate
  d = v1-v2
  V = 2*(assigned$Tr*(assigned$Y1*n_assigned$n1 - assigned$Y2*n_assigned$n2) + (1-assigned$Tr)*(assigned$Y2*n_assigned$n2 - assigned$Y1*n_assigned$n1))
  tauhat = sum(V + (1-2*assigned$Tr)*d)/sum(n_assigned$n1, n_assigned$n2)
  
  # variance estimate
  varhat = p_loop_var(assigned,v1,v2,n_assigned)
  
  return(list(tauhat = tauhat,
              varhat = varhat))
}

