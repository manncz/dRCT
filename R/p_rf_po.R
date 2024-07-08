#' Random Forest Imputation of Potential Outcomes in Paired Experiments (Individuals as Units)
#'
#' This function is used to impute potential outcomes for the \code{p_loop} function. 
#' Leaves out each pair and imputes its potential outcomes using random forests
#' on the remaining pairs. Individuals are treated as units when making predictions.
#' @param ordered A matrix of pair experimental data that his been processed by the \code{pair} function, with the treatment pair first.
#' @param assigned A matrix of pair experimental data that his been processed by the \code{pair} function.
#' @param n_assigned A matrix of pair experimental data cluster sizes that has been processed by the \code{pair} function.
#' @param weighted_imp If set to \code{TRUE}, cluster sizes will be used as weights in imputation models.
#' @export

p_rf_po = function(ordered, assigned, n_assigned, weighted_imp = F){
  
  dat1 = ordered %>% select(ends_with("1"))
  dat0 = ordered %>% select(ends_with("2"))
  
  obs1 = assigned %>% select(ends_with("1"))
  obs2 = assigned %>% select(ends_with("2"))
  
  # make variable names the same between treatment and control data
  varnames = str_sub(colnames(dat1), end = str_length(colnames(dat1))-1)
  colnames(dat1) = colnames(dat0) = colnames(obs1) = colnames(obs2) = varnames
  
  n1 = n_assigned$n1
  n2 = n_assigned$n2
  
  M = nrow(ordered)
  v1 = v2 = rep(0,M)
  
  for(i in 1:M){
    if(weighted_imp){
      rf1 = randomForest::randomForest(Y ~ . -n, dat1[-i,], weights = dat1[-i,]$n) # treatment forest
      rf0 = randomForest::randomForest(Y ~ . -n, dat0[-i,], weights = dat0[-i,]$n) # control forest
    }else{
      rf1 = randomForest::randomForest(Y ~ . -n, dat1[-i,]) # treatment forest
      rf0 = randomForest::randomForest(Y ~ . -n, dat0[-i,]) # control forest
    }
    
    that1 = predict(rf1,obs1[i,])
    chat1 = predict(rf0,obs1[i,])
    that2 = predict(rf1,obs2[i,])
    chat2 = predict(rf0,obs2[i,])
    
    v1[i] = (that1*n1[i] - chat2*n2[i])
    v2[i] = (that2*n2[i] - chat1*n1[i])
    
  }
  return(data.frame(v1,v2))
}

