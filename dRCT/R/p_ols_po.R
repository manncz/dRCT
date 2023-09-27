#' OLS Imputation of Potential Outcomes in Paired Experiments
#'
#' This function is used to impute potential outcomes for the \code{p_loop} function.
#' Leaves out each pair and imputes its potential outcomes using ordinary least squares
#' regression on the remaining pairs. Individuals are treated as units when making predictions.
#' @param ordered A matrix of pair experimental data that his been processed by the \code{pair} function, with the treatment pair first.
#' @param assigned A matrix of pair experimental data that his been processed by the \code{pair} function.
#' @param n_assigned A matrix of pair experimental data cluster sizes that has been processed by the \code{pair} function.
#' @export

p_ols_po = function(ordered, assigned, n_assigned){

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
  coefs1 <- coefs0 <- matrix(nrow=M,ncol=ncol(dat1))
  for(i in 1:M){
    lm1 = lm(Y ~ ., dat1[-i,]) # treatment model
    lm0 = lm(Y ~ ., dat0[-i,]) # control model

    that1 = predict(lm1,obs1[i,])
    chat1 = predict(lm0,obs1[i,])
    that2 = predict(lm1,obs2[i,])
    chat2 = predict(lm0,obs2[i,])

    v1[i] = (that1*n1[i] - chat2*n2[i])
    v2[i] = (that2*n2[i] - chat1*n1[i])

    coefs1[i,] <- coef(lm1)
    coefs0[i,] <- coef(lm0)

  }
  colnames(coefs1) <- colnames(coefs0) <- names(coef(lm1))

  out <- data.frame(v1,v2)
  attr(out,"coefs")=list(coefs1=coefs1,coefs0=coefs0)

  return(out)
}
