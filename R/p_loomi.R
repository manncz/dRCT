#' Leave-One-Out Mean Imputation of Potential Outcomes in Paired Experiments
#'
#' This function is used to impute potential outcomes for the \code{p_loop} function. Pairs are treated
#' as the unit. This function used with \code{p_loop} implements the LOO-MI estimator when \code{loo} is set to \code{TRUE} (the default).
#' If \code{loo} is set to \code{FALSE} and  \code{weighted_imp} is set to \code{TRUE}, (or all of the cluster sizes are the same),
#' then the point estimator returned by \code{p_loop} is equivalent to the difference in means estimator.
#' @param ordered A matrix of pair experimental data that his been processed by the \code{pair} function, with the treatment pair first.
#' @param assigned A matrix of pair experimental data that his been processed by the \code{pair} function.
#' @param n_assigned A matrix of pair experimental data cluster sizes that has been processed by the \code{pair} function.
#' @param loo If set to \code{TRUE}, leave-one-pair-out imputation should be used.
#' @param weighted_imp If set to \code{TRUE}, cluster sizes will be used as weights in imputation models.
#' @export

p_loomi <- function(ordered, assigned, n_assigned, loo = TRUE, weighted_imp = FALSE){
  
  M <- nrow(assigned)
  v1 <- v2 <- rep(NA, M)
  
  n1 = n_assigned$n1
  n2 = n_assigned$n2
  
  nt <- n1*assigned$Tr + n2*(1-assigned$Tr)
  nc <- n1*(1-assigned$Tr) + n2*assigned$Tr
  
  if(!loo){
    if(weighted_imp){
      mean_sum1 = mean_sum2 <- (sum(ordered$Y1*nt)/sum(nt) + sum(ordered$Y2*nc)/sum(nc))/2
      mean_dif1 = mean_dif2 <- (sum(ordered$Y1*nt)/sum(nt) - sum(ordered$Y2*nc)/sum(nc))/2
    }else{
      mean_sum1 = mean_sum2 <- mean((ordered$Y1 + ordered$Y2)/2)
      mean_dif1 = mean_dif2 <- mean((ordered$Y1 - ordered$Y2)/2)
    }
  }
    
  for(i in 1:M){
    
    if(loo){
      if(weighted_imp){
        mean_sum1 = mean_sum2 <- (sum((ordered$Y1*nt)[-i])/sum(nt[-i]) + sum((ordered$Y2*nc)[-i])/sum(nc[-i]))/2
        mean_dif1 = mean_dif2 <- (sum((ordered$Y1*nt)[-i])/sum(nt[-i]) - sum((ordered$Y2*nc)[-i])/sum(nc[-i]))/2
      }else{
        mean_sum1 = mean_sum2 <- mean((ordered$Y1[-i] + ordered$Y2[-i])/2)
        mean_dif1 = mean_dif2 <- mean((ordered$Y1[-i] - ordered$Y2[-i])/2)
      }
    }
    
    v1[i] <- (n1[i] - n2[i])*mean_sum1 + (n1[i] + n2[i])*mean_dif1
    v2[i] <- (n2[i] - n1[i])*mean_sum2 + (n1[i] + n2[i])*mean_dif2
    
  }
  
  return(data.frame(v1,v2))
}
