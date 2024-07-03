#' GLM Imputation of Potential Outcomes
#'
#' This function is used to impute potential outcomes for the \code{loop} function. 
#' Leaves out each observation and imputes its potential outcomes using a generalized linear model
#' regression on the remaining observations.
#' For use, be sure to specify family in the call to `loop()`
#' 
#' @param Y A vector of observed outcomes.
#' @param Tr The treatment assignment vector.
#' @param Z A matrix of pre-treatment covariates.
#' @param ... Other options passed to the `glm` function (for example, `family`) 
#' @export
#' @examples
#' ## Create Simulated Data
#' N = 30 # 30 observations
#' k = 5 # 5 covariates
#' 
#' Z = matrix(runif(N*k,-1,1)*10,N,k)
#' c = rbinom(N,1,plogis(Z %*% runif(k)))
#' t = rbinom(N,1,plogis(Z %*% runif(k)+0.5))
#' 
#' Tr = rbinom(N,1,0.5)
#' Y = ifelse(Tr == 0, c, t)
#' 
#' ## Run LOOP and compare with a Difference in Means
#' looprf.results = loop(Y, Tr, Z, loop_rf)
#' loopols.results = loop(Y, Tr, Z, loop_ols)
#' loopglm.results = loop(Y, Tr, Z, loop_glm,family=binomial(logit))
#' 
#' meandiff = mean(Y[Tr == 1]) - mean(Y[Tr == 0])
#' varhat = var(Y[Tr==1])/length(Y[Tr==1]) + var(Y[Tr==0])/length(Y[Tr==0])
#' 
#' print(c("loop rf",looprf.results))
#' print(c("loop ols",loopols.results))
#' print(c("loop logistic regression",loopglm.results))
#' print(c("Difference in Means",meandiff,varhat))


loop_glm <- function(Y,Tr,Z,...){
  N = length(Y)
  dat = data.frame(Y,Tr, Z)

  fullModT <- glm(Y~.-Tr,data=dat,subset=Tr==1,...)
  fullModC <- glm(Y~.-Tr,data=dat,subset=Tr==0,...)

  yhatTfull <- predict(fullModT,newdata=dat,type="response")
  yhatCfull <- predict(fullModC,newdata=dat,type="response")
  
  that <- numeric(N)
  chat <- numeric(N)
  
  for(i in 1:N){
    if(Tr[i]==1){
      chat[i] <- yhatCfull[i]
      that[i] <- predict(update(fullModT,data=dat[-i,]),newdata=dat[i,],type='response')
    } else{
      chat[i] <- predict(update(fullModC,data=dat[-i,]),newdata=dat[i,],type='response')
      that[i] <- yhatTfull[i]
    }
  }
  t_c=cbind(that,chat)
  return(t_c)
}
