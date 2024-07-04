#' Leave-One-Out Mean Imputation of Potential Outcomes
#'
#' This function imputes subject i's potential outcomes using the means of
#' other subjects' observed outcomes.
#' The resulting ATE estimate is equal to the Neyman difference-in-means
#' estimate; hence, this function is intended for illustrative purposes.
#'
#' @param Y A vector of observed outcomes.
#' @param Tr The treatment assignment vector.
#' @param Z A matrix of pre-treatment covariates. Ignored.
#' @param ... ignored.
#'
#' @export

loop_mean <- function(Y,Tr,Z,...){

  n <- length(Tr)
  nT <- sum(Tr)
  nC <- n-nT
  sumT <- sum(Y[Tr==1])
  sumC <- sum(Y[Tr==0])

  that <- chat <- numeric(n)
  that[Tr==0] <- sumT/nT
  chat[Tr==1] <- sumC/nC

  that[Tr==1] <- (sumT-Y[Tr==1])/(nT-1)

  chat[Tr==0] <- (sumC-Y[Tr==0])/(nC-1)

  cbind(that,chat)
}
