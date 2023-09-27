
new_loopEst <- function(x=double(),...){
  stopifnot(is.double(x))

  structure(
    x,
    ...,
    class=c("loopEst","double","numeric")
    )
}

#' @export
validate_loopEst <- function(x){
  if(length(x)!=2)
    stop("loop estimate should only contain two values, an ATE and variance estimate")

  if(x[2]<0) stop("The sampling variance is negative!")

  x
}

#' @export
loopEst <- function(x,...){
  if(is.numeric(x)) x <- as.double(x)
  if(is.list(x)) x <- as.double(unlist(x))

  x <- new_loopEst(x,...)

  validate_loopEst(x)
}

#' Print loopEst object
#'
#' @param x loopEst object
#' @return x
#' @export
print.loopEst <- function(x,digits=max(3L,getOption('digits')),...){
  df <- data.frame(ATE=x[1],SE=sqrt(x[2]))
  rownames(df) <- "Estimate:"

  print.data.frame(df,digits=digits,...)
  invisible(x)
}

#' @export

is.loopEst <- function(x) inherits(x,"loopEstxo")


#' Confidence interval for LOOP estimate
#'
#' @param x LOOP estimate
#' @export
confint.loopEst <- function(x,level=0.95,alpha=1-level,df=Inf,...){
  est <- x[1]
  se <- sqrt(x[2])
  mult <- -qt(alpha/2,df=df)
  out <- rbind(est+c(-1,1)*mult*se)
  colnames(out) <- paste(c(alpha/2*100,(1-alpha/2)*100),"%")
  rownames(out) <- "ATE"
  out
}
