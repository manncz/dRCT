
new_loopEst <- function(x=double(),contrast="Treatment vs Control",...){
  stopifnot(is.double(x))

  structure(
    x,
    contrast=contrast,
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
print.loopEst <- function(x,pAdj=NULL,pAdjMethod=NULL,digits=max(3L,getOption('digits')),
                          signif.stars = getOption("show.signif.stars"),
                          ...){
  coefs <- rbind(c(x[1],sqrt(x[2]),test.loopEst(x,...)))
  rownames(coefs) <- "ATE:"
  colnames(coefs) <- c("Estimate", "Std. Error", "t value", "Pr(>|t|)")

  if(!is.null(pAdj)){
    coefs <- cbind(coefs,pAdj)
    colnames(coefs)[5] <- paste0("Adjusted p (",pAdjMethod,")")
  }

  print(attributes(x)$contrast)
  stats::printCoefmat(coefs,digits=digits, signif.stars = signif.stars)

  invisible(x)
}

#' @export
is.loopEst <- function(x) inherits(x,"loopEst")


#' Confidence interval for LOOP estimate
#'
#' @param x LOOP estimate
#' @export
confint.loopEst <- function(x,level=0.95,alpha=1-level,...){
  est <- x[1]
  se <- sqrt(x[2])
  df <- attributes(x)$df
  if(is.null(df)) df <- Inf
  mult <- -qt(alpha/2,df=df)
  out <- rbind(est+c(-1,1)*mult*se)
  colnames(out) <- paste(round(c(alpha/2*100,(1-alpha/2)*100),3),"%")
  rownames(out) <- "ATE"
  out
}

#' Two-sided hypotheis test for LOOP estimate
#'
#' @param x LOOP estimate
#' @param mu0 null value for test (default=0)
#' @export
test.loopEst <- function(x,mu0=0,...){
  statistic <- (x[1]-mu0)/sqrt(x[2])
  df <- attributes(x)$df
  if(is.null(df)) df=Inf

  pval <- 2*pt(-abs(statistic),df=df)

  out <- c(statistic=statistic,
                    pval=pval)
  attr(out,"df")=df
  attr(out,"mu0")=mu0
  out
}


new_TrVec <- function(x=integer(),labels=c("Control","Treatment"),ind=seq_along(x),...){
  stopifnot(is.integer(x))


  structure(
    x,
    labels=labels,
    ind=ind,
    ...,
    class=c("TrVec","integer","numeric")
    )
}

#' @export
validate_TrVec <- function(x){
  stopifnot(all(x%in%c(0,1)))
  x
}

#' @export
TrVec <- function(x,...){
  if(is.numeric(x)) x <- as.integer(x)

  x <- new_TrVec(x,...)

  validate_TrVec(x)
}

is.TrVec <- function(x)
  inherits(x,"TrVec")


new_loopList <- function(x=list(),...){
  stopifnot(is.list(x))

  structure(
    x,
    ...,
    class=c("loopList","list")
  )
}

validate_loopList <- function(x){
  stopifnot(inherits(x,"loopList"))
  stopifnot(all(vapply(x, inherits,what="loopEst",TRUE)))

  x
}

loopList <- function(x,...){
  x <- new_loopList(x,...)
  validate_loopList(x)
}

print.loopList <- function(x,pAdjMethod="holm",...){
  if(length(x)==1) print(x[[1]])
  else{
    pvals <- vapply(x,
                    function(x) test.loopEst(x)['pval'],
                    0.1)

    padj <- p.adjust(pvals,method=pAdjMethod)

    for(lp in seq_along(x))
      print(x[[lp]],pAdj=padj[lp],pAdjMethod=pAdjMethod,...)
  }
}
