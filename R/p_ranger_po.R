#' Ranger-Based Random Forest Imputation of Potential Outcomes in Paired Experiments (Individuals as Units)
#'
#' This function is used to impute potential outcomes for the \code{p_loop} function.
#' It functions similarly to p_rf_po(), but instead of fitting separate random forests
#' leaving out each pair with the `randomForest` package, it fits one random forest to
#' each treatment arm, using the `ranger` package, but ensures that the out-of-box predictions
#' use trees that excluded both members of each pair.
#' @param ordered A matrix of pair experimental data that his been processed by the \code{pair} function, with the treatment pair first.
#' @param assigned A matrix of pair experimental data that his been processed by the \code{pair} function.
#' @param n_assigned A matrix of pair experimental data cluster sizes that has been processed by the \code{pair} function.
#' @param weighted_imp If set to \code{TRUE}, cluster sizes will be used as weights in imputation models.
#' @param ntree Integer specifying the number of trees to grow. Default is 500
#' @param nodesize Integer specifying the minimum node size. Default is 5.
#'
#' @importFrom ranger ranger
#' @export
p_ranger_po <- function(ordered, assigned, n_assigned, weighted_imp=FALSE,
                        ntree=500,nodesize=5,...){

    stopifnot(require(ranger))

    ## reconstruct the dataset
    stopifnot(all(ordered$P==assigned$P))
    dat1 <- ordered[,endsWith(names(ordered),"1")]
    names(dat1) <- substr(names(dat1),1,nchar(names(dat1))-1)
    dat2 <- ordered[,endsWith(names(ordered),"2")]
    names(dat2) <- substr(names(dat2),1,nchar(names(dat2))-1)

    obs1 <- assigned[,endsWith(names(assigned),"1")]
    names(obs1) <- substr(names(obs1),1,nchar(names(obs1))-1)
    obs2 <- assigned[,endsWith(names(assigned),"2")]
    names(obs2) <- substr(names(obs2),1,nchar(names(obs2))-1)


    n1 = n_assigned$n1
    n2 = n_assigned$n2
    M = nrow(ordered)
    v1 = v2 = that1=that2=chat1=chat2=rep(0, M)


    rfT <- ranger::ranger(Y ~ . - n, dat1,
                  num.trees=ntree,min.node.size=nodesize,
                     keep.inbag=TRUE, case.weights = dat1$n,...)

    rfC <- ranger::ranger(Y ~ . - n, dat2,
                  num.trees=ntree,min.node.size=nodesize,
                  keep.inbag=TRUE, case.weights = dat2$n,...)


    that1 <- rangerOOB(
        train=cbind(dat1,P=ordered$P),
        test=cbind(obs1,P=assigned$P),
        rf=rfT,id="P")
    that2 <- rangerOOB(
        train=cbind(dat1,P=ordered$P),
        test=cbind(obs2,P=assigned$P),
        rf=rfT,id="P")
    chat1 <- rangerOOB(
        train=cbind(dat2,P=ordered$P),
        test=cbind(obs1,P=assigned$P),
        rf=rfC,id="P")
    chat2 <- rangerOOB(
        train=cbind(dat2,P=ordered$P),
        test=cbind(obs2,P=assigned$P),
        rf=rfC,id="P")


    npair <- nrow(assigned)

#    return(data.frame(that1,chat1,that2,chat2))
    v1 = (that1*ordered$n1 - chat2*ordered$n2)
    v2 = (that2*ordered$n2 - chat1*ordered$n1)

    out <- data.frame(v1,v2)

    attr(out,"forests") <- forests <- list(rfT=rfT,rfC=rfC)

    return(out)
}


#' Ranger-Based Combined Random Forest Imputation of Potential Outcomes in Paired Experiments (Individuals as Units)
#'
#' This function is used to impute potential outcomes for the \code{p_loop} function.
#' It functions similarly to p_rf_po(), but instead of fitting separate random forests
#' leaving out each pair with the `randomForest` package, it fits one random forest,
#' including the treatment indicator as a predictor and using the `ranger` package,
#' but ensures that the out-of-box predictions use trees that excluded both members
#' of each pair.
#' @param ordered A matrix of pair experimental data that his been processed by the \code{pair} function, with the treatment pair first.
#' @param assigned A matrix of pair experimental data that his been processed by the \code{pair} function.
#' @param n_assigned A matrix of pair experimental data cluster sizes that has been processed by the \code{pair} function.
#' @param weighted_imp If set to \code{TRUE}, cluster sizes will be used as weights in imputation models.
#' @param ntree Integer specifying the number of trees to grow. Default is `500*exp(1)`
#' @param nodesize Integer specifying the minimum node size. Default is 5.
#'
#' @importFrom ranger ranger
#' @export
p_ranger_po_OneRF <- function(ordered, assigned, n_assigned, weighted_imp=FALSE,
                              ntree=500*exp(1),nodesize=5,...){

    stopifnot(require(ranger))

    ## reconstruct the dataset
    stopifnot(all(ordered$P==assigned$P))
    dat1 <- ordered[,endsWith(names(ordered),"1")]
    names(dat1) <- substr(names(dat1),1,nchar(names(dat1))-1)
    dat1$Tr <- 1
    dat2 <- ordered[,endsWith(names(ordered),"2")]
    names(dat2) <- substr(names(dat2),1,nchar(names(dat2))-1)
    dat2$Tr <- 0

    dat <- rbind(dat1,dat2)

    obs1 <- assigned[,endsWith(names(assigned),"1")]
    names(obs1) <- substr(names(obs1),1,nchar(names(obs1))-1)
    obs2 <- assigned[,endsWith(names(assigned),"2")]
    names(obs2) <- substr(names(obs2),1,nchar(names(obs2))-1)


    n1 = n_assigned$n1
    n2 = n_assigned$n2
    M = nrow(ordered)
    v1 = v2 = that1=that2=chat1=chat2=rep(0, M)


    rf <- ranger(Y ~ . - n, dat,
                  num.trees=ntree,min.node.size=nodesize,
                     keep.inbag=TRUE, case.weights = dat$n,...)


    that1 <- rangerOOB(
        train=cbind(dat,P=rep(ordered$P,2)),
        test=within(obs1,{
            P <- assigned$P
            Tr <- 1}),
        rf=rf,id="P")


    that1 <- rangerOOB(
        train=cbind(dat,P=rep(ordered$P,2)),
        test=within(obs2,{
            P <- assigned$P
            Tr <- 1}),
        rf=rf,id="P")


    chat1 <- rangerOOB(
        train=cbind(dat,P=rep(ordered$P,2)),
        test=within(obs1,{
            P <- assigned$P
            Tr <- 0}),
        rf=rf,id="P")


    chat2 <- rangerOOB(
        train=cbind(dat,P=rep(ordered$P,2)),
        test=within(obs2,{
            P <- assigned$P
            Tr <- 0}),
        rf=rf,id="P")

    npair <- nrow(assigned)

#    return(data.frame(that1,chat1,that2,chat2))
    v1 = (that1*ordered$n1 - chat2*ordered$n2)
    v2 = (that2*ordered$n2 - chat1*ordered$n1)

    return(data.frame(v1,v2))
}
