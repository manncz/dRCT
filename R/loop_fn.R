#' LOOP Estimator
#'
#' Covariate-adjusted estimator of the average treatment effect in randomized experiments. Given a set of experimental data (including observed outcome, treatment assignment, and covariates), uses a method \code{pred} to impute the potential outcomes for each observation, which are then used to estimate the average treatment effect.
#' @param Y A vector of experimental outcomes.
#' @param Tr The treatment assignment vector.
#' @param Z A matrix or vector of covariates.
#' @param pred The prediction algorithm used to impute potential outcomes. By default, this is \code{loop_rf}, which uses random forests. Other options include \code{loop_ols} (which uses linear regression) and \code{reloop} (which can incorporate external predictions \code{yhat} to improve precision). As implemented, \code{reloop} is very slightly biased -- if bias is a concern, \code{reloop_slow} can be used instead. User written imputation methods may be used as well.
#' @param p The probability of being assigned to treatment. Defaults to 0.5.
#' @param returnFitInfo If set to TRUE, includes results or fitted model objects from the algorithm specified in \code{pred}, as an attribute called \code{fitInfo}
#' @param ... Arguments to be passed to the imputation method. For example, \code{loop_rf} takes the argument \code{dropobs}: when set to TRUE, lowers the bootstrap sample size in the random forest by 1 when making out-of-bag predictions. By default, this is set to TRUE if the sample size is less than or equal to 30 and FALSE otherwise. For \code{reloop} and \code{reloop_slow}, external predictions \code{yhat} must be specified. This is a vector of predictions of the outcome for each participant that is obtained using an external data source.
#' @return An object of the class \code{loopEst} which is a numeric vector with two components: the estimated average treatment effect and the estimated sampling variance. It also includes the following attributes:
#'
#' * \code{df} The estimated degrees of freedom for confidence intervals and hypothesis tests
#' * \code{call} The function call
#' * \code{pred} The prediction algorithm used
#' * \code{p} The probability of being assigned to the treatment condition
#' * \code{ITE} A vector of unbiased individual treatment effect estimates
#' * If \code{returnFitInfo=TRUE} information from the \code{pred} algorithm
#'
#' @importFrom stats predict
#' @export
#' @examples
#' ## Create Simulated Data
#' N = 30 # 30 observations
#' k = 5 # 5 covariates
#'
#' Z = matrix(runif(N*k)*10,N,k)
#' c = Z %*% runif(k) + runif(N)*5
#' t = ifelse(rowMeans(Z) < 3, c + runif(N)*5, c + runif(N)*10)
#'
#' Tr = rbinom(N,1,0.5)
#' Y = ifelse(Tr == 0, c, t)
#'
#' ## Run LOOP and compare with a Difference in Means
#' looprf.results = loop(Y, Tr, Z, loop_rf)
#' loopols.results = loop(Y, Tr, Z, loop_ols)
#'
#' meandiff = mean(Y[Tr == 1]) - mean(Y[Tr == 0])
#' varhat = var(Y[Tr==1])/length(Y[Tr==1]) + var(Y[Tr==0])/length(Y[Tr==0])
#'
#' print(c("loop rf",looprf.results))
#' print(c("loop ols",loopols.results))
#' print(c("Difference in Means",meandiff,varhat))
#'
#'

loop.default = function(Y, Tr, Z=NULL,pred = loop_rf, p = 0.5, returnFitInfo=FALSE,
               P,data,...) {

  if(is.numeric(Tr) & NCOL(Tr)==1 & all(Tr==1 | Tr==0))
    return(loop.fit(Y=Y,Tr=TrVec(Tr),Z=Z,pred=pred,p=p,returnFitInfo=returnFitInfo,...))

  TrList <- format_Tr(Tr)

  if(length(TrList)==1)
    return(loop.fit(Y=Y,Tr=TrList[[1]],Z=Z,pred=pred,p=p,returnFitInfo=returnFitInfo,...))

  out <- loopList(
    lapply(TrList,
           function(tr){
             ind <- attributes(tr)$ind
             loop.fit(Y=Y[ind],Tr=tr,Z=Z[ind,],pred=pred,p=p,returnFitInfo=returnFitInfo)
           })
  )
}

#' @export
#'
loop <- function(x,...){
  UseMethod("loop")
}

loop.formula <- function(formula,covariates,pair,data,...){
  if(missing(data)) stop("Formula provided without data.")

  mf=model.frame(formula=formula,data=data,na.action=na.fail)
  if(ncol(mf)>2) stop("Formula should be of the form ",format(outcome~treatment)," including only two variables, the treatment indicator and the response/outcome")

  Y=model.response(mf)
  Tr=mf[,2]

  if(missing(covariates)) covariates <- NULL
  if(missing(pair)) pair <- NULL

  if(!is.null(covariates)){
    if(inherits(covariates, "formula")){
      covariates <- model.frame(covariates,data)
    } else if(is.vector(covariates)){
      covariates <- data[,covariates]
    }
  } else if(is.matrix(covariates)|is.data.frame(covariates)){
    stopifnot(nrow(covariates)==nrow(data))
  }

  if(!is.null(pair)){
    if(inherits(pair, "formula")){
      if(length(pair)>2) stop("Pair formula must be one-sided with a single pair-identifier, e.g. ", format(~P))
      pair <- model.frame(pair,data)
    } else if(length(pair)==1 & (is.integer(pair)|is.character(pair))){
      pair <- data[,pair]
    }
  }

  loop.default(Y=Y,Tr=Tr,Z=covariates,P=pair)
}


loop.fit = function(Y, Tr, Z=NULL,pred = loop_rf, p = 0.5, returnFitInfo=FALSE, ...) {
  Y = as.matrix(Y)
  if(is.null(Z)){
    t_c = loop_mean(Y,Tr,...)
  } else{
    if(is.data.frame(Z)) Z=model.matrix(~.,data=Z)[,-1]
    Z=as.matrix(Z)
    t_c = pred(Y,Tr,Z,...)
  }
  that = t_c[,1]
  chat = t_c[,2]

  mhat = (1-p)*that + p*chat
  ITE = (1/p)*(Y-mhat)*Tr-(1/(1-p))*(Y-mhat)*(1-Tr)
  tauhat = mean(ITE)

  n_t = length(mhat[Tr == 1])
  n_c = length(mhat[Tr == 0])
  M_t = sum((that[Tr==1]-Y[Tr==1])^2)/n_t
  M_c = sum((chat[Tr==0]-Y[Tr==0])^2)/n_c
  varhat = 1/(n_t+n_c)*((1-p)/p*M_t + p/(1-p)*M_c + 2*sqrt(M_t*M_c))

  ## this is what t.test() uses (more or less)
  df=varhat^2/(M_t^2/n_t^2/(n_t-1)+M_c^2/n_c^2/(n_c-1))

  out <- new_loopEst(
    c(tauhat, varhat),
    contrast=makeContrast(Tr),
    df=df,
    call=match.call(),
    pred=pred,
    p=p,
    ITE=ITE)

  if(returnFitInfo)
    attr(out,"fitInfo") <- attributes(t_c)

  return(out)
}

Tr_2levs <- function(Tr,ind=seq_along(Tr)){
  if(is.factor(Tr)){
      ctl <- min(levels(Tr))
      trt <- max(levels(Tr))
    } else{
      ctl <- min(Tr)
      trt <- max(Tr)
    }
    if(ctl==0 & trt==1)  return(TrVec(Tr,labels=c(ctl,trt),ind=ind))
    return(TrVec(ifelse(Tr==ctl,0,1),labels=c(ctl,trt),ind=ind))
  }


format_Tr <- function(Tr){
  if(any(is.na(Tr))) stop("NA values not allowed in Tr")

  if(!is.numeric(Tr) & !is.factor(Tr) & !is.character(Tr))
    stop("Tr must be numeric, factor, or logical")

  if(is.matrix(Tr)){
    if(ncol(Tr)==1) Tr <- Tr[,1]
    else stop("Tr must be a vector")
  }

  levs <- unique(Tr)

  lu <- length(levs)
  if(lu==1) stop("Need at least 2 randomized conditions")
  if(lu==2) return(list(Tr_2levs(Tr)))

  if(lu> length(Tr)/2) stop("Too many randomized conditions relative to sample size")

  out <- list()
  for(cond1 in 1:(lu-1))
    for(cond2 in (cond1+1):lu){
      l1 <- levs[cond1]
      l2 <- levs[cond2]
      ind <- which(Tr==l1 | Tr==l2)
      out[[paste0(l1,"_",l2)]] <- Tr_2levs(Tr[ind],ind=ind)
    }
  out
}

makeContrast <- function(Tr){
  if(!is.TrVec(Tr)) return("Treatment vs Control")
  labs <- attributes(Tr)$labels

  paste(labs[2],"vs",labs[1])
}
