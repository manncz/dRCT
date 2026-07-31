#' Out-of-Bag Predictions from Ranger Random Forest Based on an ID-Based
#'
#' Computes out-of-bag (OOB) predictions from a fitted or newly trained ranger random forest,
#' with predictions stratified by an ID variable. Predictions for test units with a given ID
#' are based only on trees that did not use any training units with that same ID.
#'
#' @param train A data frame containing the training data. Must contain a column
#'   specified by the `id` argument.
#' @param test A data frame containing the test data. Must contain a column
#'   specified by the `id` argument. Can overlap with `train`.
#' @param rf An optional fitted ranger object. If provided, must have been fit with
#'   `keep.inbag = TRUE`. If `NULL`, a new RF will be trained using the remaining arguments.
#' @param form A formula specifying the model (e.g., `Y ~ X1 + X2 + X3`).
#'   Only used if `rf` is `NULL`. If missing, constructed from `outcome` and `covars`.
#' @param outcome Character string specifying the name of the outcome variable
#'   in `train`. Default is `"Y"`. Only used if `rf` is `NULL` and `form` is missing.
#' @param covars Character vector of predictor names. Default is `c("X1", "X2", "X3")`.
#'   Only used if `rf` is `NULL` and `form` is missing.
#' @param id Character string specifying the column name in both `train` and `test`
#'   that identifies groups for stratified OOB prediction. Default is `"pair"`.
#' @param ntree Integer specifying the number of trees to grow. Default is either 
#'   `500*exp(M-1)`, where `M` is the the maximum number of training observations that 
#'   are filtered out for a single prediction, or 10,000, which ever is smaller. For 
#'.  `M<=4`, this ensures that the expected number of trees used in an individual prediction 
#'  is roughly the same as standard OOB predictions from an RF with 500 trees. 
#'   Only used if `rf` is `NULL`.
#' @param nodesize Integer specifying the minimum node size. Default is 5.
#'   Only used if `rf` is `NULL`.
#' @param returnRF Logical. If `TRUE`, the fitted ranger object is attached as an
#'   attribute (`rf`) to the returned predictions. Default is `FALSE`.
#'
#' @return A numeric vector of OOB predictions with length equal to `nrow(test)`.
#'   Two attributes are attached:
#'   \describe{
#'     \item{ntreeP}{Integer vector indicating the number of trees used for each prediction
#'       (i.e., trees that were out-of-bag with respect to all training observations
#'       with the same ID value).}
#'     \item{rf}{(Only if `returnRF = TRUE`) The fitted ranger RF object.}
#'   }
#'
#' @details
#' The ID-based stratification works as follows: for each unique ID value in the test set,
#' this function identifies corresponding training observations with the same ID. A tree
#' contributes to the prediction only if it was out-of-bag for ALL of those training observations.
#' If no training observations share an ID value, all trees are used.
#'
#' This approach is useful for:
#' \itemize{
#'   \item Computing true OOB predictions when observations are grouped or clustered
#'   \item Cross-validation schemes where avoiding data leakage within groups matters
#'   \item Stratified OOB evaluation when train and test sets overlap
#' }
#'
#' If `rf = NULL`, a new ranger RF is trained on `train` with `keep.inbag = TRUE`.
#' The `train` and `test` datasets can overlap; the ID column controls OOB behavior.
#'
#' @examples
#' \dontrun{
#' # Example 1: Standard OOB predictions (should match ranger output exactly)
#' RF <- ranger::ranger(mpg ~ cyl + disp + hp + drat + wt + qsec,
#'                      data = mtcars,
#'                      keep.inbag = TRUE)
#' mtcars$id <- 1:nrow(mtcars)
#' OOB <- rangerOOB(mtcars, mtcars, rf = RF, id = "id")
#' all.equal(OOB, RF$predictions, check.attributes = FALSE)
#'
#' # Example 2: Stratified predictions by car make
#' mtcars$make <- substr(rownames(mtcars), 1, 2)
#' OOB_by_make <- rangerOOB(mtcars, mtcars, rf = RF, id = "make")
#' }
#'
#' @importFrom ranger ranger
#' @export
rangerOOB <- function(train,
                      test,
                      rf=NULL,
                      form,
                      outcome="Y",
                      covars=c("X1","X2","X3"),
                      id="pair",
                      ntree,nodesize=5,
                      returnRF=FALSE){
    stopifnot(require(ranger))
    Ptest <- test[[id]]
    Ptrain <- train[[id]]
    
    if(is.null(rf)){
        if(missing(form))
            form <- paste(outcome,"~",paste(covars,collapse="+"))|>
                as.formula()
        
        if(missing(ntree)){
            Pdrop <- Ptrain[Ptrain%in%Ptest]
            ntree <- min(10000,500*exp(max(table(Pdrop))-1))
        }
        
        
        rf <- ranger(form,
                     data=train,
                     num.trees=ntree,min.node.size=nodesize,
                     keep.inbag=TRUE)
    }
    pred <- predict(rf,data=test,predict.all=TRUE)
    inbag <- do.call("cbind",rf$inbag.counts)
    oob <- inbag==0
    mu <- ntreeP <- numeric(nrow(test))

    for(pp in unique(Ptest)){
        testrows <- which(Ptest==pp)
        trainrows <- which(Ptrain==pp)
        use <- if(length(trainrows)){
                   apply(oob[trainrows,,drop=FALSE],2,all)
               } else
                   rep(TRUE,ntree)
        ntreeP[testrows] <- sum(use)
        mu[testrows] <- rowMeans(pred$predictions[testrows,use,drop=FALSE])
    }
    out <- mu
    attr(out,"ntreeP") <- ntreeP
    if(returnRF) attr(out,"rf") <- rf
    out
}
