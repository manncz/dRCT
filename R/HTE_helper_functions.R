#' Extract unbiased individual treatment effect estimates from a LOOP object
#' 
#' @param obj LOOP estimate, output from the `loop()` function
#' @export
getITE=function(obj){
  stopifnot(is.loopEst(obj))
  attributes(obj)$ITE
}

