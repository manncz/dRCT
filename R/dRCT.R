#' Design-Based Treatment Effect Estimation for Randomized Experiment
#'
#' XX DESCRIPTION
#' 


#' @export
#'
design_estmtr <- function(x,...){
  UseMethod("design_estmtr")
}

#' @export
#' @importFrom rlang f_rhs f_lhs
design_estmtr.formula <- function(formula, covariates, data, 
                                  interp_method = NULL, aux_pred = NULL, 
                                  prob = NULL, block_ids = NULL,
                                  cluster_ids = NULL, weights = NULL, ...){
  
  
  # design inspiration from estimatr package for formula specifications
  
  if(missing(data)) stop("Formula provided without data.")
  
  if(length(f_rhs(formula)) > 1) stop("Formula should be of the form ",format(outcome ~ treatment)," including only two variables, the treatment indicator and the response/outcome")
  
  mf = model.frame(formula = formula, data = data)
  
  arg_list <- list(
    Y = model.response(mf),
    Tr = mf[, as.character(f_rhs(formula))],
    Z = NULL,
    aux_pred = NULL,
    block_ids = NULL,
    cluster_ids = NULL,
    weights = NULL
  )
  
  
  if(!missing(covariates)){
    stopifnot("`covariates` argument must be a formula"= inherits(covariates, "formula"))
    if(!is.null(f_lhs(covariates))) stop("Formula for covariates should include covariates on the right hand side only, for example:",format( ~ x1 + x2 + x3))
                                        
    arg_list[["Z"]] <- model.frame(formula = covariates, data = data)
  }
  
  if(!missing(aux_pred))    arg_list[["aux_pred"]]    <- eval(substitute(aux_pred), envir = data)
  if(!missing(block_ids))   arg_list[["block_ids"]]   <- eval(substitute(block_ids), envir = data)
  if(!missing(cluster_ids)) arg_list[["cluster_ids"]] <- eval(substitute(cluster_ids), envir = data)
  if(!missing(weights))     arg_list[["weights"]]     <- eval(substitute(weights), envir = data)

  
  design_estmtr.default(Y = arg_list[["Y"]],
                        Tr = arg_list[["Tr"]],
                        Z = arg_list[["Z"]],
                        interp_method = interp_method,
                        prob = prob,
                        block_ids = arg_list[["block_ids"]],
                        cluster_ids = arg_list[["cluster_ids"]],
                        weights = arg_list[["weights"]])
  
}
  
  
#' @export
#'
design_estmtr.default <- function(Y, Tr, Z = NULL, aux_pred = NULL, 
                                  interp_method = NULL, 
                                  prob = NULL, block_ids = NULL,
                                  cluster_ids = NULL, weights = NULL, ...){

  
  if(!is.null(cluster_ids) & !is.null(weights)) stop("Either cluster ids or weights can be provided (not both).")
  
  # aggregate data if clusters given
  if(!is.null(cluster_ids)){
    
    if(!is.null(block_ids)){
      
      agg_check <- data.frame(block_ids, cluster_ids) |> 
        group_by(cluster_ids) |> 
        mutate(n_block = n_distinct(block_ids)) |> 
        filter(n_block > 1) |> 
        nrow()
      
      if(agg_check != 0) stop("Units in a cluster must all be in the same block.")
      
      X <- Z
      if(is.null(Z)) X <- rep(1, length(Y))

      agg.dat <- data.frame(Y, Tr, X, block_ids, cluster_ids) |> 
        group_by(block_ids, cluster_ids) |> 
        summarize(across(everything(), mean),
                  n = n()) |> 
        ungroup()
       
      block_ids <- agg.dat$block_ids
      cluster_ids <- agg.dat$cluster_ids
      Y <- agg.dat$Y
      Tr <- agg.dat$Tr
      weights <- agg.dat$n
      
      if(!is.null(Z)){
        Z <- agg.dat |> 
          select(-Y, -Tr, -block_ids, -cluster_ids, -n) |> 
          as.matrix()
      }
      
    }else  stop("Cluster randomization designs other than paired cluster randomization are not yet implemented.")

  }
  
  
  if(!is.null(block_ids)){
    
    if(length(block_ids) != length(Y)) stop("Block ids must be the same length as outcome.")
    
    if(length(block_ids)/length(table(block_ids)) == 2) {
      design <- "paired"
    }else{
      design <- "blocked"
      stop("Analyzing block randomized trials is not supported at this time.")
    }
    
    if(!is.null(prob)) warning("Provided treatment probability ignored - assuming complete randmoization within blocks/pairs.")
    
  }else if(!is.null(prob)){
    design <- "bernoulli"
  }else stop("Complete randomization is not yet implemented. Provide a treatment probability `prob` if intending Bernoilli randomization.")
  
  
  if(design == "bernoulli"){
    if(is.null(aux_pred)){
      
      if(is.null(interp_method)){
        interp_method <- loop_rf
        message("No interpolation method provided, so defaulting to Random Forest imputation")
      }else{
        interp_input <- as.character(substitute(interp_method))
        #if(!interp_input %in% c("loop_glm","loop_ols","loop_rf", "loop_mean")) stop(interp_input, " is not an interpolation method for estimation with a Bernoilli randomized experiment. Choose between `loop_ols`, `loop_glm`, `loop_rf`, and `loop_mean`")
      }
      
      out <- loop(Y, Tr, Z, pred = interp_method, p = prob, ...)
      return(out)
      
    }else{
      
      if(!is.null(interp_method)) warning("Auxiliary predictions provided so `reloop` method is used, which does not take an interpolation method.")
      
      out <- reloop(Y, Tr, Z, aux_pred)
      return(out)
      
    }
  }else if(design == "paired"){
    if(is.null(aux_pred)){
      
      if(is.null(interp_method)){
        interp_method <- p_rf_interp
        message("No interpolation method provided, so defaulting to Random Forest imputation (`p_rf_interp`)")
      }else{
        interp_input <- as.character(substitute(interp_method))
        # if(!interp_input %in% c("p_rf_interp", "p_rf_po", "p_rf_v12", 
        #                         "p_ols_interp", "p_ols_po","p_ols_v12",
        #                         "p_loomi")) stop(interp_input, " is not an interpolation method for estimation with a pair randomized experiment. Choose between appropriate methods")
      }
      
      out <- p_loop(Y, Tr, Z, pred = interp_method, P = block_ids, n = weights, ...)
      return(out)
      
    }else stop("Reloop is not yet implemented for pair randomized experiments. You can include the auxiliary predictions in the covariate matrix.")
  }
  
#end function
}
  
  

