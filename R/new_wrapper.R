
#' wrapper_formula is a function that takes a formula and df and creates a list to be used in loop or p_loop
#' wrapper_vectors is a function that takes vectors and creates a list to be used in loop or p_loop
#' @param formula1 A formula of the form Y ~ Tr.
#' @param df A data frame.
#' @param Y An outcome vector.
#' @param formula2 A formula of covariates that gets turned into the matrix of covariates Z.
#' @param Tr A vector of assigned treatments.
#' @param Z A matrix of covariates.
#' @param clust A vector encoding the pair assignments.
#' @param n A vector encoding the cluster sizes if there are pairs of clusters rather than individuals.
#' @param block A vector encoding the block assignments.
#' @return Function should send the given data into either p_loop or loop functions depending on type of data provided (Bernoulli, paired, or cluster).



wrapper_formula <- function(formula1 = Y ~ Tr, formula2 = NULL, df, 
                            clust = NULL,
                            n = NULL,
                            block= NULL,
                            pred = NULL, ...){ 
  
  clust_name <- deparse(substitute(clust))
  block_name <- deparse(substitute(block))
  pred_name  <- deparse(substitute(pred))
  
  clust <- if(clust_name != "NULL") df[[clust_name]] else NULL
  block <- if(block_name != "NULL") df[[block_name]] else NULL
  pred  <- if(pred_name  != "NULL") df[[pred_name]]  else NULL
  
  if(!inherits(formula1, "formula")) stop("Formula1 must be a formula")
  
  if(!inherits(formula2, "formula") && !is.null(formula2)) stop("Formula2 must be a formula or NULL")
  
  if(!is.data.frame(df)) stop("Df must be a data frame")
  
  vars <- all.vars(formula1)
  Y_name <- vars[1]
  Tr_name <- vars[2]
  
  Y <- df[[Y_name]]
  Tr <- df[[Tr_name]]
  
  if(!is.vector(Y)) stop("Y must be a vector")
  
  if(!is.vector(Tr)) stop("Tr must be a vector")
  
  if(!all(Tr %in% c(0,1)) && !is.factor(Tr)) stop("Tr must be a binary vector of 0s and 1s or a categorical vector")
  
  if(!is.null(clust) && !is.vector(clust)) stop("Clust must be a vector")
  
  if(!is.null(n) && !is.vector(n)) stop("n must be a number that represents the cluster sizes or a vector of cluster sizes")
  
  if(!is.null(block) && !is.vector(block)) stop("Block must be a vector")
  
  
  if(!is.null(formula2)){
    
    Z_name <- all.vars(formula2)
    
    if(!all(Z_name %in% colnames(df))) stop("Variable in formula2 not found in df")
    
    Z <- model.matrix(formula2, data = df)
    
  } else {
    
    Z <- NULL
  }
  
  new_list <- list(Y = Y, Tr = Tr, Z = Z, df = df)
  
  if(!is.null(clust) && !is.null(block)){
    
    clust_df <- data.frame(clust = clust, block = block, Tr = Tr, Y = Y, Z = Z) |>
      group_by(clust, block, Tr) |>
      summarize(Y_agg = mean(Y),
                Z_agg = if(!is.null(Z)) mean(Z) else NULL,
                .groups = "drop")
    
    Y_agg  <- clust_df$Y_agg
    Tr_agg <- round(clust_df$Tr)
    P_agg  <- clust_df$block
    Z_agg <- clust_df$Z_agg
    
    message("Note: p_loop was used because data is in blocked clusters")
    return(p_loop(Y = Y_agg, Tr = Tr_agg, Z = Z_agg, P = P_agg, n = n, ...))
    
  } else if(!is.null(clust)&&is.null(block)){
    
    message("Note: p_loop was used because data is in clusters")
    return(p_loop(Y = Y, Tr = Tr, Z = Z, P = clust, n = n, ...))
    
  } else{
    message("Note: loop was used because no clusters or blocks were provided")
    return(loop(Y = Y, Tr = Tr, Z = Z, ...))
  }
  
  
}