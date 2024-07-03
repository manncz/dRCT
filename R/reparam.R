#' Reparameterize Paired Data
#'
#' This function takes a data frame that has been processed by the \code{pair} function and formats it such that 
#' covariates are parameterized as the within pair differences (treatment - control) and means for each covariate.
#' @param assigned A matrix of pair experimental data that his been processed by the \code{pair} function.
#' @export

# re-parameterize the covariates
reparam = function(assigned){
  
  covs_long = assigned %>%
    group_by(P) %>%
    pivot_longer(starts_with("Z"),
                 names_to = c("var", "label"),
                 values_to = "val",
                 names_pattern = "(.*)(\\d)")%>%
    mutate(dif = case_when(label == 1 ~ val,
                           TRUE ~ -val)) %>%
    group_by(P, var) %>%
    summarize(avg = mean(val), dif = sum(dif))
  
  covs = covs_long %>%
    pivot_wider(id_cols = P,
                names_from = var,
                values_from = c(avg, dif),
                names_glue = "{var}_{.value}")
  
  out = assigned %>%
    select(-starts_with("Z")) %>%
    left_join(covs, by = "P")
  
  return(out)
}

