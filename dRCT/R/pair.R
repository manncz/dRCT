#' Format Data into Pairs
#'
#' This function takes data at the unit level and returns a data frame with one pair per row. It also orients the data so
#' the treated unit is labeled as unit "1" in each row. Each row includes the pair label, the observed outcomes, and the covariates for each unit.
#' @param Y A vector of experimental outcomes.
#' @param Tr The treatment assignment vector.
#' @param Z A matrix or vector of covariates.
#' @param P A vector encoding the pair assignments. Observations with the same value of \code{P} will be assumed to be in the same pair.
#' @param n A vector encoding the cluster sizes.
#' @export

pair = function(Y,Tr,Z,P,n){
  
  data = data.frame(Y,Tr,Z,n,P) %>%
    group_by(P) %>%
    mutate(nP = n(),
           drop = is.na(P) | nP != 2)
  
  # give a warning if any units were removed due to missing the pair assignment
  if(sum(data$drop) > 0) warning(paste0("Dropped ",sum(data$drop)," units because no pair assignment was provided or pair assignment did not include exactly two units. Analysis is based on ", length(unique(data$P[!data$drop])), " pairs."))
  
  # remove units missing pair assignment or in a strata with no other or more than one other unit
  data = data %>%
    filter(!drop) %>%
    select(-nP, -drop)
  
  # first reshape as assigned
  agg = data %>%
    arrange(P) %>%
    mutate(label = row_number()) %>%
    pivot_wider(id_cols = P,
                names_from = label,
                values_from = c(Y,Tr,starts_with("Z"),n),
                names_glue = "{.value}{label}")
  
  n_assigned = agg %>%
    select(P, starts_with("n"))
  
  out = agg %>%
    select(-Tr2) %>%
    rename(Tr = Tr1)%>%
    ungroup()
  
  # reshape so that treated unit is first
  reorder = data %>%
    arrange(P, -Tr) %>% # this arranges things so that the treated unit is labeled as unit 1
    mutate(label = row_number()) %>%
    pivot_wider(id_cols = P,
                names_from = label,
                values_from = c(Y,Tr,starts_with("Z"),n),
                names_glue = "{.value}{label}") 
  
  ordered = reorder %>%
    select(!contains("Tr")) %>%
    ungroup()
  
  n_ordered = reorder %>%
    select(P, starts_with("n"))
  
  return(list(agg = out, n_assigned = n_assigned,
              ordered = ordered, n_ordered = n_ordered))
}
