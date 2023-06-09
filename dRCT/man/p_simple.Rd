% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/p_simple.R
\name{p_simple}
\alias{p_simple}
\title{Simple Unbiased Estimator for Paired Clustered Randomized Trial}
\usage{
p_simple(
  ordered,
  assigned,
  n_assigned,
  loo = T,
  weighted_imp = F,
  true_sum = NULL
)
}
\arguments{
\item{ordered}{A matrix of pair experimental data that his been processed by the \code{pair} function, with the treatment pair first.}

\item{assigned}{A matrix of pair experimental data that his been processed by the \code{pair} function.}

\item{n_assigned}{A matrix of pair experimental data cluster sizes that has been processed by the \code{pair} function.}

\item{loo}{A Boolean indicating whether leave-one-pair-out imputation should be used.}

\item{weighted_imp}{A Boolean indicating whether weighted mean imputation should be used.}
}
\description{
This function is used to impute potential outcomes for the \code{p_loop} function. Pairs are treated
as the unit, and leave one out mean imputation is used, with no covariate adjustment.
}
