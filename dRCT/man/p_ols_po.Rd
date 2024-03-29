% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/p_ols_po.R
\name{p_ols_po}
\alias{p_ols_po}
\title{OLS Imputation of Potential Outcomes in Paired Experiments}
\usage{
p_ols_po(ordered, assigned, n_assigned, weighted_imp = F)
}
\arguments{
\item{ordered}{A matrix of pair experimental data that his been processed by the \code{pair} function, with the treatment pair first.}

\item{assigned}{A matrix of pair experimental data that his been processed by the \code{pair} function.}

\item{n_assigned}{A matrix of pair experimental data cluster sizes that has been processed by the \code{pair} function.}

\item{weighted_imp}{A Boolean indicating whether to weight imputation models}
}
\description{
This function is used to impute potential outcomes for the \code{p_loop} function.
Leaves out each pair and imputes its potential outcomes using ordinary least squares
regression on the remaining pairs. Individuals are treated as units when making predictions.
}
