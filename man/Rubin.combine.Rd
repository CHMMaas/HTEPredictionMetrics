% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/Rubin_combine.R
\name{Rubin.combine}
\alias{Rubin.combine}
\title{Rubin.combine}
\usage{
Rubin.combine(est, se)
}
\arguments{
\item{est}{vector of estimates}

\item{se}{vector of standard errors}
}
\value{
The output of the Rubin.combine function is a "list" with the following components.

est

Combined estimate


se

Combined standard error
}
\description{
This function combines estimates of multiple imputations.
}
\examples{
library(HTEPredictionMetrics)
set.seed(1)
m <- 5 # number of imputations
est <- runif(m)
se <- runif(m)
Rubin.combine(est, se)
}
