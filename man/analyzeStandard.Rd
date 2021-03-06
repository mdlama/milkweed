% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/mw_ipm.R
\name{analyzeStandard}
\alias{analyzeStandard}
\title{Compute standard IPM population level metrics, such as the population
growth rate, left and right eigenvectors, and sensitivity and elasticit
matrices.}
\usage{
analyzeStandard(obj)
}
\arguments{
\item{obj}{A mwIPM model object.}
}
\value{
A mwIPM model object, with the results of the computation stored
in a list.
}
\description{
Compute standard IPM population level metrics, such as the population
growth rate, left and right eigenvectors, and sensitivity and elasticit
matrices.
}
\examples{
ipm \%<>\% analyzeStandard()
str(ipm$analysis$standard)
}
