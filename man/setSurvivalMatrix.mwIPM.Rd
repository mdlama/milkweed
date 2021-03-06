% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/mw_ipm.R
\name{setSurvivalMatrix.mwIPM}
\alias{setSurvivalMatrix.mwIPM}
\title{Create survival matrix.}
\usage{
\method{setSurvivalMatrix}{mwIPM}(obj, update = TRUE, perturb = rep(0, 4))
}
\arguments{
\item{obj}{A mwIPM model object.}

\item{update}{Update dependencies?}

\item{perturb}{Parameter perturbation vector for sensitivity analysis.}
}
\value{
A mwIPM model object.
}
\description{
Create survival matrix.
}
