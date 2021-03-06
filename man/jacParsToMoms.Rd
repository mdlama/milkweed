% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/mw_ipm.R
\name{jacParsToMoms}
\alias{jacParsToMoms}
\title{Jacobian of parameters to moments.}
\usage{
jacParsToMoms(x, type = "ident")
}
\arguments{
\item{x}{Vector of parameters.}

\item{type}{Distribution type (lnorm, gamma, or identity)}
}
\value{
Jacobian matrix
}
\description{
Jacobian of parameters to moments.
}
