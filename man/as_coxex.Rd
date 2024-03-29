% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/oop_coxex.R
\name{as_coxex}
\alias{as_coxex}
\title{Create a coxex object model.}
\usage{
as_coxex(cox_model, data)
}
\arguments{
\item{cox_model}{a coxph model.}

\item{data}{a data frame used for model construction.}
}
\value{
a coxex object.
}
\description{
Generates a coxex model given a coxph object model and the
modeling data frame (used for model construction).
}
