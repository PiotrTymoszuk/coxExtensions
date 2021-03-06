% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/oop.R
\name{residuals.coxex}
\alias{residuals.coxex}
\title{Extended residuals of a coxex object.}
\usage{
\method{residuals}{coxex}(object, ...)
}
\arguments{
\item{object}{a coxex object.}

\item{...}{additional arguments passed to \code{\link{get_cox_qc}}.}
}
\value{
a data frame with the predicted values, residuals and expected normal
values for the residuals.
}
\description{
Retrieves model residuals and predicted values.
A wrapper around \code{\link[broom]{augment}}.
}
