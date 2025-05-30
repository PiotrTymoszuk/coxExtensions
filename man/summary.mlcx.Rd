% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/oop_mlcox.R
\name{summary.mlcx}
\alias{summary.mlcx}
\title{Inference, fit statistics, and assumptions for \code{mlcx} objects.}
\usage{
\method{summary}{mlcx}(
  object,
  type = c("inference", "fit", "assumptions"),
  collapse = TRUE,
  ci_type = c("percentile", "bca"),
  ...
)
}
\arguments{
\item{object}{a coxex object.}

\item{type}{statistic type, 'inference' by default.}

\item{collapse}{a logical, should the results be collapsed by
cross-validation folds (mean, SD, and 95\% confidence intervals)?
Refers only to objects generated with \code{\link{cvCox}} and
\code{type = 'inference'} or \code{type = 'fit'}, and ignored otherwise.}

\item{ci_type}{type of confidence intervals: percentile (default) or Efron's
BCA computed by \code{\link[coxed]{bca}}.}

\item{...}{extra arguments passed to \code{\link{get_cox_estimates}},
\code{\link{get_cox_stats}} or \code{\link{get_cox_assumptions}}.}
}
\value{
a data frame.
}
\description{
Generates a summary for a \code{mlcx} object.
This can be inference (\code{type = 'inference'}), fit statistics (\code{type = 'fit'})
or model assumptions (\code{type = 'assumptions'}).
The results are wrapped in a data frame.
}
\details{
A wrapper around \code{\link{summary.coxex}}.
}
