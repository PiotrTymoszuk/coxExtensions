% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/functions.R
\name{get_cox_qc}
\alias{get_cox_qc}
\title{Get the extended residual table for a CoxPH model.}
\usage{
get_cox_qc(
  cox_model,
  type.predict = "lp",
  type.residuals = "martingale",
  data = NULL,
  ...
)
}
\arguments{
\item{cox_model}{a CoxpPH model or a coxex object.}

\item{type.predict}{type of the prediction, 'lp', linear predictor score by
default. See: \code{\link[survival]{predict.coxph}} for details.}

\item{type.residuals}{type of the residuals, 'martingale' by default.
See: \code{\link[survival]{residuals.coxph}} for details.}

\item{data}{the data frame used for the model construction. Ignored,
if coxex object provided.}

\item{...}{additional arguments passed to \code{\link[broom]{augment}}.}
}
\value{
a data frame with the predicted values, residuals and expected normal
values for the residuals.
}
\description{
Retrieves model residuals and predicted values.
A wrapper around \code{\link[broom]{augment}}.
}
\references{
\itemize{
\item Therneau, T. M. & Grambsch, P. M. Modeling Survival Data: Extending
the Cox Model. (Springer Verlag, 2000).
}
}
