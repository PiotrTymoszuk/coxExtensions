% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/ml.R
\name{loocvCox}
\alias{loocvCox}
\alias{loocvCox.default}
\alias{loocvCox.coxex}
\title{Last-one-out cross-validation (LOOCV) of Cox proportional hazard models.}
\usage{
loocvCox(data, ...)

\method{loocvCox}{default}(data, formula, ...)

\method{loocvCox}{coxex}(data, ...)
}
\arguments{
\item{data}{a data frame with survival information or a \code{\link{coxex}}
model.}

\item{...}{additional arguments passed to \code{\link[survival]{coxph}}.
Note that \code{x} and \code{y} are already used.}

\item{formula}{a formula of the Cox model; it must call a \code{Surv()} object,
i.e. contain a \code{'Surv()'} term.}
}
\value{
a \code{\link{coxex}} object.
}
\description{
The function \code{loocvCox()} performs leave last-one-out cross-validation
(LOOCV) of a Cox proportional hazard model.
}
\details{
In LOOCV, a series of n Cox proportional hazard models is constructed, each
for n - 1 observations, where n stands for the total number of observations.
In other words, in each iteration, all observations except one are used to
construct a model.
For the renaming single observation, the value of the linear predictor score
of the model fit with the n - observation is generated.
The whole procedure is repeated n times, so that for each observation the
linear predictor scores are computed.
Finally, we use these leave-one-out linear predictor scores to construct a
univariable Cox proportional hazard model (an instance of \code{\link{coxex}}
class).
This cross-validated model can be subsequently evaluated with the standard
function of the package (e.g. with \code{\link{summary.coxex}} to test for
assumptions, and calculate performance statistics).
}
