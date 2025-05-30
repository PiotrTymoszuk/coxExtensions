% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/ml.R
\name{mlCox}
\alias{mlCox}
\alias{mlCox.default}
\alias{mlCox.coxex}
\title{Train a canonical Cox proportional hazard model and make predictions.}
\usage{
mlCox(data, newdata, ...)

\method{mlCox}{default}(data, newdata, formula, ...)

\method{mlCox}{coxex}(data, newdata, ...)
}
\arguments{
\item{data}{a data frame used to train the Cox model or
a \code{\link{coxex}} model.}

\item{newdata}{a data frame or a list of data frames used for model
validation.}

\item{...}{additional arguments passed to \code{\link[survival]{coxph}}.
Note that \code{x} and \code{y} are already used.}

\item{formula}{a formula of the Cox model; it must call a \code{Surv()} object,
i.e. contain a \code{'Surv()'} term.}
}
\value{
a list of \code{\link{coxex}} objects - an instance
of \code{\link{mlcx}} class.
}
\description{
A handy wrapper which fits a proportional hazard function via
\code{\link[survival]{coxph}} in a training data set (\code{data} argument) and
makes predictions in one or more validation data sets (\code{newdata}).
}
\details{
The workflow is as follows: (1) a canonical Cox proportional hazard model is
fitted in the training data set, (2) linear predictor scores are calculated
for the training and validation data sets, (3) a list of uni-variable Cox
models is generated, with each of them regressing the survival as a function
of the linear predictor scores.
The output list consists of \code{\link{coxex}} models wrapped in
a \code{\link{mlcx}} object.
Its assumptions, fit statistics, inference, and diagnostic plots can be
accessed with \code{summary()} and \code{plot()} methods called for single elements.
Alternatively, the assumptions, fit statistics, and inference can be computed
with the dedicated \code{\link{summary.mlcx}} method.
}
