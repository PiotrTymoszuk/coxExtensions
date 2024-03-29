% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/oop_coxex.R
\name{plot.coxex}
\alias{plot.coxex}
\title{Plot a coxex object.}
\usage{
\method{plot}{coxex}(
  x,
  type = c("fit", "residuals"),
  cust_theme = survminer::theme_survminer(),
  ...
)
}
\arguments{
\item{x}{a coxex object.}

\item{type}{type of the plots, 'fit' by default.}

\item{cust_theme}{customized plot theme provided by the user.}

\item{...}{extra arguments passed to \code{\link{get_cox_qc_plots}} or
\code{\link{plot_cox_fit}}.}
}
\value{
a ggplot or a list of ggplot residual plots.
}
\description{
Plots a coxex object: diagnostic graphs of residuals
(type = 'residuals') or outcome vs fitted Kaplan-Meier plot ('fit').
}
