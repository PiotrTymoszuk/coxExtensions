% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/functions.R, R/oop_coxex.R
\name{get_cox_calibration}
\alias{get_cox_calibration}
\alias{calibrate.coxex}
\title{Calculate D'Agostino - Nam calibration stats.}
\usage{
get_cox_calibration(
  cox_model,
  n = 3,
  labels = NULL,
  right = FALSE,
  use_unique = FALSE,
  ...
)

\method{calibrate}{coxex}(fit, n = 3, labels = NULL, right = FALSE, use_unique = FALSE, ...)
}
\arguments{
\item{cox_model}{a CoxpPH model or a coxex object.}

\item{n}{a single numeric defining the number of quantile intervals of
the Cox model's linear predictor score.}

\item{labels}{an optional user-provided vector of labels for
the quantile intervals.}

\item{right}{logical, indicating if the quantile intervals should be closed
on the right (and open on the left) or vice versa.}

\item{use_unique}{logical, should unique values of the Cox model's linear
predictor scores be used instead of quantile strata? This option is
non-canonical and experimental, may be utilized to check an association
between an ordinal variable and survival.}

\item{...}{additional arguments, currently none.}

\item{fit}{a CoxPH model or a coxex object.}
}
\value{
an object of class 'calibrator', whose statistic and graphical
summary can be accessed by S3 'summary' and 'plot' methods.
The object consists of the following elements:
\itemize{
\item 'lp_scores': a data frame with the linear predictor scores and its
quantile intervals;
\item 'surv_fit' an object of the 'surv_fit' class
(\code{\link[survminer]{surv_fit}}), which stores the Kaplan - Meier
survival estimates;
\item 'cox_fit': a coxph model of the survival as a function of
the linear predictor score strata;
\item 'km_estimates': a data frame with the Kaplan-Meier survival estimates
for each observation, see: \code{\link[survminer]{surv_summary}} for details;
\item 'cox_estimates': a data frame storing the number of events and
survival probability for each observation in the linear predictor
strata Cox model;
\item 'strata_calibration': a data frame storing the numbers of events and
survivors, survival probabilities, relative survival probability
(Kaplan - Meier to Cox ratio), relative risk (rr: Kaplan-Meier to Cox ratio)
and the D'Agostino-Nam chi-squared statistic (x2_dn) for each strata;
\item 'global_calibration': a data frame with the mean relative risk,
the sum D'Agostino - Nam chi-squared statistic,
its degrees of freedom (df = strata number - 2 or df = 1 for two strata)
and p value for the global linear predictor score.
\item 'squares': a data frame with times, observed event, model-predicted
survival probability and squared distance between the prediction and observed
survival.
}
}
\description{
Calculates the D'Agostino - Nam calibration for a Coxph model,
a variant of the Hosmer - Lemeshov test (DOI: 10.1016/S0169-7161(03)23001-7),
and computes square distances between the predicted survival probability and
the outcome as proposed by Graf et al.
For  the D'Agostino - Nam calibration the linear predictor score of the Cox
model is cut into n quantile strata by \code{\link{cut_quantiles}}.
The 'observed' survival odds are derived from a Kaplan-Meier
survival estimate calculated by \code{\link[survminer]{surv_fit}}.
The 'fitted' survival obs are derived from the Cox proportional hazard models
of survival as a function of the linear predictor score quantile intervals.
}
\references{
\itemize{
\item D’Agostino, R. B. & Nam, B. H. Evaluation of the Performance of
Survival Analysis Models: Discrimination and Calibration Measures.
Handb. Stat. 23, 1–25 (2003).
\item Royston, P. Tools for checking calibration of a Cox model in external
validation: Approach based on individual event probabilities.
Stata J. 14, 738–755 (2014).
\item Crowson, C. S. et al. Assessing calibration of prognostic risk
scores. Stat. Methods Med. Res. 25, 1692–1706 (2016).
\item Graf, E., Schmoor, C., Sauerbrei, W. & Schumacher, M. Assessment and
comparison of prognostic classification schemes for survival data. Stat.
Med. 18, 2529–2545 (1999).
}
}
