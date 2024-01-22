# S3 OOP for the coxex class.

#' @include imports.R

  NULL

# Constructor ------

#' Create a coxex object model.
#'
#' @description Generates a coxex model given a coxph object model and the
#' modeling data frame (used for model construction).
#' @param cox_model a coxph model.
#' @param data a data frame used for model construction.
#' @return a coxex object.
#' @export

  coxex <- function(cox_model, data) {

    ## entry control

    if(!is.data.frame(data)) {

      stop('A data frame required as the data argument.', call. = FALSE)

    }

    if(!inherits(cox_model, 'coxph')) {

      stop('A valid coxph class model required.', call. = FALSE)

    }

    ## construction

    data_call <- enexpr(data)

    structure(list(model = cox_model,
                   data = quo(!!data_call)),
              class = 'coxex')

  }

#' Create a coxex object model.
#'
#' @description Generates a coxex model given a coxph object model and the
#' modeling data frame (used for model construction).
#' @inheritParams coxex
#' @return a coxex object.
#' @export

  as_coxex <- function(cox_model, data) {

    coxex(cox_model, data)

  }

# Class testing ----

#' Test for the coxex class.
#'
#' @description Tests if an object is an instance of the coxex class.
#' @param x an object.
#' @return a logical value.
#' @export

  is_coxex <- function(x) inherits(x, 'coxex')

# Coercion -----
#'
#' Convert to a coxph class.
#'
#' @description Converts a coxex model to a coxph object.
#' @param x a coxex model.
#' @return a coxph model.
#' @export

  as_coxph <- function(x) {

    stopifnot(is_coxex(x))

    x$model

  }

# Appearance -----

#' Print a coxex model.
#'
#' @description Prints a coxex model.
#' @param x a coxex object.
#' @param ... extra arguments, none specified.
#' @return none, called for its side effects.
#' @export

  print.coxex <- function(x, ...) {

    stopifnot(is_coxex(x))

    print(x$model)

  }

# Extraction ------

#' Number of observations and events.
#'
#' @description Extracts the number of observations and events used for
#' modeling.
#' @param object a coxex object.
#' @param ... extra arguments, currently none.
#' @return a data frame with the number of total complete observations and
#' events.
#' @export

  nobs.coxex <- function(object, ...) {

    stopifnot(is_coxex(object))

    surv <- unclass(model.frame(object$model)[, 1])

    status <- surv[, 'status']

    tibble(observations = c('total', 'events'),
           n = c(nrow(surv), sum(status)))

  }

#' Modeling data.
#'
#' @description Access to the model components and data for the coxex class.
#' @param formula a coxex object.
#' @param type the output type: 'model_frame' returns a model frame,
#' 'data' return the data fame used for model construction,
#' 'surv' accesses the survival object.
#' @param ... extra arguments, currently none.
#' @return a data frame or a survival object.
#' @export model.frame.coxex
#' @export

  model.frame.coxex <- function(formula,
                                type = c('model_frame', 'data', 'surv'), ...) {

    stopifnot(is_coxex(formula))

    type <- match.arg(type[1], c('model_frame', 'data', 'surv'))

    switch(type,
           model_frame = model.frame(formula$model),
           data = eval_tidy(formula$data),
           surv = model.frame(formula$model)[, 1])

  }

#' Modeling formula.
#'
#' @description Accesses the modeling formula of a coxex object.
#' @param x a coxex object.
#' @param ... extra arguments, currently none.
#' @return a formula
#' @export

  formula.coxex <- function(x, ...) {

    stopifnot(is_coxex(x))

    formula(x$model)

  }

# Residuals and prediction -------

#' Extended residuals of a coxex object.
#'
#' @description Retrieves model residuals and predicted values.
#' A wrapper around \code{\link[broom]{augment}}.
#' @return a data frame with the predicted values, residuals and expected normal
#' values for the residuals.
#' @param object a coxex object.
#' @param ... additional arguments passed to \code{\link{get_cox_qc}}.
#' @export residuals.coxex
#' @references
#' * Therneau, T. M. & Grambsch, P. M. Modeling Survival Data: Extending
#' the Cox Model. (Springer Verlag, 2000).
#' @md
#' @export

  residuals.coxex <- function(object, ...) {

    stopifnot(is_coxex(object))

    get_cox_qc(cox_model = object, ...)

  }

#' @rdname residuals.coxex
#' @export resid.coxex
#' @export

  resid.coxex <- residuals.coxex

#' Prediction for the coxex class.
#'
#' @description Predicts survival, linear predictors etc. for a coxex model.
#' Technically, a wrapper around \code{\link[survival]{predict.coxph}}.
#' @param object a coxex object.
#' @param ... additional arguments passed
#' to \code{\link[survival]{predict.coxph}}.
#' @return a vector or matrix of predictions, or a list containing the
#' predictions.
#' @export predict.coxex
#' @export

  predict.coxex <- function(object, ...) {

    stopifnot(is_coxex(object))

    predict(object = object$model, ...)

  }

# Summary ------

#' Summary for the coxex class.
#'
#' @description Generates a summary for a coxex model. This can be inference
#' (type = 'inference'), fit statistics in the training data (type = 'fit') or
#' model assumptions (type = 'assumptions').
#' @details A wrapper around \code{\link{get_cox_estimates}} (inference),
#' \code{\link{get_cox_stats}} (goodness of fit) and
#' \code{\link{get_cox_assumptions}} (assumptions).
#' @return a data frame with the stats specified by the type argument.
#' @param object a coxex object.
#' @param type statistic type, 'inference' by default.
#' @param ... extra arguments passed to \code{\link{get_cox_estimates}},
#' \code{\link{get_cox_stats}} or \code{\link{get_cox_assumptions}}.
#'
#' @references
#' * Therneau, T. M. & Grambsch, P. M. Modeling Survival Data: Extending
#' the Cox Model. (Springer Verlag, 2000).
#' * Grambsch, P. M. & Therneau, T. M. Proportional Hazards Tests and
#' Diagnostics Based on Weighted Residuals. Biometrika 81, 515 (1994).
#' * Harrell, F. E., Lee, K. L. & Mark, D. B. Multivariable prognostic
#' models: Issues in developing models, evaluating assumptions and adequacy,
#' and measuring and reducing errors. Stat. Med. 15, 361–387 (1996).
#' * Graf, E., Schmoor, C., Sauerbrei, W. & Schumacher, M. Assessment and
#' comparison of prognostic classification schemes for survival data.
#' Stat. Med. 18, 2529–2545 (1999).
#'
#' @md
#' @export summary.coxex
#' @export

  summary.coxex <- function(object,
                            type = c('inference', 'fit', 'assumptions'), ...) {

    stopifnot(is_coxex(object))

    type <- match.arg(type[1], c('inference', 'fit', 'assumptions'))

    switch(type,
           inference = get_cox_estimates(object, ...),
           fit = get_cox_stats(object, ...),
           assumptions = get_cox_assumptions(object, ...))

  }

# Brier scores ------

#' Calculate Brier scores for the survival fit.
#'
#' @description The Brier score for survival models is computed as a squared
#' distance between the predicted survival probability and and the actual
#' 0/1-coded survival for each unique time point of the dataset.
#' @details See, \code{\link[pec]{pec}} for details.
#' The `surv_brier` is a S3 generic function.
#' @return an object of the \code{\link{brier}} class with the `plot` method
#' and multiple methods shared with traditional data frames.
#' @param fit a survival model.
#' @param data the data frame used for the model construction.
#' @param splitMethod validation method used for calculation of the Bier scores
#' as specified by the respective \code{\link[pec]{pec}} argument.
#' One of 'none' (no validation), 'cvK' (K-fold cross-validation, e.g. 'cv10'),
#' 'boot' (bootstrap), 'BootCv' (bootstrap cross-validation),
#' 'Boot632', 'Boot632plus', 'loocv' or 'NoInf', see the upstream
#' function for details.
#' @param ... extra arguments passed to \code{\link[pec]{pec}}.
#'
#' @references
#' * Graf, E., Schmoor, C., Sauerbrei, W. & Schumacher, M. Assessment and
#' comparison of prognostic classification schemes for survival data.
#' Stat. Med. 18, 2529–2545 (1999).
#'
#' @md
#' @export surv_brier.coxex
#' @export

  surv_brier.coxex <- function(fit,
                               splitMethod = 'none', ...) {

    stopifnot(is_coxex(fit))

    get_cox_pec(cox_model = fit,
                type = 'brier',
                splitMethod = splitMethod, ...)

  }

#' @rdname surv_brier.coxex
#' @export surv_brier.coxph
#' @export

  surv_brier.coxph <- function(fit,
                               data,
                               splitMethod = 'none', ...) {

    get_cox_pec(cox_model = fit,
                data = data,
                type = 'brier',
                splitMethod = splitMethod, ...)

  }

# Calibration -------

#' @rdname get_cox_calibration
#' @export calibrate.coxex
#' @export

  calibrate.coxex <- function(fit,
                              n = 3,
                              labels = NULL,
                              right = FALSE,
                              use_unique = FALSE, ...) {

    stopifnot(is_coxex(fit))

    get_cox_calibration(cox_model = fit,
                        n = n,
                        labels = labels,
                        right = right,
                        use_unique = use_unique)

  }

# Validation --------

#' Validate a 'coxex' model.
#'
#' @description Provides an access to validation stats obtained e.g.
#' by cross-validation or bootstraping via \code{\link[rms]{validate}}.
#' @details See: \code{\link[rms]{validate.cph}}.
#' @return a data frame with the following variables:
#' * `dataset`: dataset used for computation of the stats
#' * `Dxy`: Somers' DXY rank correlation
#' * `R2`: Nagelkerke R-squared
#' * `Slope`: slope shrinkage
#' * `D`: discrimination index D
#' * `U`: unreliability index
#' * `Q`: the overall quality index
#' * `g`: g-index on the log relative hazard
#' * `c_index`: Harrell's concordance index.
#'
#' @param fit a 'coxex' model.
#' @param method resampling method may be "crossvalidation",
#' "boot" (the default), ".632", or "randomization".
#' @param B number of repetitions.
#' @param bw TRUE to do fast step-down using the \code{\link[rms]{fastbw}}
#' function, for both the overall model and for each repetition.
#' \code{\link[rms]{fastbw}} keeps parameters together that represent
#' the same factor.
#' @param rule Applies if bw = TRUE. "aic" to use Akaike's information
#' criterion as a stopping rule (i.e., a factor is deleted if the chi-squared
#' falls below twice its degrees of freedom), or "p" to use p-values.
#' @param type "residual" or "individual"
#' stopping rule is for individual factors or for the residual chi-squared
#' for all variables deleted
#' @param sls significance level for a factor to be kept in a model, or for
#' judging the residual chi-squared
#' @param aics 	cutoff on AIC when rule="aic"
#' @param force see \code{\link[rms]{fastbw}}
#' @param estimates see \code{\link[rms]{print.fastbw}}
#' @param pr TRUE to print results of each repetition
#' @param ... extra arguments passed to \code{\link[rms]{validate.cph}}
#'
#' @references
#' * Harrell, F. E., Lee, K. L. & Mark, D. B. Multivariable prognostic
#' models: Issues in developing models, evaluating assumptions and adequacy,
#' and measuring and reducing errors. Stat. Med. 15, 361–387 (1996).
#'
#' @md
#' @export validate.coxex
#' @export

  validate.coxex <- function(fit,
                             method = 'boot',
                             B = 40,
                             bw = FALSE,
                             rule = 'aic',
                             type = 'residual',
                             sls = 0.05,
                             aics = 0,
                             force = NULL,
                             estimates = TRUE,
                             pr = FALSE, ...) {

    stopifnot(is_coxex(fit))

    get_cox_validation(cox_model = fit,
                       method = method,
                       B = B,
                       bw = bw,
                       rule = rule,
                       type = type,
                       sls = sls,
                       aics = aics,
                       force = force,
                       estimates = estimates,
                       pr = pr, ...)

  }

# Plotting ------

#' Plot a coxex object.
#'
#' @description Plots a coxex object: diagnostic graphs of residuals
#' (type = 'residuals') or outcome vs fitted Kaplan-Meier plot ('fit').
#' @param x a coxex object.
#' @param type type of the plots, 'fit' by default.
#' @param cust_theme customized plot theme provided by the user.
#' @param ... extra arguments passed to \code{\link{get_cox_qc_plots}} or
#' \code{\link{plot_cox_fit}}.
#' @return a ggplot or a list of ggplot residual plots.
#' @export plot.coxex
#' @export

  plot.coxex <- function(x,
                         type = c('fit', 'residuals'),
                         cust_theme = survminer::theme_survminer(), ...) {

    ## entry control

    stopifnot(is_coxex(x))

    type <- match.arg(type, c('fit', 'residuals'))

    ## plotting

    switch(type,
           residuals = get_cox_qc_plots(cox_model = x,
                                        cust_theme = cust_theme,
                                        ...),
           fit = plot_cox_fit(cox_model = x,
                              cust_theme = cust_theme,
                              ...))

  }

# END ----
