# S3 OOP for the coxex class.

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

    if(!any(class(cox_model) == 'coxph')) {

      stop('A valid coxph class model required.', call. = FALSE)

    }

    ## construction

    data_call <- rlang::enexpr(data)

    structure(list(model = cox_model,
                   data = rlang::quo(!!data_call)),
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

    coxExtensions::coxex(cox_model, data)

  }

# Class testing ----

#' Test for the coxex class.
#'
#' @description Tests if an object is an instance of the coxex class.
#' @param x an object.
#' @return a logical value.
#' @export

  is_coxex <- function(x) {

    all(class(x) == 'coxex')

  }

# Coertion -----
#'
#' Convert to a coxph class.
#'
#' @description Converts a coxex model to a coxph object.
#' @param x a coxex model.
#' @return a coxph model.
#' @export

  as_coxph <- function(x) {

    stopifnot(coxExtensions::is_coxex(x))

    x$model

  }

# Appearance -----

#' Print a coxex model.
#'
#' @description Prints a coxex model.
#' @param x a coxex object.
#' @param ... extra arguments, none specified.
#' @return none, called fot it's side effects.
#' @export

  print.coxex <- function(x, ...) {

    stopifnot(coxExtensions::is_coxex(x))

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

    stopifnot(coxExtensions::is_coxex(object))

    surv <- unclass(model.frame(object$model)[, 1])

    status <- surv[, 'status']

    tibble::tibble(observations = c('total', 'events'),
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

    stopifnot(coxExtensions::is_coxex(formula))

    type <- match.arg(type[1], c('model_frame', 'data', 'surv'))

    switch(type,
           model_frame = model.frame(formula$model),
           data = rlang::eval_tidy(formula$data),
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

    stopifnot(coxExtensions::is_coxex(x))

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
#' @export

  residuals.coxex <- function(object, ...) {

    stopifnot(coxExtensions::is_coxex(object))

    coxExtensions::get_cox_qc(cox_model = object, ...)

  }

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

    stopifnot(coxExtensions::is_coxex(object))

    survival:::predict.coxph(object = object$model, ...)

  }

# Summary ------

#' Summary for the coxex class.
#'
#' @description Generates a summary for a coxex model. This can be inference
#' (type = 'inference'), goodness of fit (type = 'fit') or
#' model assumptions (type = 'assumptions').
#' @details A wrapper around \code{\link{get_cox_estimates}} (inference),
#' \code{\link{get_cox_stats}} (goodness of fit) and
#' \code{\link{get_cox_assumptions}} (assumptions).
#' @return a data frame with the stats specified by the type argument.
#' @param object a coxex object.
#' @param type statistic type, 'inference' by default.
#' @param ... extra arguments passed to \code{\link{get_cox_estimates}},
#' \code{\link{get_cox_stats}} or \code{\link{get_cox_assumptions}}.
#' @export summary.coxex
#' @export

  summary.coxex <- function(object,
                            type = c('inference', 'fit', 'assumptions'), ...) {

    stopifnot(coxExtensions::is_coxex(object))

    type <- match.arg(type[1], c('inference', 'fit', 'assumptions'))

    switch(type,
           inference = coxExtensions::get_cox_estimates(object, ...),
           fit = coxExtensions::get_cox_stats(object, ...),
           assumptions = coxExtensions::get_cox_assumptions(object, ...))

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

    stopifnot(coxExtensions::is_coxex(x))

    type <- match.arg(type, c('fit', 'residuals'))

    ## plotting

    switch(type,
           residuals = coxExtensions::get_cox_qc_plots(cox_model = x,
                                                       cust_theme = cust_theme,
                                                       ...),
           fit = coxExtensions::plot_cox_fit(cox_model = x,
                                             cust_theme = cust_theme,
                                             ...))

  }

# END ----
