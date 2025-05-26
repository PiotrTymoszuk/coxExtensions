# Machine learning toolbox for Cox models

#' @include imports.R

  NULL

# Training of a Cox model and predictions --------

#' Train a canonical Cox model and make predictions.
#'
#' @description
#' A handy wrapper which fits a proportional hazard function via
#' \code{\link[survival]{coxph}} in a training data set (`data` argument) and
#' makes predictions in one or more validation data sets (`newdata`).
#'
#' @details
#' The workflow is as follows: (1) a canonical Cox proportional hazard model is
#' fitted in the training data set, (2) linear predictor scores are calculated
#' for the training and validation data sets, (3) a list of uni-variable Cox
#' models is generated, with each of them regressing the survival as a function
#' of the linear predictor scores.
#' The output list consists of \code{\link{coxex}} models wrapped in
#' a \code{\link{mlcx}} object.
#' Its assumptions, fit statistics, inference, and diagnostic plots can be
#' accessed with `summary()` and `plot()` methods called for single elements.
#' Alternatively, the assumptions, fit statistics, and inference can be computed
#' with the dedicated \code{\link{summary.mlcx}} method.
#'
#' @return a list of \code{\link{coxex}} objects - an instance
#' of \code{\link{mlcx}} class.
#'
#' @param data a data frame used to train the Cox model or
#' a \code{\link{coxex}} model.
#' @param newdata a data frame or a list of data frames used for model
#' validation.
#' @param formula a formula of the Cox model; it must call a `Surv()` object,
#' i.e. contain a `'Surv()'` term.
#' @param ... additional arguments passed to \code{\link[survival]{coxph}}.
#' Note that `x` and `y` are already used.
#'
#' @export

  mlCox <- function(data, newdata, ...) UseMethod('mlCox')

#' @rdname mlCox
#' @export mlCox.default
#' @export

  mlCox.default <- function(data,
                            newdata,
                            formula,
                            ...) {

    ## input control -------

    if(!is.data.frame(data)) {

      stop("'data' has to be a data frame.", call. = FALSE)

    }

    err_txt <-
      "'newdata' has to be a data frame or a named list of data frames."

    if(!is.data.frame(newdata)) {

      if(!is.list(newdata)) stop(err_txt, call. = FALSE)

    }

    data <- filter(data, complete.cases(data))

    if(!is.data.frame(newdata)) {

      if(!is.list(newdata)) stop(err_txt, call. = FALSE)

      classes <- map_lgl(newdata, is.data.frame)

      if(any(!classes)) stop(err_txt, call. = FALSE)

      if(is.null(names(newdata))) stop(err_txt, call. = FALSE)

      newdata <- map(newdata, ~filter(.x, complete.cases(.x)))

      pred_data <- c(list(train = data),
                     test = newdata)

    } else {

      newdata <- filter(newdata, complete.cases(newdata))

      pred_data <- list(train = data,
                        test = newdata)

    }

    form_terms <- as.character(formula)

    if(length(form_terms) < 3) {

      stop('A flawed formula provided.', call. = FALSE)

    }

    if(!stri_detect(form_terms[2], regex = '^Surv')) {

      stop(paste("'formula' must call a survival object, i.e.",
                 "contain a 'Surv()' term."),
           call. = FALSE)

    }

    ## training and predictions -------

    train_model <- coxph(formula = formula,
                         data = data,
                         x = TRUE,
                         y = TRUE, ...)

    lp_scores <- map(pred_data,
                     ~predict(object = train_model,
                              newdata = .x,
                              type = 'lp'))

    observation <- NULL
    lp_score <- NULL

    lp_scores <- map(lp_scores,
                     ~tibble(observation = 1:length(.x),
                             lp_score = unname(.x)))

    lp_scores <- map2(lp_scores, pred_data, cbind)

    ## univariable Cox models -------

    lp_formula <- as.formula(paste(form_terms[2], '~ lp_score'))

    lp_models <- map(lp_scores,
                     ~call2(.fn = 'coxph',
                            formula = lp_formula,
                            x = TRUE,
                            y = TRUE, ...))

    lp_models <- map(lp_models, eval)

    mlcx(map2(lp_models, lp_scores, as_coxex),
          validation = 'external')

  }

#' @rdname mlCox
#' @export mlCox.coxex
#' @export

  mlCox.coxex <- function(data,
                          newdata,
                          ...) {

    stopifnot(is_coxex(data))

    train_data <- model.frame(data, type = 'data')

    train_formula <- formula(data)

    mlCox.default(data = train_data,
                  newdata = newdata,
                  formula = train_formula, ...)

  }

# Cross-validation of the Cox models -------

#' Cross-validation of canonical Cox models.
#'
#' @description
#' Cross-validation and repeated cross-validation of Cox proportional
#' hazard models.
#'
#' @details
#' The cross-validation procedure is essentially the same as described for
#' external validation with \code{\link{mlCox}}.
#' In brief, the training portion of each re-sample or CV fold is used to fit
#' a Cox model, which is used for generation of the linear predictor scores for
#' observations in the test portion of the re-sample. Finally, for each of
#' those test portions, a uni-variable Cox model is constructed which regresses
#' survival as a function of the linear predictor scores.
#' The output list consists of \code{\link{coxex}} models wrapped in
#' a \code{\link{mlcx}} object.
#' Its assumptions, fit statistics, inference, and diagnostic plots can be
#' accessed with `summary()` and `plot()` methods called for single elements.
#' Alternatively, the assumptions, fit statistics, and inference can be computed
#' with the dedicated \code{\link{summary.mlcx}} method.
#'
#' @return a list of \code{\link{coxex}} objects - an instance
#' of \code{\link{mlcx}} class.
#'
#' @param data a data frame used to train the Cox model or
#' a \code{\link{coxex}} model.
#' @param formula a formula of the Cox model; it must call a `Surv()` object,
#' i.e. contain a `'Surv()'` term.
#' @param n_folds number of cross-validation folds.
#' @param n_repeats number of repeats of the cross-validation procedure,
#' defaults to 1.
#' @param ... additional arguments passed to \code{\link[survival]{coxph}}.
#' Note that `x` and `y` are already used.
#'
#' @export

  cvCox <- function(data, ...) UseMethod('cvCox')

#' @rdname cvCox
#' @export cvCox.default
#' @export

  cvCox.default <- function(data,
                            formula,
                            n_folds = 10,
                            n_repeats = 1, ...) {

    ## input control -------

    if(!is.data.frame(data)) {

      stop("'data' has to be a data frame.", call. = FALSE)

    }

    data <- filter(data, complete.cases(data))

    stopifnot(is.numeric(n_folds))
    stopifnot(is.numeric(n_repeats))

    n_folds <- as.integer(n_folds)
    n_repeats <- as.integer(n_repeats)

    form_terms <- as.character(formula)

    if(length(form_terms) < 3) {

      stop('A flawed formula provided.', call. = FALSE)

    }

    if(!stri_detect(form_terms[2], regex = '^Surv')) {

      stop(paste("'formula' must call a survival object, i.e.",
                 "contain a 'Surv()' term."),
           call. = FALSE)

    }

    ## survival object and CV folds ----------

    eval_envir <- env(!!!as.list(data))

    surv_quo <- rlang::parse_quo(form_terms[2], eval_envir)

    surv_object <- unclass(eval_tidy(surv_quo))

    surv_status <- surv_object[, 'status']

    fold_ids <- createMultiFolds(surv_status, k = n_folds, times = n_repeats)

    ## out-of-fold models --------

    extra_args <- list2(...)

    oof_models <- map(fold_ids,
                      ~call2('mlCox.default',
                             data = data[.x, ],
                             newdata = data[-.x, ],
                             formula = formula,
                             !!!extra_args))

    oof_models <- map(oof_models,
                      ~eval(.x)[[2]])

    mlcx(oof_models, validation = 'cv')

  }

#' @rdname cvCox
#' @export cvCox.coxex
#' @export

  cvCox.coxex <- function(data,
                          n_folds = 10,
                          n_repeats = 1, ...) {

    stopifnot(is_coxex(data))

    train_data <- model.frame(data, type = 'data')
    train_formula <- formula(data)

    cvCox.default(data = train_data,
                  formula = train_formula,
                  n_folds = n_folds,
                  n_repeats = n_repeats, ...)

  }

# END -------
