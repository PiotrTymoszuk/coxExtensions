# Machine learning toolbox for Cox proportional hazard models

#' @include imports.R

  NULL

# Training of a Cox model and predictions --------

#' Train a canonical Cox proportional hazard model and make predictions.
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

    mlCox(data = train_data,
          newdata = newdata,
          formula = train_formula, ...)

  }

# Cross-validation of the Cox models -------

#' Cross-validation of canonical Cox proportional hazard models.
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

# Last one-out cross-validation -----------

#' Last-one-out cross-validation (LOOCV) of Cox proportional hazard models.
#'
#' @description
#' The function `loocvCox()` performs leave last-one-out cross-validation
#' (LOOCV) of a Cox proportional hazard model.
#'
#' @details
#' In LOOCV, a series of n Cox proportional hazard models is constructed, each
#' for n - 1 observations, where n stands for the total number of observations.
#' In other words, in each iteration, all observations except one are used to
#' construct a model.
#' For the renaming single observation, the value of the linear predictor score
#' of the model fit with the n - observation is generated.
#' The whole procedure is repeated n times, so that for each observation the
#' linear predictor scores are computed.
#' Finally, we use these leave-one-out linear predictor scores to construct a
#' univariable Cox proportional hazard model (an instance of \code{\link{coxex}}
#' class).
#' This cross-validated model can be subsequently evaluated with the standard
#' function of the package (e.g. with \code{\link{summary.coxex}} to test for
#' assumptions, and calculate performance statistics).
#'
#' @return a \code{\link{coxex}} object.
#'
#' @param data a data frame with survival information or a \code{\link{coxex}}
#' model.
#' @param formula a formula of the Cox model; it must call a `Surv()` object,
#' i.e. contain a `'Surv()'` term.
#' @param ... additional arguments passed to \code{\link[survival]{coxph}}.
#' Note that `x` and `y` are already used.
#'
#' @export

  loocvCox <- function(data, ...) UseMethod("loocvCox")

#' @rdname loocvCox
#' @export loocvCox.default
#' @export

  loocvCox.default <- function(data,
                               formula, ...) {

    ## input control -------

    if(!is.data.frame(data)) {

      stop("'data' has to be a data frame.", call. = FALSE)

    }

    data <- filter(data, complete.cases(data))

    form_terms <- as.character(formula)

    if(length(form_terms) < 3) {

      stop("A flawed formula provided.", call. = FALSE)

    }

    if(!stri_detect(form_terms[2], regex = "^Surv")) {

      stop(paste("'formula' must call a survival object, i.e.",
                 "contain a 'Surv()' term."),
           call. = FALSE)

    }

    if(!is_tibble(data) & !is.null(rownames(data))) {

      obs_names <- rownames(data)

    } else {

      obs_names <- paste0("rep_", 1:nrow(data))

    }

    ## training and test data --------

    row_ids <- set_names(1:nrow(data), obs_names)

    train_data <- map(row_ids, ~data[-.x, , drop = FALSE])

    test_data <- map(row_ids, ~data[.x, , drop = FALSE])

    ## construction of the models ----------

    train_models <- map(train_data,
                        function(dt) coxph(formula = formula,
                                           data = dt,
                                           x = TRUE,
                                           y = TRUE, ...))

    ## predictions -------

    oof_lp <-
      pmap(list(object = train_models,
                newdata = test_data),
           predict,
           type = "lp")

    data[["lp_score"]] <- map_dbl(oof_lp, ~.x[[1]])

    ## the coxex model for the OOF predictions ------

    oof_formula <- as.formula(paste(form_terms[2], "~ lp_score"))

    oof_call <- call2(.fn = "coxph",
                      formula = oof_formula,
                      data = data,
                      x = TRUE,
                      y = TRUE, ...)

    as_coxex(eval(oof_call), data)

  }

#' @rdname loocvCox
#' @export loocvCox.coxex
#' @export

  loocvCox.coxex <- function(data, ...) {

    stopifnot(is_coxex(data))

    train_data <- model.frame(data, type = "data")

    train_formula <- formula(data)

    loocvCox(data = train_data,
             formula = train_formula, ...)

  }

# END -------
