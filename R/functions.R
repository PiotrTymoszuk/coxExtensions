# Functional tools for the CoxPH models.

#' @include imports.R

  NULL

# Inference ----

#' Get inference statistic for a CoxPH model.
#'
#' @description Retrieves inference statistic for a Cox proportional hazard
#' model. A wrapper around \code{\link[survival]{summary.coxph}}.
#' @return a data frame with the model estimates, confidence intervals
#' and p values.
#' @param cox_model a CoxpPH model or a coxex object.
#' @param trans_function function used for transformation of the estimates and
#' confidence intervals, identity() by default.
#' @param ... extra arguments, currently none.
#' @references
#' * Therneau, T. M. & Grambsch, P. M. Modeling Survival Data: Extending
#' the Cox Model. (Springer Verlag, 2000).
#' @md
#' @export

  get_cox_estimates <- function(cox_model, trans_function = identity, ...) {

    ## suppressing R-CMD check notes on NSE variables

    estimate <- parameter <- variable <- level <- NULL

    ## entry control

    stopifnot(is.function(trans_function))

    if(is_coxex(cox_model)) {

      cox_model <- cox_model$model

    }

    if(!inherits(cox_model, 'coxph')) {

      stop('Please provide a valid coxph class model.', call. = FALSE)

    }

    ## meta information

    mod_frame <- model.frame(cox_model)

    modeling_vars <- names(mod_frame)[-1]

    extr_regex <- paste(modeling_vars, collapse = '|')

    ## confidence intervals

    ci_table <- as.data.frame(confint(cox_model))

    ci_table <- set_names(ci_table, c('lower_ci', 'upper_ci'))

    ci_table <- map_dfc(ci_table, trans_function)

    ## coefficients

    coefs <- as.data.frame(summary(cox_model)$coefficients)

    coefs <- set_names(coefs[, c('coef', 'se(coef)', 'z', 'Pr(>|z|)')],
                       c('estimate', 'se', 'stat', 'p_value'))

    coefs <- rownames_to_column(coefs, 'parameter')

    coefs <- mutate(coefs,
                    estimate = trans_function(estimate),
                    se = trans_function(estimate),
                    n_complete = nrow(model.frame(cox_model)),
                    stat_name = 'z')

    ## output table

    summ_tbl <- as_tibble(cbind(coefs, ci_table))

    summ_tbl <-
      mutate(summ_tbl,
             variable = stri_extract(parameter,
                                     regex = extr_regex),
             level = stri_replace(parameter,
                                  regex = extr_regex,
                                  replacement = ''))

    ## counting the n for the

    mod_counts <- map_dfr(modeling_vars,
                          ~count_(data = mod_frame, variable = .x))

    mod_counts <- mutate(mod_counts,
                         parameter = paste0(variable, level))

    summ_tbl <- left_join(summ_tbl,
                          mod_counts[c('parameter', 'n')],
                          by = 'parameter')

    ## output

    summ_tbl[c('parameter',
               'variable',
               'level',
               'n',
               'n_complete',
               'stat_name',
               'stat',
               'estimate',
               'se',
               'lower_ci',
               'upper_ci',
               'p_value')]


  }

# Extended residual table ------

#' Get the extended residual table for a CoxPH model.
#'
#' @description Retrieves model residuals and predicted values.
#' A wrapper around \code{\link[broom]{augment}}.
#' @return a data frame with the predicted values, residuals and expected normal
#' values for the residuals.
#' @param cox_model a CoxpPH model or a coxex object.
#' @param type.predict type of the prediction, 'lp', linear predictor score by
#' default. See: \code{\link[survival]{predict.coxph}} for details.
#' @param type.residuals type of the residuals, 'martingale' by default.
#' See: \code{\link[survival]{residuals.coxph}} for details.
#' @param data the data frame used for the model construction. Ignored,
#' if coxex object provided.
#' @param ... additional arguments passed to \code{\link[broom]{augment}}.
#' @references
#' * Therneau, T. M. & Grambsch, P. M. Modeling Survival Data: Extending
#' the Cox Model. (Springer Verlag, 2000).
#' @md
#' @export

  get_cox_qc <- function(cox_model,
                         type.predict = 'lp',
                         type.residuals = 'martingale',
                         data = NULL, ...) {

    ## suppressing a R-CMD-check note on NSE variables

    .resid <- .std.resid <- NULL

    ## entry control

    if(is_coxex(cox_model)) {

      data <- model.frame(cox_model, type = 'data')

      cox_model <- cox_model$model

    }

    if(!inherits(cox_model, 'coxph')) {

      stop('Please provide a valid coxph class model.', call. = FALSE)

    }

    ## computation

    resid_tbl <- augment(x = cox_model,
                         data = data,
                         type.predict = type.predict,
                         type.residuals = type.residuals, ...)

    resid_tbl <- mutate(resid_tbl,
                        .std.resid = scale(.resid)[, 1],
                        .observation = 1:nrow(resid_tbl),
                        .sq.std.resid = .std.resid^2,
                        .candidate_missfit = ifelse(abs(.std.resid) > qnorm(0.975), 'yes', 'no'))

    calc_expected_(resid_tbl, '.std.resid')

  }

# Goodness of fit ------

#' Goodness of fit for a CoxPH model.
#'
#' @description Calculates fit statistics for a CoxPH model: AIC (Akaike
#' information criterion), BIC (Bayesian information criterion), raw_rsq
#' (unadjusted R-squared, cod, mer or mev/default, see:
#' \code{\link[survMisc]{rsq}}), MAE (mean absolute error), MSE
#' (mean squared error), RMSE (root MSE),  concordance (C) index
#' with 95% confidence intervals (normality assumption) and the
#' integrated Brier score (IBS, see: \code{\link[pec]{pec}}).
#' @return a data frame with the statistic values.
#' @param rsq_type type of R-squared statistic, see: \code{\link[survMisc]{rsq}}.
#' @param ... extra arguments passed to \code{\link{get_cox_qc}}.
#' @inheritParams get_cox_qc
#' @references
#' * Therneau, T. M. & Grambsch, P. M. Modeling Survival Data: Extending
#' the Cox Model. (Springer Verlag, 2000).
#' * Harrell, F. E., Lee, K. L. & Mark, D. B. Multivariable prognostic
#' models: Issues in developing models, evaluating assumptions and adequacy,
#' and measuring and reducing errors. Stat. Med. 15, 361–387 (1996).
#' @md
#' @export

  get_cox_stats <- function(cox_model,
                            data = NULL,
                            type.predict = 'lp',
                            type.residuals = 'martingale',
                            rsq_type = c('mev', 'mer', 'cod'), ...) {

    ## entry control ----

    if(is_coxex(cox_model)) {

      data <- model.frame(cox_model, type = 'data')

      cox_model <- cox_model$model

    }

    if(!inherits(cox_model, 'coxph')) {

      stop('Please provide a valid coxph class model.', call. = FALSE)

    }

    rsq_type <- match.arg(rsq_type[1], c('mev', 'mer', 'cod'))

    ## R-squared ----

    rsq_tbl <- rsq(cox_model)

    ## concordance ----

    c_index <- unname(summary(cox_model)$concordance)

    se <- c_index[2]

    c_tbl <- tibble(c_index = c_index[1],
                    lower_ci = c_index[1] + qnorm(0.025) * se,
                    upper_ci = c_index[1] + qnorm(0.975) * se)

    ## integrated Brier score ------

    ibs_tbl <- suppressMessages(get_cox_pec(cox_model,
                                            data = data,
                                            return_pec = FALSE))

    ## n numbers -----

    surv <- unclass(model.frame(cox_model)[, 1])

    n_num <- c('total' = nrow(surv),
               'events' = sum(surv[, 'status']))

    ## errors ------

    resid_tbl <- get_cox_qc(cox_model = cox_model,
                            type.predict = type.predict,
                            type.residuals = type.residuals,
                            data = data, ...)

    ## output -------

    res_tbl <- tibble(n_complete = n_num['total'],
                      n_events = n_num['events'],
                      aic = AIC(cox_model),
                      bic = BIC(cox_model),
                      raw_rsq = rsq_tbl[[rsq_type]],
                      mae = mean(abs(resid_tbl$.resid)),
                      mse = mean(resid_tbl$.resid^2),
                      rmse = sqrt(mean(resid_tbl$.resid^2)))

    as_tibble(cbind(res_tbl, c_tbl, ibs_tbl))

  }

# Integrated Brier scores and prediction error curves -------

#' Calculate prediction error curves and integrated Brier score (IBS).
#'
#' @description Computes prediction error curves and
#' integrated Brier score (IBS)
#' for a Cox model using \code{\link[pec]{pec}} and \code{\link[pec]{crps}}.
#' @return a numeric vector with IBS values, a \code{\link[pec]{pec}} class
#' or a \code{\link{brier}} class object.
#' @param cox_model a Cox model of the `coxph` or `coxex` class.
#' @param data the data frame used for the model construction. Ignored,
#' if coxex object provided.
#' @param type type of object to be returned: 'ibs' returns
#' the Integrated Brier Score (IBS), 'pec' returns a \code{\link[pec]{pec}}
#' object and 'brier' returns a data frame with Brier scores for unique time
#' points of the \code{\link{brier}} class.
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
#' @export

  get_cox_pec <- function(cox_model,
                          data = NULL,
                          type = c('ibs', 'pec', 'brier'),
                          splitMethod = 'none', ...) {

    ## entry control ----

    if(is_coxex(cox_model)) {

      data <- model.frame(cox_model, type = 'data')

      cox_model <- cox_model$model

    }

    if(!inherits(cox_model, 'coxph')) {

      stop('Please provide a valid coxph class model.', call. = FALSE)

    }

    type <- match.arg(type[1], c('ibs', 'pec', 'brier'))

    ## generation of the pec object ------

    ## KM formula

    km_formula <- paste(as.character(formula(cox_model))[2], '~ 1')

    km_formula <- as.formula(km_formula)

    pec_obj <- pec(object = cox_model,
                   formula = km_formula,
                   data = data,
                   splitMethod = splitMethod, ...)

    if(type == 'pec') return(pec_obj)

    if(type == 'ibs') {

      times <- pec_obj$minmaxtime

      out <- crps(object = pec_obj,
                  times = pec_obj$minmaxtime,
                  start = pec_obj$start)

      return(tibble(ibs_reference = out[1, 1],
                    ibs_model = out[2, 1]))

    }

    times <- pec_obj$time
    reference <- pec_obj$AppErr[[1]]
    training <- pec_obj$AppErr[[2]]

    if(stri_detect(splitMethod, regex = '^cv')) {

      test <- pec_obj$crossvalErr[[2]]

    } else {

      test <-
        switch(splitMethod,
               none = NULL,
               noPlan = NULL,
               boot = pec_obj$BootCvErr[[2]],
               BootCv = pec_obj$BootCvErr[[2]],
               Boot632 = pec_obj$Boot632Err[[2]],
               Boot632plus = pec_obj$Boot632plusErr[[2]],
               loocv = pec_obj$loocvErr[[2]],
               NoInf = NULL)

    }

    return(brier(times = times,
                 reference = reference,
                 training = training,
                 test = test))

  }

# Model assumptions -----

#' Check assumptions fo a CoxPH model.
#'
#' @description Checks the normality and proportional hazard assumption for
#' a CoxPH model. Technically, the normality assumption is tested with
#' Shapiro-wil test, the proportionality by \code{\link[survival]{cox.zph}}.
#' @return a data frame with the normality testing results and the
#' proportional hazard assumption testing for the model variables and the
#' global model.
#' @inheritParams get_cox_qc
#' @param ... extra arguments passed to \code{\link{get_cox_qc}}
#' @references
#' * Grambsch, P. M. & Therneau, T. M. Proportional Hazards Tests and
#' Diagnostics Based on Weighted Residuals. Biometrika 81, 515 (1994).
#' @md
#' @export

  get_cox_assumptions <- function(cox_model,
                                  type.predict = 'lp',
                                  type.residuals = 'martingale',
                                  data = NULL, ...) {

    ## suppression of R-CMD-check notes concerning global variables

    chisq <- p <- NULL

    ## entry control

    if(is_coxex(cox_model)) {

      data <- model.frame(cox_model, type = 'data')

      cox_model <- cox_model$model

    }

    if(!inherits(cox_model, 'coxph')) {

      stop('Please provide a valid coxph class model.', call. = FALSE)

    }

    ## normality

    resid_tbl <- get_cox_qc(cox_model = cox_model,
                            type.predict = type.predict,
                            type.residuals = type.residuals,
                            data = data, ...)

    tst_results <- shapiro.test(resid_tbl$.resid)

    tst_results <- tibble(variable = 'GLOBAL',
                          type = 'normality',
                          test = 'Shapiro-Wilk test',
                          stat_name = 'W',
                          stat_value = tst_results[['statistic']],
                          df = NA,
                          p_value = tst_results[['p.value']])

    ## proportional hazard

    proph <- data.frame(cox.zph(cox_model)$table)

    proph <- rownames_to_column(proph, 'variable')

    proph <- mutate(proph,
                    type = 'proportional hazard',
                    test = 'zph',
                    stat_name = 'chi-squared',
                    stat_value = chisq,
                    p_value = p)

    ## output

    as_tibble(rbind(tst_results, proph[names(tst_results)]))

  }


# Calibration by D'Agostino - Nam -------

#' Calculate D'Agostino - Nam calibration stats.
#'
#' @description Calculates the D'Agostino - Nam calibration for a Coxph model,
#' a variant of the Hosmer - Lemeshov test (DOI: 10.1016/S0169-7161(03)23001-7),
#' and computes square distances between the predicted survival probability and
#' the outcome as proposed by Graf et al.
#' For  the D'Agostino - Nam calibration the linear predictor score of the Cox
#' model is cut into n quantile strata by \code{\link{cut_quantiles}}.
#' The 'observed' survival odds are derived from a Kaplan-Meier
#' survival estimate calculated by \code{\link[survminer]{surv_fit}}.
#' The 'fitted' survival obs are derived from the Cox proportional hazard models
#' of survival as a function of the linear predictor score quantile intervals.
#'
#' @return an object of class 'calibrator', whose statistic and graphical
#' summary can be accessed by S3 'summary' and 'plot' methods.
#' The object consists of the following elements:
#'
#' * 'lp_scores': a data frame with the linear predictor scores and its
#' quantile intervals;
#' * 'surv_fit' an object of the 'surv_fit' class
#' (\code{\link[survminer]{surv_fit}}), which stores the Kaplan - Meier
#' survival estimates;
#' * 'cox_fit': a coxph model of the survival as a function of
#' the linear predictor score strata;
#' * 'km_estimates': a data frame with the Kaplan-Meier survival estimates
#' for each observation, see: \code{\link[survminer]{surv_summary}} for details;
#' * 'cox_estimates': a data frame storing the number of events and
#' survival probability for each observation in the linear predictor
#' strata Cox model;
#' * 'strata_calibration': a data frame storing the numbers of events and
#' survivors, survival probabilities, relative survival probability
#' (Kaplan - Meier to Cox ratio), relative risk (rr: Kaplan-Meier to Cox ratio)
#' and the D'Agostino-Nam chi-squared statistic (x2_dn) for each strata;
#' * 'global_calibration': a data frame with the mean relative risk,
#' the sum D'Agostino - Nam chi-squared statistic,
#' its degrees of freedom (df = strata number - 2 or df = 1 for two strata)
#' and p value for the global linear predictor score.
#' * 'squares': a data frame with times, observed event, model-predicted
#' survival probability and squared distance between the prediction and observed
#' survival.
#'
#' @param cox_model a CoxpPH model or a coxex object.
#' @param fit a CoxPH model or a coxex object.
#' @param n a single numeric defining the number of quantile intervals of
#' the Cox model's linear predictor score.
#' @param labels an optional user-provided vector of labels for
#' the quantile intervals.
#' @param right logical, indicating if the quantile intervals should be closed
#' on the right (and open on the left) or vice versa.
#' @param use_unique logical, should unique values of the Cox model's linear
#' predictor scores be used instead of quantile strata? This option is
#' non-canonical and experimental, may be utilized to check an association
#' between an ordinal variable and survival.
#' @param ... additional arguments, currently none.
#'
#' @references
#' * D’Agostino, R. B. & Nam, B. H. Evaluation of the Performance of
#' Survival Analysis Models: Discrimination and Calibration Measures.
#' Handb. Stat. 23, 1–25 (2003).
#' * Royston, P. Tools for checking calibration of a Cox model in external
#' validation: Approach based on individual event probabilities.
#' Stata J. 14, 738–755 (2014).
#' * Crowson, C. S. et al. Assessing calibration of prognostic risk
#' scores. Stat. Methods Med. Res. 25, 1692–1706 (2016).
#' * Graf, E., Schmoor, C., Sauerbrei, W. & Schumacher, M. Assessment and
#' comparison of prognostic classification schemes for survival data. Stat.
#' Med. 18, 2529–2545 (1999).
#'
#' @md
#' @export

  get_cox_calibration <- function(cox_model,
                                  n = 3,
                                  labels = NULL,
                                  right = FALSE,
                                  use_unique = FALSE, ...) {

    ## to suppress the R-CMD-check note caused by variables in NSE ------

    strata <- n.event <- .data <- surv <- prob <- NULL
    km_events <- cox_events <- km_prob <- cox_prob <- NULL
    sq_resid <- rr <- x2_dn <- label <- NULL
    df <- NULL

    ## entry control --------

    if(is_coxex(cox_model)) {

      data <- model.frame(cox_model, type = 'data')

      cox_model <- cox_model$model

    }

    if(!inherits(cox_model, 'coxph')) {

      stop('Please provide a valid coxph class model.', call. = FALSE)

    }

    stopifnot(is.numeric(n))
    stopifnot(is.logical(right))
    stopifnot(is.logical(use_unique))

    n <- as.integer(n)

    ## linear predictor score and stratification -------

    lp_score <- predict(cox_model, type = 'lp')

    if(!use_unique) {

      score_tbl <- cut_quantiles(x = lp_score,
                                 n = n,
                                 labels = labels,
                                 right = right)

      if(is.null(labels)) {

        score_tbl <- set_names(score_tbl,
                               c('lp_score',
                                 'strata'))

      } else {

        score_tbl <- set_names(score_tbl,
                               c('lp_score',
                                 'strata',
                                 'label'))

      }

    } else {

      score_tbl <-
        tibble(lp_score = lp_score,
               strata = factor(lp_score))

      if(!is.null(labels)) {

        if(length(labels) != length(levels(score_tbl[['strata']]))) {

          stop(paste('The length of the label vector must correspond',
                     'to the number of unique values of the linear',
                     'predictor score.'),
               call. = FALSE)

        }

        lab_fun <-function(x) {

          set_names(labels, levels(score_tbl[['strata']]))[x]

        }

        score_tbl <-
          mutate(score_tbl,
                 label = fct_relabel(strata, lab_fun))

      }

    }

    score_tbl <- mutate(score_tbl, strata_number = as.numeric(strata))

    surv_object <- model.frame(cox_model)[[1]]

    ## KM estimates ---------

    survfit_obj <- surv_fit(formula = surv_object ~ strata,
                            data = score_tbl)

    km_est <- suppressWarnings(surv_summary(survfit_obj))

    km_est <- as_tibble(km_est)

    ## there's obviously an error with handling some level labels
    ## by survminer::surv_summary. A patch below:

    strata_reco <- set_names(unique(as.character(km_est$strata)),
                             levels(score_tbl$strata))

    km_est <- mutate(km_est,
                     strata = fct_recode(strata,
                                         !!!strata_reco),
                     strata = factor(strata, levels(score_tbl$strata)))

    ## Cox estimates ---------

    cox_model_strata <- coxph(surv_object ~ strata_number,
                              data = score_tbl)

    cox_est <-
      tibble(n.event = predict(cox_model_strata, type = 'expected'),
             strata = score_tbl[['strata']])

    cox_est <-
      mutate(cox_est,
             prob = exp(-n.event),
             time = unclass(model.frame(cox_model_strata)[[1]])[, 1])

    ## mean KM survival and Cox risk estimates per strata ----------

    mean_km_est <- group_by(km_est, .data[['strata']])
    mean_cox_est <- group_by(cox_est, .data[['strata']])

    mean_km_est <- summarise(mean_km_est,
                             km_prob = mean(surv, na.rm = TRUE),
                             km_events = sum(n.event))

    mean_cox_est <- summarise(mean_cox_est,
                              cox_prob = mean(prob, na.rm = TRUE),
                              cox_events = sum(n.event))

    mean_km_est <- mutate(mean_km_est,
                          strata = factor(strata,
                                          levels(score_tbl$strata)))

    mean_cox_est <- mutate(mean_cox_est,
                           strata = factor(strata,
                                           levels(score_tbl$strata)))

    strata_est <- left_join(mean_km_est,
                            mean_cox_est,
                            by = 'strata')

    strata_est <-
      mutate(strata_est,
             n = survfit_obj$n,
             km_survivors = n - km_events,
             cox_survivors = n - cox_events,
             resid = km_prob - cox_prob,
             sq_resid = resid^2,
             rel_prob = km_prob/cox_prob,
             rr = (1 - km_prob)/(1 - cox_prob),
             x2_dn = sq_resid * n/(cox_prob * (1 - cox_prob)))

    global_est <- summarise(strata_est,
                            mean_sq_resid = mean(sq_resid),
                            mean_rr = mean(rr),
                            x2_dn = sum(x2_dn))

    global_est <-
      mutate(global_est,
             df = if(nrow(strata_est) > 2) nrow(strata_est) - 2 else 1,
             p_value = pchisq(q = x2_dn, df = df, lower.tail = FALSE))

    ## the observed and modeled survival without any stratification ------
    ## and squared distance between the observed and predicted survival

    ## the null model will be used to obtain square distance for a reference

    null_model <- coxph(formula = surv_object ~ 1, data = score_tbl)

    preds <- map(list(model = cox_model, reference = null_model),
                 ~tibble(n.event = predict(.x, type = 'expected'),
                         lp_score = predict(.x, type = 'lp'),
                         prob = exp(-n.event)))

    preds <- map(preds,
                 ~mutate(.x, .observation = 1:nrow(.x)))

    surv_df <- surv2df(surv_object)

    preds <- map(preds,
                 left_join,
                 surv_df, by = '.observation')

    preds <- map2(preds, c('squares_model', 'squares_reference'),
                  ~mutate(.x,
                          !!.y := ((1 - prob) - status)^2))

    preds$model <-
      preds$model[c('time', '.observation', 'lp_score',
                    'prob', 'status', 'squares_model')]

    squares <- cbind(preds$model, preds$reference['squares_reference'])

    ## output ----------

    results <- list(lp_scores = score_tbl,
                    surv_fit = survfit_obj,
                    cox_fit = cox_model_strata,
                    km_estimates = km_est,
                    cox_estimates = cox_est,
                    strata_calibration = strata_est,
                    global_calibration = global_est,
                    squares = as_tibble(squares))

    if(!is.null(labels)) {

      label_tbl <- tibble(strata = levels(score_tbl$strata),
                          label = levels(score_tbl$label))

      results[c('km_estimates',
                'cox_estimates',
                'strata_calibration')] <-
        map(results[c('km_estimates',
                      'cox_estimates',
                      'strata_calibration')],
            ~left_join(.x, label_tbl, by = 'strata'))

      results[c('km_estimates',
                'cox_estimates',
                'strata_calibration')] <-
        map(results[c('km_estimates',
                      'cox_estimates',
                      'strata_calibration')],
            mutate,
            strata = factor(strata, levels(score_tbl$strata)),
            label = factor(label, levels(score_tbl$label)))

    }

    eval(call2('calibrator', !!!results))

  }

# Validation statistics --------

#' Validation statistics for 'coxph' and 'coxex' models.
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
#' @param cox_model a 'coxph' or 'coxex' model.
#' @param data modeling data, ignored if a 'coxex' model provided.
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
#' @export

  get_cox_validation <- function(cox_model,
                                 data = NULL,
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

    ## entry control ------

    if(is_coxex(cox_model)) {

      data <- model.frame(cox_model, type = 'data')

      cox_model <- cox_model$model

    }

    if(!inherits(cox_model, 'coxph')) {

      stop('Please provide a valid coxph class model.', call. = FALSE)

    }

    dataset <- NULL
    Dxy <- NULL

    ## fitting a cph model required for validation --------

    cph_model <- cph(formula = formula(cox_model),
                     data = data,
                     x = TRUE,
                     y = TRUE)

    val_results <- validate(cph_model)

    val_results <- t(unclass(as.matrix(val_results)))

    val_results <- as.data.frame(val_results)

    val_results <- rownames_to_column(val_results, 'dataset')

    val_results <- filter(val_results,
                          dataset %in% c('index.orig',
                                         'training',
                                         'test',
                                         'index.corrected'))

    val_results <- mutate(val_results, c_index = (Dxy + 1)/2)

    as_tibble(val_results)

  }

# END -----
