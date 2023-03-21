# Functional tools for the CoxPH models.

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

    if(!any(class(cox_model) == 'coxph')) {

      stop('Please provide a valid coxph class model.', call. = FALSE)

    }

    ## meta information

    mod_frame <- model.frame(cox_model)

    modeling_vars <- names(mod_frame)[-1]

    extr_regex <- paste(modeling_vars, collapse = '|')

    ## confidence intervals

    ci_table <- as.data.frame(confint(cox_model))

    ci_table <- rlang::set_names(ci_table,
                                 c('lower_ci', 'upper_ci'))

    ci_table <- purrr::map_dfc(ci_table, trans_function)

    ## coefficients

    coefs <- as.data.frame(summary(cox_model)$coefficients)

    coefs <- rlang::set_names(coefs[, c('coef', 'se(coef)', 'z', 'Pr(>|z|)')],
                              c('estimate', 'se', 'stat', 'p_value'))

    coefs <- tibble::rownames_to_column(coefs, 'parameter')

    coefs <- dplyr::mutate(coefs,
                           estimate = trans_function(estimate),
                           se = trans_function(estimate),
                           n_complete = nrow(model.frame(cox_model)),
                           stat_name = 'z')

    ## output table

    summ_tbl <- tibble::as_tibble(cbind(coefs, ci_table))

    summ_tbl <- dplyr::mutate(summ_tbl,
                              variable = stringi::stri_extract(parameter, regex = extr_regex),
                              level = stringi::stri_replace(parameter, regex = extr_regex, replacement = ''))

    ## counting the n for the

    mod_counts <- purrr::map_dfr(modeling_vars,
                                 ~count_(data = mod_frame, variable = .x))

    mod_counts <- dplyr::mutate(mod_counts,
                                parameter = paste0(variable, level))

    summ_tbl <- dplyr::left_join(summ_tbl,
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

      data <- model.frame.coxex(cox_model, type = 'data')

      cox_model <- cox_model$model

    }

    if(!any(class(cox_model) == 'coxph')) {

      stop('Please provide a valid coxph class model.', call. = FALSE)

    }

    ## computation

    resid_tbl <- broom::augment(x = cox_model,
                                data = data,
                                type.predict = type.predict,
                                type.residuals = type.residuals, ...)

    resid_tbl <- dplyr::mutate(resid_tbl,
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

      data <- model.frame.coxex(cox_model, type = 'data')

      cox_model <- cox_model$model

    }

    if(!any(class(cox_model) == 'coxph')) {

      stop('Please provide a valid coxph class model.', call. = FALSE)

    }

    rsq_type <- match.arg(rsq_type[1], c('mev', 'mer', 'cod'))

    ## R-squared ----

    rsq_tbl <- survMisc::rsq(cox_model)

    ## concordance ----

    c_index <- unname(summary(cox_model)$concordance)

    se <- c_index[2]

    c_tbl <- tibble::tibble(c_index = c_index[1],
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

    res_tbl <- tibble::tibble(n_complete = n_num['total'],
                              n_events = n_num['events'],
                              aic = stats::AIC(cox_model),
                              bic = stats::BIC(cox_model),
                              raw_rsq = rsq_tbl[[rsq_type]],
                              mae = mean(abs(resid_tbl$.resid)),
                              mse = mean(resid_tbl$.resid^2),
                              rmse = sqrt(mean(resid_tbl$.resid^2)))

    tibble::as_tibble(cbind(res_tbl, c_tbl, ibs_tbl))

  }

# Integrated Brier scores -------

#' Calculate prediction error curves and integrated Brier score (IBS).
#'
#' @description Computes prediction error curves and
#' integrated Brier score (IBS)
#' for a Cox model using \code{\link[pec]{pec}} and \code{\link[pec]{crps}}.
#' @return a numeric vactor with IBS values or a `pec` class object.
#' @param cox_model a Cox model of the `coxph` or `coxex` class.
#' @param data the data frame used for the model construction. Ignored,
#' if coxex object provided.
#' @param return_pec logical, should a `pec` object be returned?
#' @param ... extra arguments passed to \code{\link[pec]{pec}}.
#' @import prodlim
#' @export

  get_cox_pec <- function(cox_model,
                          data = NULL,
                          return_pec = FALSE, ...) {

    ## entry control ----

    if(is_coxex(cox_model)) {

      data <- model.frame.coxex(cox_model, type = 'data')

      cox_model <- cox_model$model

    }

    if(!any(class(cox_model) == 'coxph')) {

      stop('Please provide a valid coxph class model.', call. = FALSE)

    }

    ## generation of the pec object ------

    ## KM formula

    km_formula <- paste(as.character(formula(cox_model))[2], '~ 1')

    km_formula <- as.formula(km_formula)

    pec_obj <- pec::pec(object = cox_model,
                        formula = km_formula,
                        data = data)

    if(return_pec) return(pec_obj)

    crps <- pec::crps(pec_obj)

    tibble::tibble(ibs_reference = crps[1, 1],
                   ibs_model = crps[2, 1])

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

      data <- model.frame.coxex(cox_model, type = 'data')

      cox_model <- cox_model$model

    }

    if(!any(class(cox_model) == 'coxph')) {

      stop('Please provide a valid coxph class model.', call. = FALSE)

    }

    ## normality

    resid_tbl <- get_cox_qc(cox_model = cox_model,
                            type.predict = type.predict,
                            type.residuals = type.residuals,
                            data = data, ...)

    tst_results <- shapiro.test(resid_tbl$.resid)

    tst_results <- tibble::tibble(variable = 'GLOBAL',
                                  type = 'normality',
                                  test = 'Shapiro-Wilk test',
                                  stat_name = 'W',
                                  stat_value = tst_results[['statistic']],
                                  df = NA,
                                  p_value = tst_results[['p.value']])

    ## proportional hazard

    proph <- data.frame(survival::cox.zph(cox_model)$table)

    proph <- tibble::rownames_to_column(proph, 'variable')

    proph <- dplyr::mutate(proph,
                           type = 'proportional hazard',
                           test = 'zph',
                           stat_name = 'chi-squared',
                           stat_value = chisq,
                           p_value = p)

    ## output

    tibble::as_tibble(rbind(tst_results,
                            proph[names(tst_results)]))

  }


# Calibration by D'Agostino - Nam -------

#' Calculate D'Agostino - Nam calibration stats.
#'
#' @description Calculates the D'Agostino - Nam calibration for a Coxph model,
#' a variant of the Hosmer - Lemeshov test (DOI: 10.1016/S0169-7161(03)23001-7).
#' The linear predictor score of the Cox model is cut into
#' n quantile strata by \code{\link{cut_quantiles}}.
#' The 'observed' survival odds are derived from a Kaplan-Meier
#' survival estimate calculated by \code{\link[survminer]{surv_fit}}.
#' The 'fitted' survival obs are derived from the Cox proportional hazard models
#' of survival as a function of the linear predictor score quantile intervals.
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
#'
#' @param cox_model a CoxpPH model or a coxex object.
#' @param fit a CoxpPH model or a coxex object.
#' @param n a single numeric defining the number of quantile intervals.
#' @param labels an optional user-provided vector of labels for
#' the quantile intervals.
#' @param right logical, indicating if the quantile intervals should be closed
#' on the right (and open on the left) or vice versa.
#' @param ... additional arguments, currently none.
#' @references
#' * D’Agostino, R. B. & Nam, B. H. Evaluation of the Performance of
#' Survival Analysis Models: Discrimination and Calibration Measures.
#' Handb. Stat. 23, 1–25 (2003).
#' * Royston, P. Tools for checking calibration of a Cox model in external
#' validation: Approach based on individual event probabilities.
#' Stata J. 14, 738–755 (2014).
#' * Crowson, C. S. et al. Assessing calibration of prognostic risk
#' scores. Stat. Methods Med. Res. 25, 1692–1706 (2016).
#' @md
#' @export

  get_cox_calibration <- function(cox_model,
                                  n = 3,
                                  labels = NULL,
                                  right = FALSE) {

    ## to suppress the R-CMD-check note caused by variables in NSE

    strata <- n.event <- .data <- surv <- prob <- NULL
    km_events <- cox_events <- km_prob <- cox_prob <- NULL
    sq_resid <- rr <- x2_dn <- label <- NULL

    ## entry control

    if(is_coxex(cox_model)) {

      data <- model.frame.coxex(cox_model, type = 'data')

      cox_model <- cox_model$model

    }

    if(!any(class(cox_model) == 'coxph')) {

      stop('Please provide a valid coxph class model.', call. = FALSE)

    }

    ## linear predictor score and stratification

    lp_score <- predict(cox_model, type = 'lp')

    score_tbl <- cut_quantiles(x = lp_score,
                               n = n,
                               labels = labels,
                               right = right)

    if(is.null(labels)) {

      score_tbl <- rlang::set_names(score_tbl,
                                    c('lp_score',
                                      'strata'))

    } else {

      score_tbl <- rlang::set_names(score_tbl,
                                    c('lp_score',
                                      'strata',
                                      'label'))

    }

    score_tbl <- dplyr::mutate(score_tbl,
                               strata_number = as.numeric(strata))

    surv_object <- model.frame(cox_model)[[1]]

    ## KM estimates

    survfit_obj <- survminer::surv_fit(formula = surv_object ~ strata,
                                       data = score_tbl)

    km_est <- suppressWarnings(survminer::surv_summary(survfit_obj))

    km_est <- tibble::as_tibble(km_est)

    ## there's obviously an error with handling some level labels
    ## by survminer::surv_summary. A patch below:

    strata_reco <- rlang::set_names(unique(as.character(km_est$strata)),
                                    levels(score_tbl$strata))

    km_est <- dplyr::mutate(km_est,
                            strata = forcats::fct_recode(strata,
                                                         !!!strata_reco),
                            strata = factor(strata, levels(score_tbl$strata)))

    ## Cox estimates

    cox_model_strata <- survival::coxph(surv_object ~ strata_number,
                                       data = score_tbl)

    cox_est <-
      tibble::tibble(n.event = predict(cox_model_strata, type = 'expected'),
                     strata = score_tbl[['strata']])

    cox_est <-
      dplyr::mutate(cox_est,
                    prob = exp(-n.event),
                    time = unclass(model.frame(cox_model_strata)[[1]])[, 1])

    ## mean KM survival and Cox risk estimates per strata

    mean_km_est <- dplyr::group_by(km_est, .data[['strata']])
    mean_cox_est <- dplyr::group_by(cox_est, .data[['strata']])

    mean_km_est <- dplyr::summarise(mean_km_est,
                                    km_prob = mean(surv, na.rm = TRUE),
                                    km_events = sum(n.event))

    mean_cox_est <- dplyr::summarise(mean_cox_est,
                                     cox_prob = mean(prob, na.rm = TRUE),
                                     cox_events = sum(n.event))

    mean_km_est <- dplyr::mutate(mean_km_est,
                                 strata = factor(strata,
                                                 levels(score_tbl$strata)))

    mean_cox_est <- dplyr::mutate(mean_cox_est,
                                  strata = factor(strata,
                                                  levels(score_tbl$strata)))

    strata_est <- dplyr::left_join(mean_km_est,
                                   mean_cox_est,
                                   by = 'strata')

    strata_est <-
      dplyr::mutate(strata_est,
                    n = survfit_obj$n,
                    km_survivors = n - km_events,
                    cox_survivors = n - cox_events,
                    resid = km_prob - cox_prob,
                    sq_resid = resid^2,
                    rel_prob = km_prob/cox_prob,
                    rr = (1 - km_prob)/(1 - cox_prob),
                    x2_dn = sq_resid * n/(cox_prob * (1 - cox_prob)))

    global_est <- dplyr::summarise(strata_est,
                                   mean_sq_resid = mean(sq_resid),
                                   mean_rr = mean(rr),
                                   x2_dn = sum(x2_dn))

    global_est <-
      dplyr::mutate(global_est,
                    df = if(nrow(strata_est) > 2) nrow(strata_est) - 2 else 1,
                    p_value = pchisq(q = x2_dn, df = df, lower.tail = FALSE))

    results <- list(lp_scores = score_tbl,
                    surv_fit = survfit_obj,
                    cox_fit = cox_model_strata,
                    km_estimates = km_est,
                    cox_estimates = cox_est,
                    strata_calibration = strata_est,
                    global_calibration = global_est)

    if(!is.null(labels)) {

      label_tbl <- tibble::tibble(strata = levels(score_tbl$strata),
                                  label = levels(score_tbl$label))

      results[c('km_estimates',
                'cox_estimates',
                'strata_calibration')] <-
        purrr::map(results[c('km_estimates',
                             'cox_estimates',
                             'strata_calibration')],
                  ~dplyr::left_join(.x, label_tbl, by = 'strata'))

      results[c('km_estimates',
                'cox_estimates',
                'strata_calibration')] <-
        purrr::map(results[c('km_estimates',
                             'cox_estimates',
                             'strata_calibration')],
                   dplyr::mutate,
                   strata = factor(strata, levels(score_tbl$strata)),
                   label = factor(label, levels(score_tbl$label)))

    }

    eval(rlang::call2('calibrator', !!!results))

  }

# END -----
