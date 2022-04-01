# Functional tools for the CoxPH models.

# Inference ----

#' Get inference statistic for a CoxPH model.
#'
#' @description Retrieves inference statistic for a Cox proportional hazard
#' model. A wrapper around \code{\link[survival]{summary.coxph}}.
#' @return a data frame with the model estimates, confidence intervals
#' and p values.
#' @param cox_model a CoxpPH model or a coxex object.
#' @param transf_fun function used for transformation of the estimates and
#' confidence intervals, identity() by default.
#' @param ... extra arguments, currently none.
#' @export

  get_cox_estimates <- function(cox_model, trans_function = identity, ...) {

    ## entry control

    stopifnot(is.function(trans_function))

    if(coxExtensions::is_coxex(cox_model)) {

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
                                 ~lmqc:::count_(data = mod_frame, variable = .x))

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
#' @export

  get_cox_qc <- function(cox_model,
                         type.predict = 'lp',
                         type.residuals = 'martingale',
                         data = NULL, ...) {

    ## entry control

    if(coxExtensions::is_coxex(cox_model)) {

      data <- coxExtensions::model.frame.coxex(cox_model, type = 'data')

      cox_model <- cox_model$model

    }

    if(!any(class(cox_model) == 'coxph')) {

      stop('Please provide a valid coxph class model.', call. = FALSE)

    }

    ## computation

    resid_tbl <- broom:::augment.coxph(x = cox_model,
                                       data = data,
                                       type.predict = type.predict,
                                       type.residuals = type.residuals, ...)

    resid_tbl <- dplyr::mutate(resid_tbl,
                               .std.resid = scale(.resid)[, 1],
                               .observation = 1:nrow(resid_tbl),
                               .sq.std.resid = .std.resid^2,
                               .candidate_missfit = ifelse(abs(.std.resid) > qnorm(0.975), 'yes', 'no'))

    lmqc:::calc_expected_(resid_tbl, '.std.resid')

  }

# Goodness of fit ------

#' Goodness of fit for a CoxPH model.
#'
#' @description Calculates fit statistics for a CoxPH model: AIC (Akaike
#' information criterion), BIC (Bayesian information criterion), raw_rsq
#' (unadjusted R-squared, cod, mer or mev/default, see:
#' \code{\link[survMisc]{rsq}}), MAE (mean absolute error), MSE
#' (mean squared error), RMSE (root MSE) and C-index with 95% confidence
#' intervals (normality assumption).
#' @return a data frame with the statistic values.
#' @param rsq_type type of R-quared statistic, see: \code{\link[survMisc]{rsq}}.
#' @param ... extra arguments passed to \code{\link{get_cox_qc}}.
#' @inheritParams get_cox_qc
#' @export

  get_cox_stats <- function(cox_model,
                            data = NULL,
                            type.predict = 'lp',
                            type.residuals = 'martingale',
                            rsq_type = c('mev', 'mer', 'cod'), ...) {

    ## entry control

    if(coxExtensions::is_coxex(cox_model)) {

      data <- coxExtensions::model.frame.coxex(cox_model, type = 'data')

      cox_model <- cox_model$model

    }

    if(!any(class(cox_model) == 'coxph')) {

      stop('Please provide a valid coxph class model.', call. = FALSE)

    }

    rsq_type <- match.arg(rsq_type[1], c('mev', 'mer', 'cod'))

    ## R-squared

    rsq_tbl <- survMisc::rsq(cox_model)

    ## concordance

    c_index <- unname(survival:::summary.coxph(cox_model)$concordance)

    se <- c_index[2]

    c_tbl <- tibble::tibble(c_index = c_index[1],
                            lower_ci = c_index[1] + qnorm(0.025) * se,
                            upper_ci = c_index[1] + qnorm(0.975) * se)

    ## n numbers

    surv <- unclass(model.frame(cox_model)[, 1])

    n_num <- c('total' = nrow(surv),
               'events' = sum(surv[, 'status']))

    ## errors

    resid_tbl <- coxExtensions::get_cox_qc(cox_model = cox_model,
                                           type.predict = type.predict,
                                           type.residuals = type.residuals,
                                           data = data, ...)

    ## output

    res_tbl <- tibble::tibble(n_complete = n_num['total'],
                              n_events = n_num['events'],
                              aic = stats::AIC(cox_model),
                              bic = stats::BIC(cox_model),
                              raw_rsq = rsq_tbl[[rsq_type]],
                              mae = mean(abs(resid_tbl$.resid)),
                              mse = mean(resid_tbl$.resid^2),
                              rmse = sqrt(mean(resid_tbl$.resid^2)))

    tibble::as_tibble(cbind(res_tbl, c_tbl))

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
#' @export

  get_cox_assumptions <- function(cox_model,
                                  type.predict = 'lp',
                                  type.residuals = 'martingale',
                                  data = NULL, ...) {

    ## entry controls

    ## entry control

    if(coxExtensions::is_coxex(cox_model)) {

      data <- coxExtensions::model.frame.coxex(cox_model, type = 'data')

      cox_model <- cox_model$model

    }

    if(!any(class(cox_model) == 'coxph')) {

      stop('Please provide a valid coxph class model.', call. = FALSE)

    }

    ## normality

    resid_tbl <- coxExtensions::get_cox_qc(cox_model = cox_model,
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

