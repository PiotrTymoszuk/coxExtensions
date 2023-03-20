# Plotting functions

# Residual plots ------

#' Residuals plots for a CoxPH model.
#'
#' @description Generates standard residual plots for a CoxPH model: residuals
#' vs fitted, stabdardized residuals vs fitted, squared standardized residuals
#' vs fitted and a quantile-quantile plot of the standardized residuals.
#' @return a list of ggplot objects.
#' @param cox_model a CoxpPH model or a coxex object.
#' @param cust_theme a custom ggplot theme.
#' @param data the data frame used for the model construction. Ignored,
#' if coxex object provided.
#' @param ... extra arguments passed to \code{\link{get_cox_qc}}.
#' @export

  get_cox_qc_plots <- function(cox_model,
                               data = NULL,
                               cust_theme = ggplot2::theme_classic(), ...) {

    ## entry control

    if(coxExtensions::is_coxex(cox_model)) {

      data <- coxExtensions::model.frame.coxex(cox_model, type = 'data')

      cox_model <- cox_model$model

    }

    if(!any(class(cox_model) == 'coxph')) {

      stop('Please provide a valid coxph class model.', call. = FALSE)

    }

    stopifnot(any(class(cust_theme) == 'theme'))

    ## residuals

    qc_tbl <- coxExtensions::get_cox_qc(cox_model = cox_model,
                                        data = data, ...)

    ## plotting list

    qc_plotting_lst <-
      list(x_var = c('.fitted', '.fitted', '.fitted', '.expect.norm'),
           y_var = c('.resid', '.std.resid', '.sq.std.resid', '.std.resid'),
           plot_title = c('Residuals vs. fitted',
                          'Standardized residuals vs. fitted',
                          'Sqared residuals vs. fitted',
                          'QQ standardized residuals vs expected normal'),
           method = c('loess', 'loess', 'loess', 'lm'),
           smooth = c(TRUE, TRUE, TRUE, TRUE))

    plot_names <- c('resid_fitted',
                    'std.resid_fitted',
                    'sq.resid_fitted',
                    'qq.std.resid')

    ## plotting

    qc_plots <- purrr::pmap(qc_plotting_lst,
                            point_plot_,
                            data = qc_tbl,
                            cust_theme = cust_theme)

    rlang::set_names(qc_plots, plot_names)

  }

# Fit plots ------

#' Plot outcome and expected survival.
#'
#' @description Plots the actual and Cox model-predicted survival as a
#' Kaplan-Meier plot.
#' @inheritParams get_cox_qc_plots
#' @param palette color palette.
#' @param cust_theme custom ggplot theme.
#' @param ... extra arguments passed to
#' \code{\link[survminer]{ggsurvplot_combine}}.
#' @return a ggplot object.
#' @export

  plot_cox_fit <- function(cox_model,
                           data = NULL,
                           palette = c('steelblue', 'firebrick'),
                           cust_theme = survminer::theme_survminer(), ...) {

    ## suppression of a R-CMD-check note on NSE variables

    raw_rsq <- c_index <- lower_ci <- upper_ci <- NULL

    ## entry control

    if(coxExtensions::is_coxex(cox_model)) {

      data <- coxExtensions::model.frame.coxex(cox_model, type = 'data')

      cox_model <- cox_model$model

    }

    if(!any(class(cox_model) == 'coxph')) {

      stop('Please provide a valid coxph class model.', call. = FALSE)

    }

    stopifnot(any(class(cust_theme) == 'theme'))

    ## survival fits

    mod_frame <- model.frame(cox_model)

    surv <- mod_frame[, 1]

    null_formula <- as.formula(surv ~ 1)

    outcome_fit <- survminer::surv_fit(null_formula, data = data)

    model_fit <- survminer::surv_fit(formula = cox_model, data = data)

    ## meta information

    fit_stats <- coxExtensions::get_cox_stats(cox_model = cox_model,
                                              data = data)

    fit_stats <- dplyr::mutate(fit_stats,
                               plot_cap = paste0('R\u00B2 = ', signif(raw_rsq, 2),
                                                 ', C = ',
                                                 signif(c_index, 2), ' [',
                                                 signif(lower_ci, 2), ' - ',
                                                 signif(upper_ci, 2), ']'))

    ## plotting

    surv_plot <- survminer::ggsurvplot_combine(fit = list(outcome = outcome_fit,
                                                          fitted = model_fit),
                                               data = data,
                                               ggtheme = cust_theme,
                                               palette = palette,
                                               legend.labs = c('outcome', 'fitted'),
                                               legend.title = '', ...)

    surv_plot$plot <- surv_plot$plot +
      ggplot2::labs(tag = paste0('total: n = ', fit_stats$n_complete[1],
                                 ', events: n = ', fit_stats$n_events[1]),
                    subtitle = fit_stats$plot_cap[1]) +
      ggplot2::theme(plot.tag.position = 'bottom')

    surv_plot

  }

# END ------
