# Plotting functions

#' @include imports.R

  NULL

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

    if(is_coxex(cox_model)) {

      data <- model.frame(cox_model, type = 'data')

      cox_model <- cox_model$model

    }

    if(!inherits(cox_model, 'coxph')) {

      stop('Please provide a valid coxph class model.', call. = FALSE)

    }

    stopifnot(is.theme(cust_theme))

    ## residuals

    qc_tbl <- get_cox_qc(cox_model = cox_model, data = data, ...)

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

    qc_plots <- pmap(qc_plotting_lst,
                     point_plot_,
                     data = qc_tbl,
                     cust_theme = cust_theme)

    set_names(qc_plots, plot_names)

  }

# Fit plots ------

#' Plot outcome and expected survival.
#'
#' @description
#' Plots the actual and Cox model-predicted survival as a
#' Kaplan-Meier plot.
#'
#' @inheritParams get_cox_qc_plots
#' @param palette color palette.
#' @param cust_theme custom ggplot theme.
#' @param ... extra arguments passed to
#' \code{\link[survminer]{ggsurvplot_combine}}.
#'
#' @return a ggplot object.
#'
#' @export

  plot_cox_fit <- function(cox_model,
                           data = NULL,
                           palette = c('steelblue', 'firebrick'),
                           cust_theme = survminer::theme_survminer(), ...) {

    ## suppression of a R-CMD-check note on NSE variables

    raw_rsq <- c_index <- lower_ci <- upper_ci <- NULL

    ## entry control

    if(is_coxex(cox_model)) {

      data <- model.frame(cox_model, type = 'data')

      cox_model <- cox_model$model

    }

    if(!inherits(cox_model, 'coxph')) {

      stop('Please provide a valid coxph class model.', call. = FALSE)

    }

    stopifnot(is.theme(cust_theme))

    ## survival fits

    mod_frame <- model.frame(cox_model)

    surv <- mod_frame[, 1]

    null_formula <- as.formula(surv ~ 1)

    outcome_fit <- surv_fit(null_formula, data = data)

    model_fit <- surv_fit(formula = cox_model, data = data)

    ## meta information

    fit_stats <- get_cox_stats(cox_model = cox_model, data = data)

    fit_stats <- mutate(fit_stats,
                        plot_cap = paste0('R\u00B2 = ', signif(raw_rsq, 2),
                                          ', C = ',
                                          signif(c_index, 2), ' [',
                                          signif(lower_ci, 2), ' - ',
                                          signif(upper_ci, 2), ']'))

    ## plotting

    surv_plot <- ggsurvplot_combine(fit = list(outcome = outcome_fit,
                                               fitted = model_fit),
                                    data = data,
                                    ggtheme = cust_theme,
                                    palette = palette,
                                    legend.labs = c('outcome', 'fitted'),
                                    legend.title = '', ...)

    surv_plot$plot <- surv_plot$plot +
      labs(tag = paste0('total: n = ', fit_stats$n_complete[1],
                        ', events: n = ', fit_stats$n_events[1]),
           subtitle = fit_stats$plot_cap[1]) +
      theme(plot.tag.position = 'bottom')

    surv_plot

  }

# Plots of squared errors --------

#' Plot squared errors between the outcome and observed survival.
#'
#' @description
#' Generates a series of plots of averaged squared errors as a function of
#' unique time points, observations, and unique values of the linear predictor
#' scores.
#'
#' @details
#' Intended for internal use.
#' Mean values are presented as points connected with lines.
#' Optional error region is shown as a ribbon.
#'
#' @param cal_object a 'calibrator' class object.
#' @param palette color palette.
#' @param show_reference logical, should averaged squared errors be displayed
#' for the reference (NULL model)? Defaults to TRUE.
#' @param error_stat error statistic to be presented as a ribbon, defaults
#' to none.
#' @param ribbon_alpha alpha of the ribbon representing the requested error.
#' @param line_width width of the lines connecting the data points.
#' @param point_size size of the data points.
#' @param cust_theme custom ggplot theme.
#' @param ... extra arguments, currently none.
#'
#' @return a list of ggplot objects.

  plot_squares <- function(cal_object,
                           palette = c(reference = 'gray60',
                                       model = 'coral3'),
                           show_reference = TRUE,
                           error_stat = c('none', 'se', '2se', '95perc'),
                           ribbon_alpha = 0.25,
                           line_width = 0.5,
                           point_size = 2,
                           cust_theme = ggplot2::theme_classic(), ...) {

    ## control of the inputs is done by an upstream function

    error_stat <- match.arg(error_stat[1], c('none', 'se', '2se', '95perc'))

    ## plotting data -------

    plot_data <-
      map(c(time = 'time',
            observation = 'observation',
            lp = 'lp'),
          ~summary(cal_object,
                   type = 'squares',
                   by = .x))

    if(show_reference) {

      type <- NULL

      plot_data <-
        map(plot_data,
            function(dat) map(c(model = 'model', reference = 'reference'),
                              ~select(dat,
                                      any_of(c('time', '.observation', 'lp_score')),
                                      ends_with(.x))))

      plot_data <-
        map(plot_data,
            map,
            ~set_names(.x,
                       stri_replace(names(.x),
                                    regex = '_(model|reference)$',
                                    replacement = '')))

      plot_data <-
        map(plot_data,
            ~map2_dfr(.x, names(.x),
                      ~mutate(.x, type = .y)))

    } else {

      plot_data <- map(plot_data,
                       select,
                       any_of(c('time', '.observation', 'lp_score')),
                       ends_with('_model'))

      plot_data <-
        map(plot_data,
            ~set_names(.x,
                       stri_replace(names(.x),
                                    fixed = '_model',
                                    replacement = '')))

      plot_data <- map(plot_data,
                       mutate, type = 'model')

    }

    plot_data <- map(plot_data,
                     mutate,
                     type = factor(type, c('model', 'reference')))

    ## errors to be plotted --------

    lower_error <- upper_error <- NULL
    lower_quant <- upper_quant <- NULL

    if(error_stat == 'none') {

      y_lab <- c('Mean squared error')

    } else if(error_stat == 'se') {

      plot_data <- map(plot_data,
                       mutate,
                       lower_error = mean - sem,
                       upper_error = mean + sem)

      y_lab <- 'Mean squared error \u00B1 SEM'

    } else if(error_stat == '2se') {

      plot_data <- map(plot_data,
                       mutate,
                       lower_error = mean - 2 * sem,
                       upper_error = mean + 2 * sem)

      y_lab <- 'Mean squared error \u00B1 2\u00D7SEM'

    } else {

      plot_data <- map(plot_data,
                       mutate,
                       lower_error = lower_quant,
                       upper_error = upper_quant)

      y_lab <- 'Mean squared error, 95% quantile region'

    }

    plot_data <- map(plot_data,
                     select,
                     any_of(c('time', '.observation', 'lp_score')),
                     type,
                     mean,
                     any_of(c('upper_error', 'lower_error')))

    ## plotting ---------

    sq_plots <- list()

    time <- .observation <- lp_score <- NULL

    sq_plots$time <-
      ggplot(plot_data$time,
             aes(x = time,
                 y = mean,
                 color = type))

    sq_plots$observation <-
      ggplot(filter(plot_data$observation,
                    type == 'model'),
             aes(x = reorder(.observation, mean),
                 y = mean,
                 color = type))

    sq_plots$lp <-
      ggplot(plot_data$lp,
             aes(x = lp_score,
                 y = mean,
                 color = type))

    if(error_stat != 'none') {

      sq_plots <-
        map(sq_plots,
            ~.x +
              geom_ribbon(aes(ymin = lower_error,
                              ymax = upper_error,
                              fill = type,
                              group = type),
                          color = NA,
                          alpha = ribbon_alpha,
                          show.legend = FALSE))

    }

    sq_plots[c('time', 'observation')] <-
      map(sq_plots[c('time', 'observation')] ,
          ~.x +
            geom_line(aes(group = type),
                      linewidth = line_width) +
            geom_point(shape = 16,
                       size = point_size))

    sq_plots$lp <- sq_plots$lp +
      geom_point(shape = 16,
                 size = point_size)

    ## styling of the plots --------

    sq_plots <-
      pmap(list(x = sq_plots,
                y = paste('Squared errors and',
                          c('unique time points',
                            'observations',
                            'linear predictor score')),
                z = c('Survival time',
                      'Sorted observation',
                      'Linear predictor score')),
           function(x, y, z) x +
             scale_color_manual(values = palette,
                                name = '') +
             scale_fill_manual(values = palette,
                               name = '') +
             cust_theme +
             labs(title = y,
                  y = y_lab,
                  x = z))

    sq_plots$observation <- sq_plots$observation +
      theme(axis.text.x = element_blank(),
            axis.ticks.x = element_blank())

    sq_plots

  }

# END ------
