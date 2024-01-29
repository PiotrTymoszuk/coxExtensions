# S3 OOP for the class calibrator

#' @include imports.R

  NULL

# Constructor ------

#' Build a calibration object.
#'
#' @description Generates a calibration object,
#' used internally by the \code{\link{get_cox_calibration}} function.
#' Contains all data used in calculation of Cox model calibration with the
#' D'Agostino-Nam method (DOI: 10.1016/S0169-7161(03)23001-7) and visualization
#' of the results.
#'
#' @param lp_scores a data frame containing linear predictor scores.
#' @param surv_fit a 'surv_fit' class instance storing the Kaplan-Meier survival
#' estimates for the linear predictor score strata
#' (see: \code{\link[survminer]{surv_fit}}).
#' @param cox_fit a Cox proportional hazard model of the survival as a function
#' of the linear predictor score strata.
#' @param km_estimates a data frame with Kaplan-Meier survival estimates.
#' @param cox_estimates a data frame with cox survival estimates.
#' @param strata_calibration a data frame with the calibration statistics
#' for each linear predictor score strata.
#' @param global_calibration a data frame with the global
#' calibration statistics.
#' @param squares an optional data frame with the time points, predicted and
#' observed survival and the squared distance between the predicted and observed
#' survival.
#'
#' @return an instance of 'calibrator' class.
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
#'
#' @md
#' @export

  calibrator <- function(lp_scores,
                         surv_fit,
                         cox_fit,
                         km_estimates,
                         cox_estimates,
                         strata_calibration,
                         global_calibration,
                         squares = NULL) {

    ## entry control

    if(!is.data.frame(lp_scores)) {

      stop('lp_scores must be a data frame.', call. = FALSE)

    }

    if(!all(c('lp_score', 'strata') %in% names(lp_scores))) {

      stop("lp_score needs to contain 'lp_score' and 'strata' variables.",
           call. = FALSE)

    }

    if(!inherits(surv_fit, 'survfit')) {

      stop("surv_fit needs to be an instance of the 'surv_fit' class.",
           call. = TRUE)

    }

    if(!inherits(cox_fit, 'coxph')) {

      stop("cox_fit needs to be a valid 'coxph' model.",
           call. = TRUE)

    }

    classes <- map_lgl(list(km_estimates,
                                   cox_estimates,
                                   strata_calibration,
                                   global_calibration),
                              is.data.frame)

    if(!all(classes)) {

      stop(paste("The arguments: 'km_estimates', 'cox_estimates',",
                 "'strata_calibration' and 'global_calibration'",
                 "must be data frames."),
           call. = FALSE)

    }

    ## object

    structure(list(lp_scores = lp_scores,
                   surv_fit = surv_fit,
                   cox_fit = cox_fit,
                   km_estimates = km_estimates,
                   cox_estimates = cox_estimates,
                   strata_calibration = strata_calibration,
                   global_calibration = global_calibration,
                   squares = squares),
              class = 'calibrator')

  }

# Class testing ------

#' Test for the calibrator class.
#'
#' @description Tests if an object is an instance of the calibrator class.
#' @param x an object.
#' @return a logical value.
#' @export

  is_calibrator <- function(x) inherits(x, 'calibrator')

# Appearance ------

#' Print a calibrator object.
#'
#' @description Prints a coxex model.
#' @param x a calibrator object.
#' @param ... extra arguments, none specified.
#' @return none, called for its side effects.
#' @export

  print.calibrator <- function(x, ...) {

    stopifnot(is_calibrator(x))

    print(x$global_calibration)

  }

# Extraction ------

#' Number of observations and events.
#'
#' @description Extracts the number of observations and events in
#' each linear predictor score strata.
#' @param object a calibrator object.
#' @param ... extra arguments, currently none.
#' @return a data frame with the number of total complete observations and
#' events for each linear predictor score strata.
#' @export

  nobs.calibrator<- function(object, ...) {

    stopifnot(is_calibrator(object))

    if('label' %in% names(object$strata_calibration)) {

      return(object$strata_calibration[c('strata', 'label', 'km_events', 'n')])

    } else {

      return(object$strata_calibration[c('strata', 'km_events', 'n')])

    }

  }

# Summary ------

#' Statistic summary of the D'Agostino-Nam calibration and square errors.
#'
#' @description
#' The function returns two types of data:
#'
#' * a data frame storing the numbers of events and survivors,
#' survival probabilities,
#' relative survival probability (Kaplan - Meier to Cox ratio),
#' relative risk (rr: Kaplan-Meier to Cox ratio) and
#' the D'Agostino-Nam chi-squared statistic (x2_dn) for each strata
#' (type = 'strata') or for the entire model (type = 'global').
#'
#' * a data frame with mean, SEM, median, and 95% percentile
#' range for square errors averaged over either unique time points or
#' observations or unique values of the Cox model's linear predictor score.
#'
#' @return a data frame with the calibration statistics.
#'
#' @param object a calibrator object.
#' @param type type of the calibration statistic:
#' 'global' (default) and 'strata' return the output of the D'Agostino - Nam
#' algorithm, while 'squares' returns a data frame with summary statistics of
#' square errors.
#' @param by type of averaging of the square errors by unique time points
#' ('time', default), observations or unique values of the Cox model linear
#' predictor score ('lp').
#' @param ... additional arguments, currently none.
#'
#' @references
#' * D’Agostino, R. B. & Nam, B. H. Evaluation of the Performance of
#' Survival Analysis Models: Discrimination and Calibration Measures.
#' Handb. Stat. 23, 1–25 (2003).
#'
#' * Royston, P. Tools for checking calibration of a Cox model in external
#' validation: Approach based on individual event probabilities.
#' Stata J. 14, 738–755 (2014).
#'
#' * Crowson, C. S. et al. Assessing calibration of prognostic risk
#' scores. Stat. Methods Med. Res. 25, 1692–1706 (2016).
#'
#' * Graf, E., Schmoor, C., Sauerbrei, W. & Schumacher, M. Assessment and
#' comparison of prognostic classification schemes for survival data.
#' Stat. Med. 18, 2529–2545 (1999).
#'
#' @md
#' @export summary.calibrator
#' @export

  summary.calibrator <- function(object,
                                 type = c('global',
                                          'strata',
                                          'squares'),
                                 by = c('time',
                                        'observation',
                                        'lp'), ...) {

    stopifnot(is_calibrator(object))

    type <- match.arg(type[1], c('global', 'strata', 'squares'))

    by <- match.arg(by[1], c('time', 'observation', 'lp'))

    ## D'Agostino - Nam algorithm output

    if(type %in% c('global', 'strata')) {

      return(switch(type,
                    global = object$global_calibration,
                    strata = object$strata_calibration))

    }

    ## square errors

    mean_model <- sem_model <- median_model <- NULL
    lower_quant_model <- upper_quant_model <- NULL
    squares_model <- NULL

    mean_reference <- sem_reference <- median_reference <- NULL
    lower_quant_reference <- upper_quant_reference <- NULL
    squares_reference <- NULL

    by_var <- switch(by,
                     time = 'time',
                     observation = '.observation',
                     lp = 'lp_score')

    summ_tbl <- group_by(object$squares, .data[[by_var]])

    summarise(summ_tbl,
              mean_model = mean(squares_model),
              sem_model = sem(squares_model),
              median_model = median(squares_model,
                                    na.rm = TRUE),
              lower_quant_model = quantile(squares_model, 0.025,
                                           na.rm = TRUE),
              upper_quant_model = quantile(squares_model, 0.975,
                                           na.rm = TRUE),
              mean_reference = mean(squares_reference),
              sem_reference = sem(squares_reference),
              median_reference = median(squares_reference,
                                        na.rm = TRUE),
              lower_quant_reference = quantile(squares_reference, 0.025,
                                               na.rm = TRUE),
              upper_quant_reference = quantile(squares_reference, 0.975,
                                               na.rm = TRUE))

  }

# Plotting ------

#' Plot the observed and fitted survival for the linear predictor score strata.
#'
#' @description Draws two types of plots:
#'
#' * a Kaplan-Meier plots of the observed survival (solid line)
#' and the Cox model-predicted survival in each linear predictor score strata.
#' The global calibration measures are presented in plot caption. The legend
#' contains the number of observations assigned to each linear predictor score
#' strata. The total numbers of observations and events are shown in the plot
#' tag.
#'
#' * a list of plots of square errors averaged over unique time points,
#' observations and unique values of the linear predictor score
#'
#' @details The Kaplan-Meier plot is generated with
#' \code{\link[survminer]{ggsurvplot}}.
#'
#' @return a ggplot object or a list og ggplot objects.
#'
#' @param x a calibrator object.
#' @param type type of the plot or plots:
#'
#' * 'strata' generates a Kaplan-Meier plot of survival in the score strata
#'
#' * 'squares' draws a series of plots of squared errors averaged over unique
#' time points, observations and unique values of the linear predictor score
#'
#' @param palette a vector of color names corresponding to the strata number.
#' @param cust_theme custom ggplot theme.
#' @param KM_size size of the Kaplan-Meier plot line, ignored
#' if `type = 'squares'`.
#' @param show_cox logical, should the Cox model-predicted survival for
#' the strata be plotted? Ignored if `type = 'squares'`.
#' @param cox_size size of the Cox survival plot line.
#' Ignored if `show_cox = FALSE` or `type = 'squares'`.
#' @param cox_linetype type of the Cox survival plot line, dashed by default.
#' Ignored if `show_cox = FALSE` or `type = 'squares'`.
#' @param color_seed seed for the random color generator, ignored
#' if palette provided.
#' @param signif_digits significant digits for rounding the calibration
#' statistics displayed in the plot caption. Ignored if `type = 'squares'`.
#' @param show_reference logical, should averaged squared errors be displayed
#' for the reference (NULL model)? Defaults to TRUE. Used only if
#' `type = 'squares'`.
#' @param error_stat error statistic to be presented as a ribbon, defaults
#' to none. Used only if `type = 'squares'`.
#' @param ... extra arguments passed to \code{\link[survminer]{ggsurvplot}}
#' (if `type = 'strata'`) or to \code{\link{plot_squares}}.
#' @export plot.calibrator
#' @export

  plot.calibrator <- function(x,
                              type = c('strata', 'squares'),
                              palette = NULL,
                              cust_theme = survminer::theme_survminer(),
                              KM_size = 0.5,
                              show_cox = TRUE,
                              cox_size = 0.5,
                              cox_linetype = 'dashed',
                              color_seed = 1234,
                              signif_digits = 2,
                              show_reference = TRUE,
                              error_stat = c('none', 'se', '2se', '95perc'), ...) {

    ## suppression of a R-CMD-check note on NSE variables

    x2_dn <- p_value <- strata <- prob <- NULL
    df <- time <- NULL

    ## entry control --------

    stopifnot(is_calibrator(x))

    type <- match.arg(type[1], c('strata', 'squares'))

    if(!is.theme(cust_theme)) {

      stop('Please provide a valid ggplot2 theme object.', call. = FALSE)

    }

    stopifnot(is.numeric(KM_size))
    stopifnot(is.numeric(cox_size))

    n_strata <- nrow(x$strata_calibration)

    if(is.null(palette)) {

      if(!is.null(color_seed)) set.seed(color_seed)

      av_colors <- grDevices::colors()

      if(type == 'strata') {

        palette <- sample(av_colors, size = n_strata, replace = FALSE)

      } else {

        palette <- c(reference = 'gray60',
                     model = 'coral3')

      }

    }

    if(type == 'strata') {

      if(length(palette) != n_strata) {

        stop(paste0('Number of colors provided in palette must be equal',
                    ' to the strata number (n = ',
                    n_strata, ')'),
             call. = FALSE)

      }

    } else {

      if(length(palette) < 2) {

        stop(paste("The 'palette' argument has to provide names",
                   "for exactly two colors."),
             call. = FALSE)

      }

    }

    stopifnot(is.logical(show_reference))

    error_stat <- match.arg(error_stat[1], c('none', 'se', '2se', '95perc'))

    ## plots of squared errors ---------

    if(type == 'squares') {

      return(plot_squares(cal_object = x,
                          palette = palette,
                          show_reference = show_reference,
                          error_stat = error_stat, ...))

    }

    ## KM plot for the strata ----------

    ## labels and titles

    if('label' %in% names(x$strata_calibration)) {

      strata_tags <- paste(x$strata_calibration$label,
                           x$strata_calibration$n,
                           sep = '\nn = ')

    } else {

      strata_tags <- paste(x$strata_calibration$strata,
                           x$strata_calibration$n,
                           sep = '\nn = ')

    }

    strata_tags <- set_names(strata_tags, x$strata_calibration$strata)

    plot_tag <- paste0('\ntotal: n = ',
                       sum(x$strata_calibration$n),
                       ', events: n = ',
                       sum(x$strata_calibration$km_events))

    plot_subtitle <- mutate(x$global_calibration,
                            plot_lab = paste0('\u03C7\u00B2DN(',
                                              df,
                                              ') = ',
                                              signif(x2_dn,
                                                     signif_digits),
                                              ', p = ',
                                              signif(p_value,
                                                     signif_digits)))

    plot_subtitle <- plot_subtitle$plot_lab

    ## base KM plot

    base_plot <- ggsurvplot(fit = x$surv_fit)

    ## to patch the problem with some strata labels
    ## handled by survminer in a wrong way:

    strata_reco <-
      set_names(unique(as.character(base_plot$plot$data$strata)),
                levels(x$lp_scores$strata))

    plot_data <-
      mutate(base_plot$plot$data,
             strata = fct_recode(strata, !!!strata_reco),
             strata = factor(strata, levels(x$lp_scores$strata)))

    km_plot <- ggsurvplot(fit = plot_data,
                          palette = palette,
                          legend.labs = strata_tags,
                          ggtheme = cust_theme,
                          size = KM_size, ...)

    km_plot <- km_plot +
      labs(subtitle = plot_subtitle,
           tag = plot_tag)

    if(!show_cox) return(km_plot)

    ## adding the Cox predictions

    cox_preds <- split(x$cox_estimates[c('prob', 'time', 'strata')],
                       f = x$cox_estimates[['strata']])

    cox_preds <- map2(cox_preds,
                      names(cox_preds),
                      ~rbind(tibble(prob = 1,
                                    time = 0,
                                    strata = .y),
                             .x))

    for(i in 1:length(cox_preds)) {

      km_plot <- km_plot +
        geom_line(data = cox_preds[[i]],
                  aes(x = time,
                      y = prob),
                  color = palette[[i]],
                  show.legend = FALSE,
                  linetype = cox_linetype,
                  size = cox_size)

    }

    return(km_plot)

  }

# END ------
