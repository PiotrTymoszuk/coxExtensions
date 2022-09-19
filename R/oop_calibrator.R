# S3 OOP for the class calibrator

# Constructor ------

#' Builds a calibration object.
#'
#' @description Generates a calibration object,
#' used internally by the \code{\link{get_cox_calibration}} function.
#' Contains all data used in calculation of Cox model calibration with the
#' D'Agostino-Nam method (DOI: 10.1016/S0169-7161(03)23001-7) and visualization
#' of the results.
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
#' @return an instance of 'calibrator' class.
#' @export

  calibrator <- function(lp_scores,
                         surv_fit,
                         cox_fit,
                         km_estimates,
                         cox_estimates,
                         strata_calibration,
                         global_calibration) {

    ## entry control

    if(!is.data.frame(lp_scores)) {

      stop('lp_scores must be a data frame.', call. = FALSE)

    }

    if(!all(c('lp_score', 'strata') %in% names(lp_scores))) {

      stop("lp_score needs to contain 'lp_score' and 'strata' variables.",
           call. = FALSE)

    }

    if(!'survfit' %in% class(surv_fit)) {

      stop("surv_fit needs to be an instance of the 'surv_fit' class.",
           call. = TRUE)

    }

    if(!'coxph' %in% class(cox_fit)) {

      stop("cox_fit needs to be a valid 'coxph' model.",
           call. = TRUE)

    }

    classes <- purrr::map_lgl(list(km_estimates,
                                   cox_estimates,
                                   strata_calibration,
                                   global_calibration),
                              is.data.frame)

    if(!all(classes)) {

      stop("The arguments: 'km_estimates', 'cox_estimates', 'strata_calibration' and 'global_calibration' must be data frames.",
           call. = FALSE)

    }

    ## object

    structure(list(lp_scores = lp_scores,
                   surv_fit = surv_fit,
                   cox_fit = cox_fit,
                   km_estimates = km_estimates,
                   cox_estimates = cox_estimates,
                   strata_calibration = strata_calibration,
                   global_calibration = global_calibration),
              class = 'calibrator')

  }

# Class testing ------

#' Test for the calibrator class.
#'
#' @description Tests if an object is an instance of the calibrator class.
#' @param x an object.
#' @return a logical value.
#' @export

  is_calibrator <- function(x) {

    all(class(x) == 'calibrator')

  }

# Appearance ------

#' Print a calibrator object.
#'
#' @description Prints a coxex model.
#' @param x a calibrator object.
#' @param ... extra arguments, none specified.
#' @return none, called fot it's side effects.
#' @export

  print.calibrator <- function(x, ...) {

    stopifnot(coxExtensions::is_calibrator(x))

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

    stopifnot(coxExtensions::is_calibrator(object))

    if('label' %in% names(object$strata_calibration)) {

      return(object$strata_calibration[c('strata', 'label', 'km_events', 'n')])

    } else {

      return(object$strata_calibration[c('strata', 'km_events', 'n')])

    }

  }

# Summary ------

#' Statistic summary of the D'Agostino-Nam calibration.
#'
#' @description Returns a a data frame storing the numbers of events
#' and survivors, survival probabilities,
#' relative survival probability (Kaplan - Meier to Cox ratio),
#' relative risk (rr: Kaplan-Meier to Cox ratio) and
#' the D'Agostino-Nam chi-squared statistic (x2_dn) for each strata
#' (type = 'strata') or for the entire model (type = 'global').
#' @return a data frame with the calibration statistics.
#' @param object a calibrator object.
#' @param type type of the calibration statistic: 'global' (default) or 'strata'.
#' @param ... additional arguments, currently none.
#' @export summary.calibrator
#' @export

  summary.calibrator <- function(object, type = c('global', 'strata'), ...) {

    stopifnot(coxExtensions::is_calibrator(object))

    type <- match.arg(type[1], c('global', 'strata'))

    switch(type,
           global = object$global_calibration,
           strata = object$strata_calibration)

  }

# Plotting ------

#' Plot the observed and fitted survival for the linear predictor score strata.
#'
#' @description Draws a Kaplan-Meier plots of the observed survival (solid line)
#' and the Cox model-predicted survival in each linear predictor score strata.
#' The global calibration measures are presented in plot caption. The legend
#' contains the number of observations assigned to each linear predictor score
#' strata.
#' The total numbers of observations and events are shown in the plot tag.
#' @details The Kaplan-Meier plot is generated with
#' \code{\link[survminer]{ggsurvplot}}.
#' @return a ggplot object.
#' @param x a calibrator object.
#' @param palette a vector of color names corresponding to the strata number.
#' @param plot_title plot title text.
#' @param x_lab X axis label.
#' @param cust_theme custom ggplot theme.
#' @param KM_size size of the Kaplan-Meier plot line.
#' @param cox_size size of the Cox survival plot line.
#' @param cox_linetype type of the Cox survival plot line, dashed by default.
#' @param color_seed seed for the random color generator, ignored
#' if palette provided.
#' @param signif_digits significant digits for rounding the calibration
#' statistics displayed in the plot caption.
#' @export

  plot.calibrator <- function(x,
                              palette = NULL,
                              plot_title = NULL,
                              x_lab = 'Survival',
                              cust_theme = survminer::theme_survminer(),
                              KM_size = 0.5,
                              cox_size = 0.7,
                              cox_linetype = 'dashed',
                              color_seed = 1234,
                              signif_digits = 2) {

    ## entry control

    stopifnot(coxExtensions::is_calibrator(x))

    if(!ggplot2::is.theme(cust_theme)) {

      stop('Please provide a valid ggplot2 theme object.', call. = FALSE)

    }

    stopifnot(is.numeric(KM_size))
    stopifnot(is.numeric(cox_size))

    n_strata <- nrow(x$strata_calibration)

    if(is.null(palette)) {

      if(!is.null(color_seed)) set.seed(color_seed)

      av_colors <- grDevices::colors()

      palette <- sample(av_colors, size = n_strata, replace = FALSE)

    }

    if(length(palette) != n_strata) {

      stop(paste0('Number of colors provided in palette must be equal to the strata number (n = ',
                  n_strata, ')'),
           call. = FALSE)

    }

    ## plot caption, legend labs and tag

    if('label' %in% names(x$strata_calibration)) {

      strata_tags <- paste(x$strata_calibration$label,
                           x$strata_calibration$n,
                           sep = '\nn = ')

    } else {

      strata_tags <- paste(x$strata_calibration$strata,
                           x$strata_calibration$n,
                           sep = '\nn = ')

    }

    strata_tags <- rlang::set_names(strata_tags,
                                    x$strata_calibration$strata)

    plot_tag <- paste0('\ntotal: n = ',
                       sum(x$strata_calibration$n),
                       ', events: n = ',
                       sum(x$strata_calibration$km_events))

    plot_subtitle <- dplyr::mutate(x$global_calibration,
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

    base_plot <- survminer::ggsurvplot(fit = x$surv_fit)

    ## to patch the problem with some strata labels
    ## handled by survminer in a wrong way:

    strata_reco <-
      rlang::set_names(unique(as.character(base_plot$plot$data$strata)),
                       levels(x$lp_scores$strata))

    plot_data <-
      dplyr::mutate(base_plot$plot$data,
                    strata = forcats::fct_recode(strata,
                                                 !!!strata_reco),
                    strata = factor(strata,
                                    levels(x$lp_scores$strata)))

    km_plot <- survminer::ggsurvplot(fit = plot_data,
                                     palette = palette,
                                     legend = 'right',
                                     legend.labs = strata_tags,
                                     title = plot_title,
                                     xlab = x_lab,
                                     ggtheme = cust_theme)

    km_plot <- km_plot +
      ggplot2::labs(subtitle = plot_subtitle,
                    tag = plot_tag)

    ## adding the Cox predictions

    cox_preds <- plyr::dlply(x$cox_estimates[c('prob', 'time', 'strata')],
                             'strata')

    cox_preds <- purrr::map2(cox_preds,
                             names(cox_preds),
                             ~rbind(tibble::tibble(prob = 1,
                                                   time = 0,
                                                   strata = .y),
                                    .x))

    for(i in 1:length(cox_preds)) {

      km_plot <- km_plot +
        ggplot2::geom_line(data = cox_preds[[i]],
                           ggplot2::aes(x = time,
                                        y = prob),
                           color = palette[[i]],
                           show.legend = FALSE,
                           linetype = cox_linetype,
                           size = cox_size)


    }

    return(km_plot)

  }



# END ------
