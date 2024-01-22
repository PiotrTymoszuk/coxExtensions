# Exported and non-exported utils.

#' @include imports.R

  NULL

# Numeric vector stratification ------

#' Cut vector in n-quantile intervals.
#'
#' @description A numeric vector is cut in n-intervals corresponding
#' to quantiles (e.g: for n = 3, the vector is cut by tertiles).
#' A data frame is returned which contains the initial value, the interval
#' assignment and, optionally, the user-provided labeling.
#' @details the minimum and maximum values are always included in the first
#' and the last interval, respectively. NAs are silently removed.
#' @return A tibble with the initial values in the order as they were provided,
#' the interval assignment and the user-provided labeling.
#' @param x a numeric vector.
#' @param n a single numeric defining the number of quantile intervals.
#' @param labels an optional user-provided vector of labels.
#' @param right logical, indicating if the intervals should be closed on the
#' right (and open on the left) or vice versa.
#' @export

  cut_quantiles <- function(x, n, labels = NULL, right = FALSE) {

    ## entry control

    if(!is.numeric(x)) {

      stop('The x vector needs to be numeric.', call. = FALSE)

    }

    n <- as.integer(n)

    if(!is.numeric(n) | length(n) > 1) {

      stop('The n paramater needs to be a single integer numeric.',
           call. = FALSE)

    }

    if(!is.null(labels)) {

      if(length(labels) != n) {

        stop('The label vector length needs to be equal to n.', call. = FALSE)

      }

    }

    stopifnot(is.logical(right))

    x <- x[!is.na(x)]

    ## finding the quantiles

    prob_vec <- seq(0, 1, by = 1/n)

    quant_vector <- unname(quantile(x, probs = prob_vec))

    quant_vector[1] <- -Inf
    quant_vector[n + 1] <- Inf

    ## cutting

    strata <- NULL ## to suppress a R-CMD-check note

    cut_vec <- cut(x, breaks = quant_vector)

    cut_tbl <- tibble(x = x, strata = cut_vec)

    if(!is.null(labels)) {

      levs <- levels(cut_tbl$strata)

      label_tbl <- tibble(strata = levs,
                          label = factor(labels, levels = labels))

      cut_tbl <- left_join(cut_tbl,
                           label_tbl,
                           by = 'strata')

      cut_tbl <- mutate(cut_tbl, strata = factor(strata, levels = levs))

    }

    cut_tbl

  }

# Expected values ------

#' Calculate expected values from normal distribution
#'
#' @description Computes expected normal distribution values for a sequence of
#' probability points (\code{\link[stats]{ppoints}})
#' @param data a data frame.
#' @param observed name of the variable storing the observed values.
#' @return the input data frame with the variable `.expect.norm` storing the
#' expected normal distribution variables.

  calc_expected_ <- function(data, observed) {

    mutate(data[order(data[[observed]]), ],
           .expect.norm = qnorm(ppoints(nrow(data))))

  }

# Counting variables -------

#' Count levels of a variable
#'
#' @description Count factor variables or returns the number of complete cases
#' for numeric features.
#' @param data a data frame.
#' @param variable name of the variable of interest.
#' @return a data frame with the counts.

  count_ <- function(data, variable) {

    .data <- NULL

    data <- filter(data, !is.na(.data[[variable]]))

    if(is.numeric(data[[variable]])) {

      tibble(variable = variable,
             level = NA_character_,
             n = nrow(data))

    } else {

      count_tbl <- count(data, .data[[variable]])

      count_tbl <- mutate(count_tbl,
                          level = .data[[variable]],
                          variable = variable)

      count_tbl[c('variable',
                  'level',
                  'n')]

    }

  }

# Scatter plots ------

#' Generate a customized point plot.
#'
#' @description Draws a simple point plot for model diagnostic purposes.
#' @param data data frame.
#' @param x_var name of the variable to be shown in the x axis.
#' @param y_var name of the variable to be shown in the y axis.
#' @param x_lab x axis title.
#' @param y_lab y axis title.
#' @param plot_title plot title.
#' @param plot_subtitle plot subtitle.
#' @param plot_tag plot tag.
#' @param smooth logical, should a trend line be displayed.
#' @param point_wjitter width of the point jittering, defaults to 0.
#' @param point_hjitter height of the point jittering, defaults to 0.
#' @param silent logical, display warnings?
#' @param cust_theme custom ggplot2 theme.
#' @param ... extra arguments passed to geom_smooth().
#' @return a ggplot graphic

  point_plot_ <- function(data, x_var,
                          y_var,
                          x_lab = x_var,
                          y_lab = y_var,
                          plot_title = NULL,
                          plot_subtitle = NULL,
                          plot_tag = NULL,
                          smooth = TRUE,
                          point_wjitter = 0,
                          point_hjitter = 0,
                          silent = TRUE,
                          cust_theme = ggplot2::theme_classic(), ...) {

    ## clearing the issue with variable names and NSE

    .candidate_missfit <- .rownames <- .observation <- .data <- misslab <- NULL

    ## table for plotting

    if('.candidate_missfit' %in% names(data)) {

      if('.rownames' %in% names(data)) {

        data <- mutate(data,
                       misslab = ifelse(.candidate_missfit == 'yes',
                                        .rownames,
                                        NA))

      } else {

        data <- mutate(data,
                       misslab = ifelse(.candidate_missfit == 'yes',
                                        .observation,
                                        NA))

      }

      fill_colors <- c(no = 'cornflowerblue',
                       yes = 'firebrick4')

      point_plot <- ggplot(data,
                           aes(x = .data[[x_var]],
                               y = .data[[y_var]],
                               fill = .candidate_missfit)) +
        geom_point(size = 2,
                   shape = 21) +
        geom_text_repel(aes(label = misslab),
                        show.legend = FALSE) +
        scale_fill_manual(values = fill_colors,
                          name = 'Candidate outlier')

    } else {

      point_plot <- ggplot(data,
                           aes(x = .data[[x_var]],
                               y = .data[[y_var]])) +
        geom_point(size = 2,
                   shape = 21,
                   fill = 'steelblue',
                   position = position_jitter(width = point_wjitter,
                                              height = point_hjitter))

    }

    ## point plot

    point_plot <- point_plot +
      labs(x = x_lab,
           y = y_lab,
           title = plot_title) +
      cust_theme

    if(smooth) {

      if(silent) {

        suppressWarnings(point_plot <- point_plot +
                           geom_smooth(show.legend = FALSE,
                                       color = 'black',
                                       fill = 'dodgerblue2', ...))

      } else {

        point_plot <- point_plot +
          geom_smooth(show.legend = FALSE,
                      color = 'black',
                      fill = 'dodgerblue2', ...)

      }

    }

    return(point_plot)

  }

# Conversion of survival objects to data frames ----------

#' Convert a survival object to a data frame with censored survival data.
#'
#' @description
#' Converts a survival object (class 'Surv') to a data frame with unique time
#' points and event indicators. The function accounts for censoring, i.e.
#' observations censored at a particular time point i do not appear at the
#' time point i + 1.
#'
#' @param surv_object and instance of the 'Surv' class.
#'
#' @return
#' a data frame with the following columns:
#'
#' * 'time': unique time point inferred from the survival object.
#'
#' * '.observation': index of the observation in the input survival object.
#'
#' * 'status': event index as indicated in the survival object
#'
#' @export

  surv2df <- function(surv_object) {

    ## data frame with the survival times and censoring indexes
    ## unique data points and observations

    if(!inherits(surv_object, 'Surv')) {

      stop(!"'surv_object' has to be an instance of the 'Surv' class.",
           call. = FALSE)

    }

    .observation <- NULL
    time <- NULL

    surv_df <- as_tibble(unclass(surv_object))

    surv_df <- mutate(surv_df,
                      .observation = 1:nrow(surv_df))

    unique_times <- sort(unique(surv_df$time))

    ## a list with survival status for unique time points

    surv_res <-
      map(unique_times,
          ~filter(surv_df, time >= .x))

    surv_res <- map2(surv_res, unique_times,
                     ~mutate(.x, status = ifelse(time == .y, status, 0)))

    surv_res <- map(surv_res, ~.x[c('.observation', 'status')])

    surv_res <- map2_dfr(surv_res, unique_times,
                         ~mutate(.x, time = .y))

    relocate(surv_res, time)

  }

# Standard error of the mean ---------

#' Standard error of the mean.
#'
#' @description
#' Computes standard error of the mean (SEM) for a numeric vactor. Any NAs are
#' silently removed.
#'
#' @param x a numeric vector.
#'
#' @return a single numeric value.

  sem <- function(x) {

    x <- x[!is.na(x)]

    sd(x)/sqrt(length(x))

  }

# END ------
