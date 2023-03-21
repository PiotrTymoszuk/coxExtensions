# Non-exported.

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

    cut_tbl <- tibble::tibble(x = x,
                              strata = cut_vec)

    if(!is.null(labels)) {

      levs <- levels(cut_tbl$strata)

      label_tbl <- tibble::tibble(strata = levs,
                                  label = factor(labels, levels = labels))

      cut_tbl <- dplyr::left_join(cut_tbl,
                                  label_tbl,
                                  by = 'strata')

      cut_tbl <- plyr::mutate(cut_tbl,
                              strata = factor(strata, levels = levs))

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

    dplyr::mutate(data[order(data[[observed]]), ],
                  .expect.norm = stats::qnorm(stats::ppoints(nrow(data))))

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

    data <- dplyr::filter(data, !is.na(.data[[variable]]))

    if(is.numeric(data[[variable]])) {

      tibble::tibble(variable = variable,
                     level = NA_character_,
                     n = nrow(data))

    } else {

      count_tbl <- dplyr::count(data, .data[[variable]])

      count_tbl <- dplyr::mutate(count_tbl,
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

        data <- dplyr::mutate(data, misslab = ifelse(.candidate_missfit == 'yes',
                                                     .rownames,
                                                     NA))

      } else {

        data <- dplyr::mutate(data, misslab = ifelse(.candidate_missfit == 'yes',
                                                     .observation,
                                                     NA))

      }

      fill_colors <- c(no = 'cornflowerblue',
                       yes = 'firebrick4')

      point_plot <- ggplot2::ggplot(data,
                                    ggplot2::aes(x = .data[[x_var]],
                                                 y = .data[[y_var]],
                                                 fill = .candidate_missfit)) +
        ggplot2::geom_point(size = 2,
                            shape = 21) +
        ggrepel::geom_text_repel(ggplot2::aes(label = misslab),
                                 show.legend = FALSE) +
        ggplot2::scale_fill_manual(values = fill_colors,
                                   name = 'Candidate outlier')

    } else {

      point_plot <- ggplot2::ggplot(data,
                                    ggplot2::aes(x = .data[[x_var]],
                                                 y = .data[[y_var]])) +
        ggplot2::geom_point(size = 2,
                            shape = 21,
                            fill = 'steelblue',
                            position = ggplot2::position_jitter(width = point_wjitter,
                                                                height = point_hjitter))

    }

    ## point plot

    point_plot <- point_plot +
      ggplot2::labs(x = x_lab,
                    y = y_lab,
                    title = plot_title) +
      cust_theme

    if(smooth) {

      if(silent) {

        suppressWarnings(point_plot <- point_plot +
                           ggplot2::geom_smooth(show.legend = FALSE,
                                                color = 'black',
                                                fill = 'dodgerblue2', ...))

      } else {

        point_plot <- point_plot +
          ggplot2::geom_smooth(show.legend = FALSE,
                               color = 'black',
                               fill = 'dodgerblue2', ...)

      }

    }

    return(point_plot)

  }

# END ------
