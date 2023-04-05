# Object oriented programming for the brier class ------

# Constructor and class inheritance --------

#' Build a brier class object
#'
#' @description Generates a `brier` class object based on a list
#' of unique times, Bier scores per timepoint e.g. calculated with
#' \code{\link[pec]{pec}} for the reference, training and test data.
#' @return an instance of the `brier` class with the `plot()` method.
#' The `brier` object is a data frame (`time`, `reference`, `training`
#' and `test` variables) bundling the unique times points with
#' their Brier scores obtained for the reference, training and test data. The
#' `brier` class inherits many of traditional data frame methods, e.g. `filter`
#' or `group_by` provided by the `dplyr` package.
#' @param times a numeric vector of unique time points.
#' @param reference a numeric vector of Brier scores for the reference survival.
#' @param training a numeric vector of Brier scores for the modeled survival in
#' the training dataset.
#' @param test a numeric vector of Brier scores for the modeled survival
#' in the test dataset. Defaults to NULL, which means that no validation errors
#' are provided.
#'
#' @references
#' * Graf, E., Schmoor, C., Sauerbrei, W. & Schumacher, M. Assessment and
#' comparison of prognostic classification schemes for survival data.
#' Stat. Med. 18, 2529–2545 (1999).
#'
#' @md
#' @export

  brier <- function(times,
                    reference,
                    training,
                    test = NULL) {

    ## entry control ------

    if(!is.numeric(times)) {

      stop("'times' must be a numeric vector.", call. = FALSE)

    }

    if(!is.numeric(reference)) {

      stop("'reference' must be a numeric vector.", call. = FALSE)

    }

    if(!is.numeric(training)) {

      stop("'training' must be a numeric vector.", call. = FALSE)

    }

    if(!is.null(test)) {

      if(!is.numeric(test)) {

        stop("'test' must be a numeric vector")

      }

      if(length(test) != length(times)) {

        stop("Lengths of 'test' and 'times' vectors must be equal.",
             call. = FALSE)

      }

    }

    if(length(reference) != length(times) |
       length(training) != length(times)) {

      stop("Lengths of 'reference', 'training' and 'times' must be equal.",
           call. = FALSE)

    }

    ## object structure --------

    brier_obj <-
      tibble::tibble(time = times,
                     reference = reference,
                     training = training)

    if(!is.null(test)) {

      brier_obj <- dplyr::mutate(brier_obj,
                                 test = test)

    } else {

      brier_obj <- dplyr::mutate(brier_obj,
                                 test = NA)

    }

    structure(brier_obj,
              class = c('brier', class(brier_obj)))

  }

#' Test for the brier class.
#'
#' @description Tests if an object is an instance of the `brier` class.
#' @param x an object.
#' @return a logical value.
#' @export

  is_brier <- function(x) {

    inherits(x, 'brier')

  }

# Plotting of the Brier objects -------

#' Plot a 'brier' class object.
#'
#' @description Plots Brier scores as a function of unique time points.
#' @return a single `ggplot` graphics (if `one-plot` is TRUE) or a list of
#' `ggplot` plots for the Brier scores obtained for reference,
#' training and test data each.
#' @param x a \code{\link{brier}} class object.
#' @param one_plot logical, should Brier scores for all datasets
#' be presented in one plot? Defaults to TRUE.
#' @param palette defines colors of the Bier score curves.
#' @param linewidth line size.
#' @param show_reference logical, should the Bier score curve for the
#' reference be plotted? Defaults to TRUE.
#' @param cust_theme custom `ggplot` theme.
#' @param ... extra arguments, currently none.
#' @importFrom rlang .data
#'
#' @references
#' * Graf, E., Schmoor, C., Sauerbrei, W. & Schumacher, M. Assessment and
#' comparison of prognostic classification schemes for survival data.
#' Stat. Med. 18, 2529–2545 (1999).
#'
#' @md
#' @export plot.brier
#' @export

  plot.brier <- function(x,
                         one_plot = TRUE,
                         palette = c(reference = 'gray60',
                                     training = 'steelblue',
                                     test = 'coral3'),
                         linewidth = 0.5,
                         show_reference = TRUE,
                         cust_theme = ggplot2::theme_classic(), ...) {

    ## entry control ------

    stopifnot(is_brier(x))
    stopifnot(is.logical(one_plot))
    stopifnot(is.numeric(linewidth))
    stopifnot(is.logical(show_reference))

    if(!ggplot2::is.theme(cust_theme)) {

      stop("'cust_theme' has to be a valid ggplot theme object.",
           call. = FALSE)

    }



    ## plot metadata -------

    plot_variables <- names(x)[names(x) != 'time']

    if(!show_reference) {

      plot_variables <- plot_variables[plot_variables != 'reference']

    }

    if(any(is.na(x$test))) {

      plot_variables <- plot_variables[plot_variables != 'test']

    }

    brier_score <- NULL
    variable <- NULL

    ## plotting: a single plot -------

    if(one_plot) {

      plot_data <- tidyr::pivot_longer(data = x,
                                       cols = plot_variables,
                                       names_to = 'variable',
                                       values_to = 'brier_score')

      brier_plot <-
        ggplot2::ggplot(plot_data,
                        ggplot2::aes(x = time,
                                     y = brier_score,
                                     color = variable)) +
        ggplot2::geom_path(linewidth = linewidth) +
        ggplot2::scale_color_manual(values = palette,
                                    name = 'Dataset') +
        cust_theme +
        ggplot2::labs(x = 'Time',
                      y = 'Brier score',
                      title = 'Prediction error')

      return(brier_plot)

    }

    ## a list of plots --------

    plot_lst <-
      list(var = plot_variables,
           y = palette[plot_variables],
           z = c(reference = 'Reference',
                 training = 'Training',
                 test = 'Test')[plot_variables])

    plot_lst <-
      purrr::pmap(plot_lst,
                  function(var, y, z) ggplot2::ggplot(x,
                                                      ggplot2::aes(x = time,
                                                                   y = .data[[var]])) +
                    ggplot2::geom_path(linewidth = linewidth,
                                       color = y) +
                    cust_theme +
                    ggplot2::labs(x = 'Time',
                                  y = 'Brier score',
                                  title = z))

    return(rlang::set_names(plot_lst,
                            plot_variables))

  }

# END ------
