# Non-exported.

# Numeric vector stratification

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
