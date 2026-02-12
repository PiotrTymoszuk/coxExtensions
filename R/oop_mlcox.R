# S3 object-oriented programming interface or the `mlcx` class.

#' @include imports.R

  NULL

# Constructor --------

#' Create a `mlcx` object.
#'
#' @description
#' Creates an instance of `mlcx` class. Technically, it is a list of
#' \code{\link{coxex}} models, storing regression of survival as a function of
#' linear predictor scores in cross-validation folds or external data sets.
#'
#' @details
#' Intended for internal use.
#'
#' @return a `mlcx` object described above with a dedicated
#' \code{\link{summary.mlcx}} method.
#'
#' @param x a list of \code{\link{coxex}} models.
#' @param validation a character specifying validation type: cross-validation
#' (`type = 'cv'`) or external validation (`type = 'external'`). This
#' information will be stored as an attribute of the object.
#' @param ... extra arguments, currently none.

  mlcx <- function(x, validation = c('cv', 'external'), ...) {

    ## input control

    err_text <- "'x' has to be a named list of 'coxex' objects."

    if(!is.list(x)) stop(err_text, call. = FALSE)

    classes <- map_lgl(x, is_coxex)

    if(any(!classes)) stop(err_text, call. = FALSE)

    validation <- match.arg(validation[1], c('cv', 'external'))

    ## output

    structure(x,
              class = 'mlcx',
              validation = validation)

  }

# Class inheritance ------

#' Check inheritance for `mlclass`.
#'
#' @description
#' Checks if an object inherits from `mlcx` class.
#'
#' @param x an object.
#'
#' @return a logical value.
#'
#' @export

  is_mlcx <- function(x) inherits(x, 'mlcx')

# Summary method --------

#' Inference, fit statistics, and assumptions for `mlcx` objects.
#'
#' @description
#' Generates a summary for a `mlcx` object.
#' This can be inference (`type = 'inference'`), fit statistics (`type = 'fit'`)
#' or model assumptions (`type = 'assumptions'`).
#' The results are wrapped in a data frame.
#'
#' @details
#' A wrapper around \code{\link{summary.coxex}}.
#'
#' @return a data frame.
#'
#' @inheritParams summary.coxex
#' @param collapse a logical, should the results be collapsed by
#' cross-validation folds (mean, SD, and 95% confidence intervals)?
#' Refers only to objects generated with \code{\link{cvCox}} and
#' `type = 'inference'` or `type = 'fit'`, and ignored otherwise.
#' @param ci_type type of confidence intervals: percentile (default) or Efron's
#' BCA.
#'
#' @export summary.mlcx
#' @export

  summary.mlcx <- function(object,
                            type = c('inference', 'fit', 'assumptions'),
                            collapse = TRUE,
                            ci_type = c('percentile', 'bca'), ...) {


    stopifnot(is_mlcx(object))

    type <- match.arg(type[1], c('inference', 'fit', 'assumptions'))

    ci_type <- match.arg(ci_type[1], c('percentile', 'bca'))

    output <- map(object, summary, type)

    if(type == 'assumptions') {

      variable <- NULL

      output <- map(output, filter, variable == 'GLOBAL')

    }

    if(attr(object, 'validation') == 'external') {

      out_name <- 'data_set'

    } else {

      out_name <- 'resample'

    }

    output <- map2_dfr(output, names(output),
                       ~mutate(.x, !!out_name := .y))

    output <- relocate(output, .data[[out_name]])

    if(attr(object, 'validation') == 'external') return(output)

    if(type == 'assumptions') return(output)

    if(!collapse) return(output)

    output <-
      output[!names(output) %in% c('data_set', 'resample', 'stat_name',
                                   'parameter', 'variable', 'level',
                                   'n', 'n_complete', 'n_events',
                                   'lower_ci', 'upper_ci', 'p_value', 'se')]

    ci_fun <-
      switch(ci_type,
             percentile = function(x) quantile(x, c(0.025, 0.975), na.rm = TRUE),
             bca = function(x) bca(na.omit(x)))

    means <- map_dbl(output, mean, na.rm = TRUE)
    sds <- map_dbl(output, sd, na.rm = TRUE)
    lower_ci <- map_dbl(output, function(x) ci_fun(x)[1])
    upper_ci <- map_dbl(output, function(x) ci_fun(x)[2])

    summ_tbl <- as.data.frame(rbind(means, sds, lower_ci, upper_ci))

    summ_tbl <- set_names(summ_tbl, names(output))

    summary_stat <- NULL

    summ_tbl <- mutate(summ_tbl,
                       summary_stat = c('mean', 'SD', 'lower_ci', 'upper_ci'))

    as_tibble(relocate(summ_tbl, summary_stat))

  }

# END ------
