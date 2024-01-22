# Package S3 generics

#' @include imports.R

  NULL

# Brier scores ----------

#' Calculate Brier scores for survival models.
#'
#' @rdname surv_brier.coxex
#' @export

  surv_brier <- function(fit, ...) UseMethod('surv_brier')

# END ------
