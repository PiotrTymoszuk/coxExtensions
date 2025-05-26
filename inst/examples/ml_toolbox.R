# Functionality of the ML toolbox for canonical Cox models

# tools ------

  library(survival)
  library(coxExtensions)
  library(tidyverse)
  library(rlang)

# Analysis data -------

  ## the NSW collective serves as a training cohort
  ## for development of Cox models, which are subsequently validated in the
  ## remaining Australian states. Adults only

  mass_aids <- MASS::Aids2

  mass_aids <- mass_aids %>%
    filter(age >= 18) %>%
    mutate(state = factor(state, c('NSW', 'VIC', 'QLD', 'Other')),
           os_days = death - diag,
           death = as.numeric(status) - 1,
           age_strata = cut(age,
                            c(-Inf, 30, 65, Inf),
                            c('young adult', 'middle-aged', 'elderly'))) %>%
    map_dfc(function(x) if(is.factor(x)) droplevels(x) else x) %>%
    filter(complete.cases(.)) %>%
    as_tibble

  mass_aids <- split(mass_aids, mass_aids$state)

# Training the models -------

  ## age-stratified model. to find age-independent markers of overall survival
  ## metaprogramming: to make sure, that the data is always kept with the model

  full_model <-
    call2('coxph',
          formula = Surv(os_days, death) ~ strata(age_strata) + diag + T.categ + sex,
          data = mass_aids$NSW,
          x = TRUE,
          y = TRUE) %>%
    eval %>%
    as_coxex(data = mass_aids$NSW)

  opt_model1 <-
    call2('coxph',
          formula = Surv(os_days, death) ~ strata(age_strata) + diag + T.categ,
          data = mass_aids$NSW,
          x = TRUE,
          y = TRUE) %>%
    eval %>%
    as_coxex(data = mass_aids$NSW)

  opt_model2 <-
    call2('coxph',
          formula = Surv(os_days, death) ~ strata(age_strata) + diag,
          data = mass_aids$NSW,
          x = TRUE,
          y = TRUE) %>%
    eval %>%
    as_coxex(data = mass_aids$NSW)

  null_model <-
    call2('coxph',
          formula = Surv(os_days, death) ~ strata(age_strata),
          data = mass_aids$NSW,
          x = TRUE,
          y = TRUE) %>%
    eval %>%
    as_coxex(data = mass_aids$NSW)

  ## comparison of the nested models by analysis of deviance table
  ## the diagnosis date-only model seems to do the job

  anova(full_model,
        opt_model1,
        opt_model2,
        null_model)

# Assumption testing, fit stats, and inference -----------

  ## assumption test by Grambsch et al.

  opt_model2 %>%
    summary('assumptions')

  ## performance stats

  opt_model2 %>%
    summary('fit')

  ## inference

  opt_model2 %>%
    summary('inference')

  ## cross-validated inference and fit stats

  cox_cv <- cvCox(opt_model2, n_folds = 10, n_repeats = 5)

  cox_cv %>%
    summary('inference')

  cox_cv %>%
    summary('fit', ci_type = 'bca')

  ## last-one-out cross-validation (takes a while)

  cox_loocv <- loocvCox(opt_model2)

  cox_loocv %>%
    summary("inference")

  cox_loocv %>%
    summary("fit")


# Predictions --------

  cox_predictions <-
    mlCox(opt_model2,
          newdata = mass_aids[c("VIC", "QLD", "Other")])

  cox_predictions %>%
    summary('fit')

# END ----
