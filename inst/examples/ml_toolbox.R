# Functionality of the ML toolbox for canonical Cox models:
#
# Cox proportional hazard models of AIDS mortality in Australia's states.
# The models are trained with the New south Wales data and validated in the
# remaining states.

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
    mutate(state = factor(state, c("NSW", "VIC", "QLD", "Other")),
           os_days = death - diag,
           death = as.numeric(status) - 1,
           age_strata = cut(age,
                            c(-Inf, 30, 65, Inf),
                            c("young adult", "middle-aged", "elderly"))) %>%
    map_dfc(function(x) if(is.factor(x)) droplevels(x) else x) %>%
    filter(complete.cases(.)) %>%
    as_tibble

  mass_aids <- split(mass_aids, mass_aids$state)

# Training the models -------

  ## the models are trained with the AIDS data in New South Wales (NSW).
  ## A series of models is constructed with nestdd sets of predictor,
  ## the models are age-stratifies to find age-independent markers of overall survival
  ##
  ## The models are constructed metaprogramming: to make sure, that
  ## the data is always kept with the model

  aids_models <-
    list(full = Surv(os_days, death) ~ strata(age_strata) + diag + T.categ + sex,
         option1 = Surv(os_days, death) ~ strata(age_strata) + diag + T.categ,
         option2 = Surv(os_days, death) ~ strata(age_strata) + diag,
         null = Surv(os_days, death) ~ strata(age_strata)) %>%
    map(~call2(.fn = "coxph",
               formula = .x,
               data = mass_aids$NSW,
               x = TRUE,
               y = TRUE)) %>%
    map(eval) %>%
    map(as_coxex, data = mass_aids$NSW)

  ## comparison of the nested models by analysis of deviance table

  anova(aids_models$full, aids_models$null)

  anova(aids_models$full,
        aids_models$option1,
        aids_models$option2,
        aids_models$null) ## all non-null models add to the prediction

  ## the full model is analyzed further

  full_models <- list()

  full_models$NSW <- aids_models$full

# The full model: assumption testing, fit stats, and inference -----------

  ## assumption test by Grambsch et al.: potentially violated by
  ## the `T.categ` variable

  full_models$NSW %>%
    summary("assumptions")

  ## performance stats

  full_models$NSW %>%
    summary("fit")

  ## inference

  full_models$NSW %>%
    summary("inference")

# Cross-validation of the full model --------

  ## cross-validated fit stats

  full_models_cv <- cvCox(full_models$NSW,
                          n_folds = 10,
                          n_repeats = 5)

  full_models_cv %>%
    summary("fit", ci_type = "bca")

  ## last-one-out cross-validation (takes a while)

  full_models_loocv <- loocvCox(full_models$NSW)

  full_models_loocv %>%
    summary("fit", ci_type = "bca")

# Predictions of survival for the remaining states --------

  full_models[c("NSW", "VIC", "QLD", "Other")] <-
    mlCox(full_models$NSW,
          newdata = mass_aids[c("VIC", "QLD", "Other")])

  ## assumptions and fit stats

  full_models %>%
    map(summary, "assumptions")

  full_models %>%
    map(summary, "fit")

# END ----
