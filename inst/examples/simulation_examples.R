# Package tests at the development stage

  library(tidyverse)
  library(trafo)

  library(survival)
  library(coxExtensions)

  library(patchwork)

# simulated testing data -------

  set.seed(1234)

  survs <- sample(1:102, 500, replace = TRUE)

  events <- as.double(sample(0:1, 500, replace = TRUE))

  cox_data <- tibble(surv_time = sort(survs),
                     event_var = events,
                     var1 = rnorm(500),
                     var2 = rnorm(500, 15, sd = 5),
                     var3 = factor(c(rep("A", 150),
                                     rep("B", 200),
                                     rep("C", 150)),
                                   c("C", "B", "A")))

  ## random re-shuffling,
  ## definition of the training and test subset

  set.seed(123243)

  cox_data <- cox_data[sample(1:nrow(cox_data), nrow(cox_data)), ]

  train_data <- cox_data[1:350, ]

  test_data <- cox_data[-1:-350, ]

# Univariable and multi-variable models -------

  ## direct declaration

  cox_ph_models <-
    list(uni = coxph(Surv(surv_time, event_var) ~ var1,
                     data = train_data,
                     x = TRUE,
                     y = TRUE),
         multi = coxph(Surv(surv_time, event_var) ~ var1 + var2 + var3,
                       data = train_data,
                       x = TRUE,
                       y = TRUE)) %>%
    map(coxex, data = train_data)

  ## working with meta-programming via call
  ## to ensure that the data and formula are kept with the model

  cox_ph_models <-
    list(uni = Surv(surv_time, event_var) ~ var1,
         multi = Surv(surv_time, event_var) ~ var1 + var2 + var3) %>%
    map(~call2(.fn = "coxph",
               formula = .x,
               data = train_data,
               x = TRUE,
               y = TRUE)) %>%
    map(eval) %>%
    map(coxex, data = train_data)

# Methods for the coxex class ----------

  ## the proportional hazard assumption,
  ## fit statistics (such as C-index of integrated Brier score)
  ## inference

  cox_ph_models %>%
    map(summary, "assumptions")

  cox_ph_models %>%
    map(summary, "fit")

  cox_ph_models %>%
    map(summary, "inference")

  ## numbers of observations and extraction of modeling data

  cox_ph_models %>%
    map(nobs)

  cox_ph_models %>%
    map(model.frame)

  cox_ph_models %>%
    map(model.frame, type = "surv")

  ## quality control: residuals and QC plots

  cox_ph_models %>%
    map(residuals,
        type.predict = "survival",
        type.residuals = "martingale")

  cox_ph_models %>%
    map(plot,
        type = "fit")

  cox_ph_models %>%
    map(plot,
        type = "residuals",
        type.predict = "survival",
        type.residuals = "martingale")

# cross-validation with the RMS tool set --------

  ### cross-validation with the RMS tool set:
  ### simple 10-fold cross-validation

  set.seed(3131)

  rms_validation_stats <- cox_ph_models %>%
    map(validate,
        method = "crossvalidation",
        B = 10)

  ### cross-validation with the coxExtension tool set:
  ### 10-repeats 10-fold cross-validation

  set.seed(3131)

  coxex_validation_models <- cox_ph_models %>%
    map(cvCox,
        n_folds = 10,
        n_repeats = 10)

  coxex_validation_stats <- coxex_validation_models %>%
    map(summary, type = "fit")

# Prediction of survival -----

  ## linear predictor scores for the test subset of the data
  ## and construction of models to validate the predictions

  test_lp_scores <- cox_ph_models %>%
    map(predict,
        newdata = test_data) %>%
    map(~tibble(lp_score = .x)) %>%
    map(cbind, test_data[, c("surv_time", "event_var")])

  cox_ph_models[c("test_uni", "test_multi")] <-
    list(coxph(formula = Surv(surv_time, event_var) ~ lp_score,
               data = test_lp_scores$uni,
               x = TRUE,
               y = TRUE),
         coxph(formula = Surv(surv_time, event_var) ~ lp_score,
               data = test_lp_scores$multi,
               x = TRUE,
               y = TRUE)) %>%
    map2(test_lp_scores, coxex)

  ## assumptions and fit stats

  cox_ph_models %>%
    map(summary, "assumptions")

  cox_ph_models %>%
    map(summary, "fit")

# Assessment of confidence and calibration ---------

  ## confidence: Brier scores for unique time points

  brier_scores <- cox_ph_models %>%
    map(surv_brier)

  brier_scores %>%
    map(plot)

  ## the Brier scores for the null model are
  ## stored in the `reference` column.
  ## the Brier scores for the model of interest are
  ## stored in the `training` column.
  ## we merge them into a single data frame

  brier_scores <- brier_scores %>%
    map(~.x[c("time", "training")]) %>%
    c(list(reference =
             set_names(brier_scores[[1]][, c("time", "reference")],
                       c("time", "training")))) %>%
    compress(names_to = "model_type") %>%
    mutate(model_type = factor(model_type,
                               c("reference", names(cox_ph_models))))

  ## plot of the Brier scores

  brier_plot <- brier_scores %>%
    ggplot(aes(x = time,
               y = training,
               color = model_type)) +
    geom_path() +
    scale_color_manual(values = c(reference = "gray60",
                                  uni = "firebrick2",
                                  multi = "steelblue2",
                                  test_uni = "firebrick4",
                                  test_multi = "steelblue4"),
                       labels = c(reference = "null model",
                                  uni = "univariable, training",
                                  multi = "multivariable, training",
                                  test_uni = "univairiable, test",
                                  test_multi = "multivariable, training"),
                       name = "model type, data subset") +
    theme_classic() +
    labs(title = "Model confidence",
         x = "survival time",
         y = "model error, Brier score")

  ## calibration check with the D'Agostino-Nam method:
  ## model fit in the quartiles of the linear predictor scores

  calibrator_obj <- cox_ph_models %>%
    map(calibrate,
        n = 4,
        labels = paste0("Q", 1:4))

  calibrator_plots <-
    list(x = calibrator_obj,
         title = c(uni = "Univariable, training",
                   multi = "Multivariable, training",
                   test_uni = "Univairiable, test",
                   test_multi = "Multivariable, training")) %>%
    pmap(plot,
         palette = c("steelblue4",
                     "steelblue2",
                     "firebrick2",
                     "firebrick4"),
         cust_theme = theme_classic() +
           theme(plot.tag.position = "bottom"))

  calibrator_plots$uni +
    calibrator_plots$test_uni +
    calibrator_plots$multi +
    calibrator_plots$test_multi

  ### plots of Brier squares

  square_plots <-
    list(x = calibrator_obj,
         title = c(uni = "Univariable, training",
                   multi = "Multivariable, training",
                   test_uni = "Univairiable, test",
                   test_multi = "Multivariable, training")) %>%
    pmap(plot,
         type = "squares")

# END ---------
