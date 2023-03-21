# Package tests at the development stage

  library(survival)
  library(coxExtensions)

  survs <- sample(1:102, 500, replace = TRUE)

  events <- as.double(sample(0:1, 500, replace = TRUE))

  test_vars <- rnorm(500)

  test_data <- tibble::tibble(surv_time = sort(survs),
                              event_var = events,
                              var1 = test_vars,
                              var2 = rnorm(500, 15, sd = 5),
                              var3 = c(rep('A', 150),
                                       rep('B', 200),
                                       rep('C', 150)))

  test_uni <- coxph(Surv(surv_time, event_var) ~ var1,
                    data = test_data,
                    x = TRUE)

  test_model <- coxph(Surv(surv_time, event_var) ~ var1 + var2 + var3,
                      data = test_data,
                      x = TRUE)

  test_obj <- coxex(test_model, test_data)

  test_uni <- coxex(test_uni, data = test_data)

  summary(test_uni)


  print(test_obj)

  nobs(test_obj)

  model.frame(test_obj, type = 'surv')

  get_cox_estimates(test_obj, trans_function = exp)

  get_cox_qc(test_obj,
             type.predict = 'survival',
             type.residuals = 'martingale')

  get_cox_qc_plots(test_obj,
                   type.predict = 'survival',
                   type.residuals = 'martingale')

  get_cox_stats(test_obj)

  get_cox_assumptions(test_obj)

  predict(test_obj)

  summary(test_obj, type = 'assumptions')

  plot(test_obj, type = 'residuals')

  plot_cox_fit(test_obj)

  plot(test_obj,
       type = 'fit',
       conf.int = TRUE,
       conf.int.alpha = 0.15)

  predict(test_obj)

  resid(test_obj)

  summary(test_obj, 'fit')

  test_cal <- calibrate(test_obj, n = 4, labels = paste0('Q', 1:4))
  test_uni_cal <- get_cox_calibration(test_uni, n = 3)

  plot(test_cal, show_cox = FALSE)
  plot(test_uni_cal)

  summary(test_cal, 'strata')
  summary(test_uni_cal, 'strata')

  plot(test_obj, 'residuals')

  get_cox_pec(cox_model = test_obj)




