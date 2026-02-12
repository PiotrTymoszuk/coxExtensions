[![R](https://github.com/PiotrTymoszuk/coxExtensions/actions/workflows/r.yml/badge.svg)](https://github.com/PiotrTymoszuk/coxExtensions/actions/workflows/r.yml)

<img src="https://user-images.githubusercontent.com/80723424/226474877-98e6eb4e-daee-495e-90a4-ab110e281e08.png" width="12%" height="12%" align = "right">

# coxExtensions
Extended Inference and Goodness of Fit for Cox Models

## Summary

Extends the toolbox for Cox proportional hazard models (packages `survival`, `SurvMisc`, `rms`) 
with inference statistic, goodness of fit,  diagnostic functions, calibration tests, and 
cross-validation pipelines compatible with the tidyverse format.

## Installation

You may fetch the package easily with `devtools`: 

```r

devtools::install_github('PiotrTymoszuk/coxExtensions')

```

## Terms of use

The package is available under a [GPL-3 license](https://github.com/PiotrTymoszuk/coxExtensions/blob/main/LICENSE).

## Contact

The package maintainer is [Piotr Tymoszuk](mailto:piotr.s.tymoszuk@gmail.com).

## Acknowledgements

`coxExtensions` uses tools provided by the  packages
[survival](https://cran.r-project.org/web/packages/survival/index.html), 
[rms](https://cran.r-project.org/web/packages/rms/index.html), 
[survminer](https://github.com/kassambara/survminer), 
[survMisc](https://cran.r-project.org/web/packages/survMisc/index.html), 
[rlang](https://rlang.r-lib.org/), 
[pec](https://cran.r-project.org/web/packages/pec/index.html), 
[tidyverse](https://www.tidyverse.org/), 
[stringi](https://stringi.gagolewski.com/),  
[ggrepel](https://github.com/slowkow/ggrepel), 
[caret](https://topepo.github.io/caret/), and 
R code for calculation of BCA confidence intervals by 
Jonathan Kropko (https://github.com/jkropko/coxed/blob/master/R/bca.R).

Many thanks to their Authors and Contributors.

## Basic usage: Cox proportional hazard modeling

### Example's data

In the following example of we will model survival in a simulated data set. 
The survival response consists of 0/1-coded survival event and normally distributed 
survival time. 
There are three explanatory variables: normally distributed numeric variables `var1` and 
`var2`, and a categorical variable `var3` with three levels `A`, `B` and `C`. 
Randomly chosen 350 observations in the data set will be used for training of Cox proportional 
hazard models (training subset named `train_data`), while the remaining 150 observations will be 
use only for validation of the models (test subset, `test_data`). 

```r

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

```

```r

> head(train_data)

# A tibble: 6 × 5
  surv_time event_var     var1  var2 var3 
      <int>     <dbl>    <dbl> <dbl> <fct>
1        14         1  0.175   14.0  A    
2        91         0  0.0198  21.9  C    
3       102         0 -0.260   14.5  C    
4       102         0  0.00800 22.3  C    
5       100         1 -0.911    8.82 C    
6        13         1 -0.605   12.2  A    

> head(test_data)

# A tibble: 6 × 5
  surv_time event_var    var1  var2 var3 
      <int>     <dbl>   <dbl> <dbl> <fct>
1        36         0 -1.20   12.0  B    
2        49         0 -1.55   22.5  B    
3        48         0 -0.511   4.62 B    
4        54         0  0.872  11.7  B    
5        48         1 -0.0719 14.5  B    
6        80         1  0.605  17.8  C   

```

### Construction of Cox proportional hazard models

Two models will be trained in the training subset of the data:  
a univariable model with `var1` as the sole explanatory variable, and 
a multi-variable model with explanatory variables `var1`, `var2` and `var3`. 
Below we store the models in a list. 
The function `coxex()` (or `as_coxex()`) from our `coxExtensions` package creates 
a `coxex` class object which stores and the Cox proportional hazard model together with 
the modeling data and enables a seamless use of the package's tools. 
During execution of the code, you may encounter a convergence warning. 
It does not affect the example's results substantially - but in real life cases, 
it's more than recommended to pay attention to such warnings! 
Of importance, it's highly recommended set arguments `x` and `y` to `TRUE` 
during construction of Cox proportional hazard models in R as they enable 
access to the modeling data by multiple tools of the `survival` and `rms` packages.
```r

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

```

```r

> cox_ph_models

$uni
COXEX object with formula:
Surv(surv_time, event_var) ~ var1
<environment: 0x000001d114a69048>

$multi
COXEX object with formula:
Surv(surv_time, event_var) ~ var1 + var2 + var3
<environment: 0x000001d114b6d610>

```

With `summary()` method called for `coxex` models: 

* with `summary(type = "assumptions")`, we can check the proportional 
hazard assumption with the method proposed by Grambsch and Therneau (DOI 10.1093/biomet/81.3.515)

* with `summary(type = "fit")`, we can retrieve fit statistics such as Harrell's concordance index 
(C-index, DOI: 10.1002/(SICI)1097-0258(19960229)15:4<361::AID-SIM168>3.0.CO;2-4 ) or Graf's 
Integrated Brier Score (IBS, DOI: 10.1002/(SICI)1097-0258(19990915/30)18:17/18<2529::AID-SIM274>3.0.CO;2-5) 
for our model and a null/reference model

* with `summary(type = "inference")` we compute model's coefficients with confidence intervals, 
test statistics and p values. 

As you can see in the snippets below, the proportional hazard assumption holds for both the 
univariable and multivariable models. 
Yet the univariable model fits poorly with C = 0.52 and IBS hardly better than the null model. 
By contrast, the multivariable model seems to accurately predict the survival in the training 
data subset; mainly because of the highly significant predictor `var3`.

```r

>   cox_ph_models %>%
+     map(summary, "assumptions")
$uni
# A tibble: 3 × 7
  variable type                test              stat_name   stat_value    df       p_value
  <chr>    <chr>               <chr>             <chr>            <dbl> <dbl>         <dbl>
1 GLOBAL   normality           Shapiro-Wilk test W                0.952    NA 0.00000000284
2 var1     proportional hazard zph               chi-squared      0.806     1 0.369        
3 GLOBAL   proportional hazard zph               chi-squared      0.806     1 0.369        

$multi
# A tibble: 5 × 7
  variable type                test              stat_name   stat_value    df       p_value
  <chr>    <chr>               <chr>             <chr>            <dbl> <dbl>         <dbl>
1 GLOBAL   normality           Shapiro-Wilk test W               0.949     NA 0.00000000118
2 var1     proportional hazard zph               chi-squared     0.0949     1 0.758        
3 var2     proportional hazard zph               chi-squared     1.53       1 0.216        
4 var3     proportional hazard zph               chi-squared     0.663      2 0.718        
5 GLOBAL   proportional hazard zph               chi-squared     2.35       4 0.671  

```
```r

>   cox_ph_models %>%
+     map(summary, "fit")
$uni
# A tibble: 1 × 13
  n_complete n_events   aic   bic raw_rsq   mae   mse  rmse c_index lower_ci upper_ci ibs_reference ibs_model
       <dbl>    <dbl> <dbl> <dbl>   <dbl> <dbl> <dbl> <dbl>   <dbl>    <dbl>    <dbl>         <dbl>     <dbl>
1        350      158 1555. 1558.  0.0017 0.535 0.447 0.669   0.520    0.470    0.571         0.160     0.160

$multi
# A tibble: 1 × 13
  n_complete n_events   aic   bic raw_rsq   mae   mse  rmse c_index lower_ci upper_ci ibs_reference ibs_model
       <dbl>    <dbl> <dbl> <dbl>   <dbl> <dbl> <dbl> <dbl>   <dbl>    <dbl>    <dbl>         <dbl>     <dbl>
1        350      158 1216. 1228.    0.83 0.525 0.435 0.660   0.840    0.818    0.863         0.160    0.0575

``` 

```r

>   cox_ph_models %>%
+     map(summary, "inference")
$uni
# A tibble: 1 × 12
  parameter variable level     n n_complete stat_name  stat estimate     se lower_ci upper_ci p_value
  <chr>     <chr>    <chr> <int>      <int> <chr>     <dbl>    <dbl>  <dbl>    <dbl>    <dbl>   <dbl>
1 var1      var1     ""       NA        350 z         0.657   0.0564 0.0564   -0.112    0.224   0.511

$multi
# A tibble: 4 × 12
  parameter variable level     n n_complete stat_name     stat estimate      se   lower_ci   upper_ci      p_value
  <chr>     <chr>    <chr> <int>      <int> <chr>        <dbl>    <dbl>   <dbl>      <dbl>      <dbl>        <dbl>
1 var1      var1     ""       NA        350 z          0.744     0.0632  0.0632    -0.103     0.230   0.457       
2 var2      var2     ""       NA        350 z         -1.45     -0.0233 -0.0233    -0.0549    0.00826 0.148       
3 var3B     var3     "B"     139        350 z          5.53      5.68    5.68       3.66      7.69    0.0000000319
4 var3A     var3     "A"     104        350 z          0.00988  27.2    27.2    -5372.     5426.      0.992   

```

### Cross-validation of Cox proportional hazard models

Bootstrapping or cross-validation are handy methods to check if and how our models may predict 
survival in an independent data set. 
For details please consult the paper by Simon and colleagues (DOI: 10.1093/bib/bbr001). 
Our package offers two cross-validation interfaces for `coxex` models: 

* Frank Harrell's RMS package validation tool set accessible with `validate()` function

* our native cross-validation interface accessible with `cvCox()` function. Of note it 
enables also repeated cross-validation. Additionally, with `loocvCox()`, we have a possibility 
to assess the model's performance in so called 
[least-one-out cross-validation](https://machinelearningmastery.com/loocv-for-evaluating-machine-learning-algorithms/), 
which is particularly interesting for small size data sets

By calling `validate()` for our models, we compute a bunch of fit statistics for the entire training data set (`training`)
and out-of-fold survival predictions (`test`) including the C-index as a measure of model's accuracy and 
R-square as a metric of explained variance: 

```r

  ### cross-validation with the RMS tool set:
  ### simple 10-fold cross-validation

  set.seed(3131)

  rms_validation_stats <- cox_ph_models %>%
    map(validate,
        method = "crossvalidation",
        B = 10)

```

```r

> rms_validation_stats
$uni
# A tibble: 4 × 9
  dataset            Dxy       R2  Slope         D         U         Q      g c_index
  <chr>            <dbl>    <dbl>  <dbl>     <dbl>     <dbl>     <dbl>  <dbl>   <dbl>
1 index.orig      0.0408  0.00125 1      -0.000365 -0.00129   0.000922 0.0607   0.520
2 training        0.0509  0.00540 1       0.000564 -0.00131   0.00187  0.105    0.525
3 test            0.0245  0.00125 0.0209 -0.000365  0.000747 -0.00111  0.0607   0.512
4 index.corrected 0.0145 -0.00290 0.0209 -0.00129   0.000769 -0.00206  0.0168   0.507

$multi
# A tibble: 4 × 9
  dataset           Dxy    R2 Slope     D        U     Q     g c_index
  <chr>           <dbl> <dbl> <dbl> <dbl>    <dbl> <dbl> <dbl>   <dbl>
1 index.orig      0.681 0.634 1     0.221 -0.00129 0.223  6.88   0.840
2 training        0.692 0.634 1     0.223 -0.00130 0.224  7.59   0.846
3 test            0.680 0.632 0.878 0.220  0.00169 0.218  6.38   0.840
4 index.corrected 0.669 0.632 0.878 0.219  0.00170 0.217  5.68   0.834

```

`cvCox` works in a slightly different manner. 
For each iteration of the cross-validation algorithm, a Cox proportional hazard 
model is fitted for the training portion of the data and linear predictor scores 
are calculated for the test portion. 
The linear predictor scores are then used to fit another Cox proportional hazard 
model in the test portion, where the linear predictor score is the only explanatory 
factor. 
The `cvCox` function returns a list of such "out-of-fold" models, whose performance 
can be easily assessed with `summary(type = "fit")` as illustrated below. 
The disadvantage of this approach is that the "out-of-fold" models with low numbers 
of observations may suffer from convergence problems (we're working on a remedy!). 
The performance statistics are returned only for test portions of the data as means 
with standard deviations and confidence intervals. 

```r
  ### cross-validation with the coxExtension tool set:
  ### 10-repeats 10-fold cross-validation

  set.seed(3131)

  coxex_validation_models <- cox_ph_models %>%
    map(cvCox,
        n_folds = 10,
        n_repeats = 10)
        
  coxex_validation_stats <- coxex_validation_models %>%
    map(summary, type = "fit")

```

```r

>   coxex_validation_stats

$uni
# A tibble: 4 × 10
  summary_stat   aic   bic   raw_rsq    mae    mse   rmse c_index ibs_reference ibs_model
  <chr>        <dbl> <dbl>     <dbl>  <dbl>  <dbl>  <dbl>   <dbl>         <dbl>     <dbl>
1 mean          84.9  85.6 0.0315    0.503  0.402  0.631   0.559         0.154     0.151 
2 SD            16.1  16.3 0.0432    0.0601 0.0772 0.0619  0.0609        0.0167    0.0170
3 lower_ci      56.6  57.0 0.0000536 0.384  0.266  0.516   0.461         0.122     0.118 
4 upper_ci     113.  114.  0.14      0.586  0.532  0.729   0.683         0.183     0.180 

$multi
# A tibble: 4 × 10
  summary_stat   aic   bic raw_rsq    mae    mse   rmse c_index ibs_reference ibs_model
  <chr>        <dbl> <dbl>   <dbl>  <dbl>  <dbl>  <dbl>   <dbl>         <dbl>     <dbl>
1 mean          54.4  55.1  0.783  0.452  0.333  0.574   0.836         0.154    0.0542 
2 SD            11.0  11.2  0.0642 0.0637 0.0664 0.0583  0.0389        0.0186   0.00909
3 lower_ci      32.5  32.7  0.654  0.318  0.207  0.455   0.762         0.120    0.0389 
4 upper_ci      74.6  75.7  0.875  0.575  0.473  0.688   0.910         0.190    0.0715 

```

Independently of the validation interface, we can see that the performance of our 
Cox proportional hazard models in the training data does not differ much from 
cross-validation. 
This means that they likely make generalizable predictions. 
In the next section we will challenge this finding by testing predictions in the test subset 
of the data not used for training of the models.

### Survival predictions

Predictions of `coxex` Cox proportional hazard models, e.g. in form of linear predictor scores, can be 
made with `predict()` method. 
Below, we will compute linear predictor scores for the test subset of our simulated data set and 
use the linear predictor scores to construct Cox proportional hazard models used later for #
evaluation of the predictions. 
These validation models are added to our list storing the models for the test subsets of the data: 

```r

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

```

As above, we assess the proportional hazard assumption and fit statistics 
with `summary()`. 
A prompt inspection of the later shows that the multi-parameter model is highly 
accurate and confident both in the training and the test subset, as inferred from 
high C-indexes and low IBS values:

```r

  ## assumptions and fit stats

  cox_ph_models %>%
    map(summary, "assumptions")

  cox_ph_models %>%
    map(summary, "fit")

```

```r

>   cox_ph_models %>%
+     map(summary, "fit")

$uni
# A tibble: 1 × 13
  n_complete n_events   aic   bic raw_rsq   mae   mse  rmse c_index lower_ci upper_ci ibs_reference ibs_model
       <dbl>    <dbl> <dbl> <dbl>   <dbl> <dbl> <dbl> <dbl>   <dbl>    <dbl>    <dbl>         <dbl>     <dbl>
1        350      158 1555. 1558.  0.0017 0.535 0.447 0.669   0.520    0.470    0.571         0.160     0.160

$multi
# A tibble: 1 × 13
  n_complete n_events   aic   bic raw_rsq   mae   mse  rmse c_index lower_ci upper_ci ibs_reference ibs_model
       <dbl>    <dbl> <dbl> <dbl>   <dbl> <dbl> <dbl> <dbl>   <dbl>    <dbl>    <dbl>         <dbl>     <dbl>
1        350      158 1216. 1228.    0.83 0.525 0.435 0.660   0.840    0.818    0.863         0.160    0.0575

$test_uni
# A tibble: 1 × 13
  n_complete n_events   aic   bic raw_rsq   mae   mse  rmse c_index lower_ci upper_ci ibs_reference ibs_model
       <dbl>    <dbl> <dbl> <dbl>   <dbl> <dbl> <dbl> <dbl>   <dbl>    <dbl>    <dbl>         <dbl>     <dbl>
1        150       66  540.  542.  0.0013 0.529 0.428 0.654   0.501    0.408    0.594         0.162     0.162

$test_multi
# A tibble: 1 × 13
  n_complete n_events   aic   bic raw_rsq   mae   mse  rmse c_index lower_ci upper_ci ibs_reference ibs_model
       <dbl>    <dbl> <dbl> <dbl>   <dbl> <dbl> <dbl> <dbl>   <dbl>    <dbl>    <dbl>         <dbl>     <dbl>
1        150       66  423.  425.    0.75 0.559 0.429 0.655   0.835    0.798    0.872         0.162    0.0695

```

Our package offers another tool to assess model confidence. 
With `surv_brier()` function, Brier scores, i.e. squared differences between the 
predicted risk and 0/1-coded survival event, are computed for all unique time points. 
As such these Brier scores may be interpreted as squared errors of the model and be used 
to check the overall and time point-specific confidence or accuracy of the model. 
For such analysis, a plot of Brier scores for our models and a null model is particularly useful.

```r

  ## confidence: Brier scores for unique time points

  brier_scores <- cox_ph_models %>%
    map(surv_brier)
    
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

```

```r

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

```

### Assessment of calibration with the D'Agostino-Nam method

Finally, we will assess the model's calibration with so called D'Agostino - Nam 
method (DOI: 10.1016/S0169-7161(03)23001-7). 
This approach resembles the general calibration check with the Lemeshov - Hosmer method 
and compares the predicted and actual outcome for quantiles of the linear predictor scores. 
Of note, this method may suffer from several sources of bias - most evidently, it delivers 
different results for varying numbers of quantiles - and is no longer recommended. 
Yet, the methodology allows for a fast visual check of the model fit and "eye-balled" identification 
of possible risk groups.

For our Cox proportional hazard models, we will use the D'Agostino - Nam method to compare 
and visualize survival between quartiles of the linear predictor scores. 
To this end, we will use function `calibrate()` which takes a `coxex` model, 
number of quantiles of the linear predictor score (in our case four), and, optionally, 
labels of the quantiles, which will be displayed in the plots. 

```r

  ## calibration check with the D'Agostino-Nam method:
  ## model fit in the quartiles of the linear predictor scores

  calibrator_obj <- cox_ph_models %>%
    map(calibrate,
        n = 4,
        labels = paste0("Q", 1:4))

```

Kaplan-Meier plots of survival in the quartiles will be generated with `plot()`. 
By default, in the plots, the Kaplan-Meier estimates of survival in the quartiles of 
the linear predictor scores will be displayed as solid lines, predictions of survival 
by the Cox proportional hazard models will be shown as thin dashed lines. 
In the subtitles, statistics of the D'Agostino - Nam tests are presented. 
In particular, high values of the chi-square metric indicate poor concordance of the 
predicted and observed survival in the quantiles of the linear predictor scores, and, 
hence, poor calibration. 
As you can appreciate in the plot panel, the multivariable models, although highly 
accurate and confident, suffers from sub-optimal calibration. 
As inferred from the Kaplan-Meier plots, this phenomenon is primarily attributed 
to poor differences in survival between the second and third quartile of the linear predictor score.

```r

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

```

```r

  calibrator_plots$uni +
    calibrator_plots$test_uni +
    calibrator_plots$multi +
    calibrator_plots$test_multi

```
