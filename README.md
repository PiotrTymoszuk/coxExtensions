[![R](https://github.com/PiotrTymoszuk/coxExtensions/actions/workflows/r.yml/badge.svg)](https://github.com/PiotrTymoszuk/coxExtensions/actions/workflows/r.yml)

<img src="https://user-images.githubusercontent.com/80723424/226474877-98e6eb4e-daee-495e-90a4-ab110e281e08.png" width="12%" height="12%" align = "right">

# coxExtensions
Extended Inference and Goodness of Fit for Cox Models

## Summary

Extends the toolbox for Cox proportional hazard models (packages `survival`, `SurvMisc`, `rms`) with inference statistic, goodness of fit,  diagnostic functions, calibration tests, and cross-validation pipelines compatible with the tidyverse format.

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

`coxExtensions` uses tools provided by the 
[survival](https://cran.r-project.org/web/packages/survival/index.html), [rms](https://cran.r-project.org/web/packages/rms/index.html), 
[survminer](https://github.com/kassambara/survminer), 
[survMisc](https://cran.r-project.org/web/packages/survMisc/index.html), 
[rlang](https://rlang.r-lib.org/), 
[pec](https://cran.r-project.org/web/packages/pec/index.html), 
[tidyverse](https://www.tidyverse.org/), 
[stringi](https://stringi.gagolewski.com/) and 
[ggrepel](https://github.com/slowkow/ggrepel) packages. 
Many thanks to their Authors and Contributors.
