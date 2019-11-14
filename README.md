
<!-- README.md is generated from README.Rmd. Please edit that file -->
REGRM
=====

<!-- badges: start -->
<!-- badges: end -->
The goal of REGRM is to provide convenient tools for robust estimation of a generalized ratio model, *y*<sub>*i*</sub> = *β**x*<sub>*i*</sub> + *x*<sub>*i*</sub><sup>*γ*</sup>*ε*<sub>*i*</sub>, where *ε* is a homoscedastic error term.

The conventional ratio model *y*<sub>*i*</sub> = *β**x*<sub>*i*</sub> + *ϵ*<sub>*i*</sub> has the heteroscedastic error term *ϵ* and this is the obstacle to have robust estimation by means fo M-estimator.

We re-defined the model with homoscedastic error term *ε*, and then generalized and robustified using the iteratively re-weighted least squares (IRLS) algorithm.

The corresponding estimator is shown as follows:

$$
\\hat{\\beta}\_{rob}=\\frac{\\sum^n\_{i=1}w\_i y\_i x\_i^{1-2\\gamma}}{\\sum^n\_{i=1} w\_i x\_i^{2(1-\\gamma)}},
$$
 where *w*<sub>*i*</sub> is obtained by a weight function. The conventional ratio model is obtained when *γ* = 1/2.

We have two different weight functions, two different scale prameters, and three choices of the gamma value. There are 12 functions according to these settings.

Functions included
------------------

An integrated function REGRM calls an appropriate child function shown in the table internally. REGRM only accepts the tuning parameter c2=4, 6 or 8 and the value will be automatically converted according to the internal 12 founctions.

| Weight function                       | *γ* = 1  | *γ* = 1/2 | *γ* = 0  |
|---------------------------------------|----------|-----------|----------|
| Tukey's biweight function & AAD scale | RrTa.aad | RrTb.aad  | RrTc.aad |
| Tukey's biweight function & MAD scale | RrTa.mad | RrTb.mad  | RrTc.mad |
| Huber's weight function & AAD scale   | RrHa.aad | RrTb.aad  | RrTc.aad |
| Hukey's weight function & MAD scale   | RrHa.mad | RrHb.mad  | RrHc.mad |

If a user wish to select other values of the tuning parameter, any child functions can be called directly. In this case, the tuning parameter c1 should be selected appropriately according to the list below.

-   List of the recommended tuning parameter for REGRM and the 12 functions

| function      | Very robust | ...  | Less robust | Default |
|---------------|-------------|------|-------------|---------|
| REGRM         | 4           | 6    | 8           | 8       |
| Tukey AAD     | 4           | 6    | 8           | 8       |
| Tukey MAD(SD) | 5.01        | 7.52 | 10.03       | 10.03   |
| Huber AAD     | 1.15        | 1.72 | 2.30        | 2.30    |
| Huber MAD(SD) | 1.44        | 2.16 | 2.88        | 2.88    |

Please note Smaller figure provides more robust estimation. As for the function REGRM, tuning parameter c2 accepts only 4, 6 and 8. It will be converted appropriately when REGRAM internally calls an appropriate functions.

Installation
------------

You can install the development version from [GitHub](https://github.com/) with:

``` r
# install.packages("devtools")
devtools::install_github("kazwd2008/REGRM")
```