
<!-- README.md is generated from README.Rmd. Please edit that file -->

# powertools

<!-- badges: start -->

<!-- badges: end -->

The goal of `powertools` is to provide functions for performing power
and sample size calculations for a variety of different study designs
and outcomes. This package accompanies the book “Power and Sample Size
in R” (2025) by Catherine M. Crespi, CRC Press/Chapman & Hall. This
package provide a variety of functions for sample size and power calculations.

## Installation

The latest stable version of powertools can be installed from CRAN using
`install.packages()`:

``` r
install.packages("powertools")
```

The current development version can be installed using
`devtools::install_github()`:

``` r
devtools::install_github("powerandsamplesize/powertools")
```

## Example

A study hopes to show whether or not a new experimental therapy is
promising. A treatment is considered promising if at least 20% of
participants respond well, and the researchers believe that the true
response proportion is 30%. This is a basic example that shows how to
calculate the sample size needed for this study to achieve 80% power,
following Example 6.1 in the textbook.

``` r
library(powertools)
prop.1samp(N = NULL, p0 = 0.2, pA = 0.3, power = 0.8, sides = 1)
#> [1] 129.8337
```

`prop.1samp` returns `N`, the required sample size for a one proportion
test.
