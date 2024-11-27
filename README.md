
<!-- README.md is generated from README.Rmd. Please edit that file -->

# powertools

<!-- badges: start -->

[![R-CMD-check](https://github.com/powerandsamplesize/powertools/actions/workflows/R-CMD-check.yaml/badge.svg)](https://github.com/powerandsamplesize/powertools/actions/workflows/R-CMD-check.yaml)
<!-- badges: end -->

Power and sample size calculations for a variety of study designs and
outcomes. Methods include t tests, ANOVA (including tests for
interactions, simple effects and contrasts), proportions, categorical
data (chi-square tests and proportional odds), linear, logistic and
Poisson regression, alternative and coprimary endpoints, power for
confidence intervals, correlation coefficient tests, cluster randomized
trials, individually randomized group treatment trials, multisite
trials, treatment-by-covariate interaction effects and nonparametric
tests of location. Utilities are provided for computing various effect
sizes.

This package is a companion to the book, Crespi (2025) Power and Sample
Size in R, Chapman and Hall/CRC. That book uses this package as well as
many other R packages.

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
