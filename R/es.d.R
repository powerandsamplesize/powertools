#' Cohen d effect size calculation for one or two means
#'
#' @description
#' This function takes as inputs a difference between two means and a standard deviation
#' and outputs the d effect size as defined by Cohen (1988),
#' also called the standardized mean difference.
#'
#' @details
#' Cohen J (1988) Statistical Power Analysis for the Behavioral Sciences, 2nd edition.
#' Lawrence Erlbaum Associates, Hillsdale, New Jersey
#'
#'
#'
#' @param delta If one mean: muA (the true mean) - mu0 (the mean under the null). If two means: DeltaA (the true difference mu1 - mu2) - Delta0 (the difference under the null).
#' @param sd The estimated standard deviation; defaults to 1.
#'
#' @return A list of the arguments and the d effect size.
#' @export
#'
#' @examples
#' es.d(delta = 6.5 - 5.7, sd = 0.4)

es.d <- function (delta = NULL, sd = 1) {

  # Check if the arguments are specified correctly
  check.param(delta, "req"); check.param(delta, "num")
  check.param(sd, "req"); check.param(sd, "pos")

  # Calculate d
  d <- abs(delta) / sd

  # Print output as a power.htest object
  METHOD <- "Cohen's d effect size calculation for one or two means"
  structure(list(delta = delta, sd = sd, d = d,
                 method = METHOD), class = "power.htest")

}
