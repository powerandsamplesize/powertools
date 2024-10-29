#' Cohen h effect size calculation for two proportions
#'
#' @description
#' Takes as input the outcome proportions in two groups and returns the h effect size
#' as defined by Cohen (1988).
#'
#' @details
#' Cohen J (1988) Statistical Power Analysis for the Behavioral Sciences, 2nd edition.
#' Lawrence Erlbaum Associates, Hillsdale, New Jersey
#'
#'
#' @param p1 The outcome proportion in group 1.
#' @param p2 The outcome proportion in group 2.
#'
#' @return A list of the arguments and the h effect size.
#' @export
#'
#' @examples
#' es.h(p1 = 0.8, p2 = 0.6)

es.h <- function (p1 = NULL, p2 = NULL) {

  # Check if the arguments are specified correctly
  check.param(p1, "req"); check.param(p1, "unit")
  check.param(p2, "req"); check.param(p2, "unit")

  # Calculate h
  h <- 2 * asin(sqrt(p1)) - 2 * asin(sqrt(p2))

  # Print output as a power.htest object
  METHOD <- "Cohen h effect size calculation for two proportions"
  structure(list(p1 = p1, p2 = p2, h = h,
                 method = METHOD), class = "power.htest")

}
