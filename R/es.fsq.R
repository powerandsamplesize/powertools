#' Cohen f-squared effect size for overall F test in multiple linear regression
#'
#' @description
#' Computes the f-squared (r^2) effect size for an overall F test in a multiple linear regression model
#' based on the model R^2 (Rsq). Based on Cohen (1988).
#'
#' @details
#' Cohen J (1988) Statistical Power Analysis for the Behavioral Sciences, 2nd edition.
#' Lawrence Erlbaum Associates, Hillsdale, New Jersey
#'
#' @param Rsq The squared sample multiple correlation coefficient.
#'
#' @return A list of the arguments and the f^2 effect size.
#' @export
#'
#' @examples
#' es.fsq(Rsq = 0.02)

es.fsq <- function (Rsq = 0.02) {

  # Check if the arguments are specified correctly
  check.param(Rsq, "req"); check.param(Rsq, "uniti")

  # Calculate f^2
  fsq <- Rsq / (1 - Rsq)

  # Print output as a power.htest object
  METHOD <- "Cohen f^2 effect size calculation for overall F test"
  structure(list(Rsq = Rsq, fsq = fsq,
                 method = METHOD), class = "power.htest")

}
