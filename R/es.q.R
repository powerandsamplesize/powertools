#' Cohen q effect size calculation for two correlation coefficients
#'
#' @description
#' Calculates the q effect size for comparing two correlation coefficients. Based on Cohen (1988).
#'
#' @details
#' Cohen J (1988) Statistical Power Analysis for the Behavioral Sciences, 2nd edition.
#' Lawrence Erlbaum Associates, Hillsdale, New Jersey
#'
#'
#' @param rho1 The correlation coefficient in group 1.
#' @param rho2 The correlation coefficient in group 2.
#'
#' @return A list of the arguments and the q effect size.
#' @export
#'
#' @examples
#' es.q(rho1 = 0.3, rho2 = 0.1)

es.q <- function (rho1 = NULL, rho2 = NULL) {

  # Check if the arguments are specified correctly
  check.param(rho1, "req"); check.param(rho1, "unit")
  check.param(rho2, "req"); check.param(rho2, "unit")

  # Calculate q
  rhoprime1 <- 0.5 * log((1 + rho1)/(1 - rho1))
  rhoprime2 <- 0.5 * log((1 + rho2)/(1 - rho2))
  q <- abs(rhoprime1 - rhoprime2)

  # Print output as a power.htest object
  METHOD <- "Cohen q effect size calculation for two correlation coefficients"
  structure(list(rho1 = rho1, rho2 = rho2, q = q,
                 method = METHOD), class = "power.htest")

}
