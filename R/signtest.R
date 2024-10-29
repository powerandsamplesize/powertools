#' Power calculation for sign test
#'
#' @description
#' The sign test is a nonparametric one-sample test of location, specifically,
#' a test of whether the median
#' equals (or is less than or greater than) zero. Its typical use is in place of
#' a one-sample or paired t test when the normality assumption is violated.
#' This function performs power and sample size calculations for a sign test
#' using the normal approximation to the binomial distribution, based on
#' Noether (1987) Sample size determination for some common nonparametric tests. JASA 82(398):645-647.
#'
#' Power calculation
#' for an exact sign test using the exact binomial test can be performed
#' using the power_binom_test function from the MESS package;
#' see Crespi CM (2025) Power and Sample Size in R. Routledge, New York, NY.
#'
#'
#' @details
#' When solving for p, two values, p and 1 - p, are returned.
#' For a two-sided test, due to symmetry, the power for p is equal to the power for 1 - p.
#' For a one-sided upper-tailed test (rejecting null hypothesis when median > 0), select the higher value.
#' For a one-sided lower-tailed test (rejecting null hypothesis when median < 0), select the lower value.
#'
#'
#'
#' @param N The sample size.
#' @param p The probability of a positive difference when the alternative hypothesis is true.
#' @param alpha The significance level (type 1 error rate); defaults to 0.05.
#' @param power The specified level of power.
#' @param sides Either 1 or 2 (default) to specify a one- or two-sided hypothesis test.
#' @param v Either TRUE for verbose output or FALSE (default) to output computed argument only.
#'
#' @return A list of the arguments (including the computed one).
#' @export
#'
#' @examples
#' signtest(N = 40, p = 0.7, power = NULL, alpha = 0.05, sides = 1)


signtest <- function (N = NULL, p = NULL, alpha = 0.05, power = NULL,
                      sides = 2, v = FALSE) {

  # Check if the arguments are specified correctly
  check.many(list(N, alpha, power, p), "oneof")
  check.param(N, "pos"); check.param(N, "min", min = 2)
  check.param(alpha, "unit")
  check.param(power, "unit")
  check.param(p, "unit")
  check.param(sides, "req"); check.param(sides, "vals", valslist = c(1, 2))
  check.param(v, "req"); check.param(v, "bool")

  # power equation
  p.body <- quote({
    stats::pnorm(stats::qnorm(alpha / sides) + 2 * sqrt(N) * abs(p - 0.5))
  })

  # Use stats::uniroot function to calculate missing argument
  if (is.null(power)) {
    power <- eval(p.body)
    if (!v) return(power)
  }
  else if (is.null(N)) {
    N <- stats::uniroot(function(N) eval(p.body) - power, c(2, 1e+07))$root
    if (!v) return(N)
  }
  else if (is.null(p)) {
    p <- stats::uniroot(function(p) eval(p.body) - power, c(0.5, 1 - 1e-10))$root
    p <- c(p, 1 - p)
    if (!v) return(p)
  }
  else if (is.null(alpha)) {
    alpha <- stats::uniroot(function(alpha) eval(p.body) - power, c(1e-10, 1 - 1e-10))$root
    if (!v) return(alpha)
  }
  else stop("internal error")

  # Generate output text
  METHOD <- "One-sample sign test power calculation (normal approximation)"

  # Print output as a power.htest object
  structure(list(N = N, p = p, alpha = alpha,
                 power = power, sides = sides,
                 method = METHOD), class = "power.htest")
}
