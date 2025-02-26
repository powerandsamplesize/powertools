#' Power calculation for signed-rank test
#'
#' @description
#' The signed-rank test is a nonparametric alternative to a one-sample or paired t test. This function
#' performs power and sample size calculations for the signed-rank test using Noether's approximation;
#' see Noether (1987) Sample size determination for some common nonparametric tests. JASA 82(398):645-647.
#'
#' @details
#' Due to symmetry, the power for p is equal to the power for 1 - p. Therefore,
#' when solving for p, two values, p and 1 - p, are returned.
#'
#'
#' @param N The sample size; number of observations or paired differences.
#' @param ps The probability that the sum of two values exceeds zero when the alternative hypothesis is true.
#' @param alpha The significance level or type 1 error rate; defaults to 0.05.
#' @param power The specified level of power.
#' @param sides Either 1 or 2 (default) to specify a one- or two-sided hypothesis test.
#' @param v Either TRUE for verbose output or FALSE to output computed argument only.
#'
#' @return A list of the arguments (including the computed one).
#' @export
#'
#' @examples
#' signedrank(N = 20, ps = 0.87, power = NULL, sides = 2)


signedrank <- function (N = NULL, ps = NULL, alpha = 0.05, power = NULL,
                        sides = 2, v = FALSE) {

  # Check if the arguments are specified correctly
  check.many(list(N, alpha, power, ps), "oneof")
  check.param(N, "pos"); check.param(N, "min", min = 2)
  check.param(alpha, "unit")
  check.param(power, "unit")
  check.param(ps, "unit")
  check.param(sides, "req"); check.param(sides, "vals", valslist = c(1, 2))
  check.param(v, "req"); check.param(v, "bool")

  # power equation
  p.body <- quote({
    stats::pnorm(stats::qnorm(alpha / sides) + sqrt(3) * sqrt(N) * abs(ps - 0.5))
  })

  # Use safe.uniroot function to calculate missing argument
  if (is.null(power)) {
    power <- eval(p.body)
    if (!v) return(power)
  }
  else if (is.null(N)) {
    N <- safe.uniroot(function(N) eval(p.body) - power, c(2, 1e+07))$root
    if (!v) return(N)
  }
  else if (is.null(ps)) {
    ps <- safe.uniroot(function(ps) eval(p.body) - power, c(0.5, 1 - 1e-10))$root
    ps <- c(ps, 1 - ps)
    if (!v) return(ps)
  }
  else if (is.null(alpha)) {
    alpha <- safe.uniroot(function(alpha) eval(p.body) - power, c(1e-10, 1 - 1e-10))$root
    if (!v) return(alpha)
  }
  else stop("internal error")

  # Generate output text
  METHOD <- "Signed-rank test power calculation (normal approximation)"

  # Print output as a power.htest object
  structure(list(N = N, ps = ps, alpha = alpha,
                 power = power, sides = sides,
                 method = METHOD), class = "power.htest")
}
