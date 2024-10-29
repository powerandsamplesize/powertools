#' Power calculation for one-sample proportion test
#'
#' @description
#' Performs power and sample size calculations for a one-sample test of a proportion.
#' Calculations use the large-sample normal approximation to the binomial. Can solve
#' for power, N or alpha.
#'
#'
#' @param N The sample size.
#' @param pA The true proportion.
#' @param p0 The proportion under the null.
#' @param alpha The significance level (type 1 error rate); defaults to 0.05.
#' @param power The specified level of power.
#' @param sides Either 1 or 2 (default) to specify a one- or two- sided hypothesis test.
#' @param v Either TRUE for verbose output or FALSE (default) to output computed argument only.
#'
#' @return A list of the arguments (including the computed one).
#' @export
#'
#' @examples
#' prop.1samp(N = NULL, p0 = 0.2, pA = 0.3, power = 0.8, sides = 1)
#' prop.1samp(N = NULL, p0 = 0.4, pA = 0.5, power = 0.8, sides = 1)
#'

prop.1samp <- function (N = NULL, p0 = NULL, pA = NULL, alpha = 0.05,
                        power = NULL, sides = 2, v = FALSE) {

  # Check if the arguments are specified correctly
  check.many(list(N, alpha, power), "oneof")
  check.param(N, "pos")
  check.param(p0, "req"); check.param(p0, "unit")
  check.param(pA, "req"); check.param(pA, "unit")
  check.param(alpha, "unit")
  check.param(power, "unit")
  check.param(sides, "req"); check.param(sides, "vals", valslist = c(1, 2))
  check.param(v, "req"); check.param(v, "bool")

  # Calculate test statistic
  p.body <- quote({
    d <- abs(pA - p0)
    (stats::pnorm(sqrt(N) * d / sqrt(pA * (1 - pA)) -
                  stats::qnorm(alpha / sides, lower = FALSE)))
  })

  # Use stats::uniroot function to calculate missing argument
  if (is.null(power)) {
    power <- eval(p.body)
    if (!v) return(power)
  }
  else if (is.null(N)) {
    N <- stats::uniroot(function(N) eval(p.body) - power, c(2 + 1e-10, 1e+09))$root
    if (!v) return(N)
  }
  else if (is.null(alpha)) {
    alpha <- stats::uniroot(function(alpha) eval(p.body) - power, c(1e-10, 1 - 1e-10))$root
    if (!v) return(alpha)
  }
  else stop("internal error")

  # Generate output text
  METHOD <- "One sample proportion test power calculation"
  p <- c(p0, pA)

  # Print output as a power.htest object
  structure(list(N = N, `p0, pA` = p, alpha = alpha,
                 power = power, sides = sides,
                 method = METHOD), class = "power.htest")

}

