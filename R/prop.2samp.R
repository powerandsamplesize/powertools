#' Power calculation for two-sample proportion test
#'
#' @description
#' Performs power and sample size calculations for two-sample tests of proportions
#' using normal approximation to the binomial. Can solve for power, n1, n.ratio
#' or alpha.
#'
#' @details
#' For a noninferiority or superiority by a margin test, a one-sided test should be used. See Crespi (2025)
#' for more guidance. For an equivalence test for two proportions, see the prop.test.equiv.
#'
#'
#' @param n1 The sample size for group 1.
#' @param n.ratio The ratio n2/n1 between the sample sizes of two groups; defaults to 1 (equal group sizes).
#' @param p1 The proportion in group 1.
#' @param p2 The proportion in group 2.
#' @param margin The margin of noninferiority or superiority; defaults to 0.
#' @param alpha The significance level (type 1 error rate); defaults to 0.05.
#' @param power The specified level of power.
#' @param sides Either 1 or 2 (default) to specify a one- or two- sided hypothesis test.
#' @param v Either TRUE for verbose output or FALSE (default) to output computed argument only.
#'
#' @return A list of the arguments (including the computed one).
#' @export
#'
#' @examples
#' prop.2samp(n1 = NULL, p1 = 0.6, p2 = 0.8, alpha = 0.025, power = 0.9, sides = 1)
#' prop.2samp(n1 = NULL, p1 = 0.25, p2 = 0.25, margin = 0.1, alpha = 0.025, power = 0.8, sides = 1)

prop.2samp <- function (n1 = NULL, n.ratio = 1, p1 = NULL, p2 = NULL, margin = 0,
                        alpha = 0.05, power = NULL, sides = 2, v = FALSE) {

  # Check if the arguments are specified correctly
  check.many(list(n1, n.ratio, alpha, power), "oneof")
  check.param(n1, "pos")
  check.param(n.ratio, "pos")
  check.param(p1, "req"); check.param(p1, "unit")
  check.param(p2, "req"); check.param(p2, "unit")
  check.param(margin, "req"); check.param(margin, "num")
  check.param(alpha, "unit")
  check.param(power, "unit")
  check.param(sides, "req"); check.param(sides, "vals", valslist = c(1, 2))
  check.param(v, "req"); check.param(v, "bool")

  # Calculate n1
  p.body <- quote({
    d <- abs(p1 - p2) - margin
    q1 <- 1 - p1
    q2 <- 1 - p2
    ((stats::qnorm(alpha / sides) + stats::qnorm(1 - power))^2 *
    (n.ratio * p1 * q1 + p2 * q2) / (n.ratio * d^2))
  })

  # Use stats::uniroot function to calculate missing argument
  if (is.null(n1)) {
    n1 <- eval(p.body)
    if (!v) return(n1)
  }
  else if (is.null(power)) {
    power <- stats::uniroot(function(power) eval(p.body) - n1, c(1e-05, 0.99999))$root
    if (!v) return(power)
  }
  else if (is.null(n.ratio)) {
    n.ratio <- stats::uniroot(function(n.ratio) eval(p.body) - n1, c(2/n1, 1e+07))$root
    if (!v) return(n.ratio)
  }
  else if (is.null(alpha)) {
    alpha <- stats::uniroot(function(alpha) eval(p.body) - n1, c(1e-10, 1 - 1e-10))$root
    if (!v) return(alpha)
  }
  else stop("internal error")

  # Generate output text
  METHOD <- "Two sample comparison of proportions power calculation"
  n <- c(n1, n1 * n.ratio)
  p <- c(p1, p2)

  # Print output as a power.htest object
  structure(list(`n1, n2` = n, `p1, p2` = p, margin = margin,
                 alpha = alpha, power = power, sides = sides,
                 method = METHOD), class = "power.htest")
}

