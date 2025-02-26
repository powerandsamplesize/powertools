#' Power calculation for relative risk
#'
#' @description
#' Performs power and sample size calculations for a test of relative risk, p2/p1.
#' Can solve for power, N1, n.ratio or alpha.
#'
#'
#' @param n1 The sample size for group 1.
#' @param n.ratio The ratio n2/n1 between the sample sizes of two groups; defaults to 1 (equal group sizes).
#' @param p1 The proportion in group 1.
#' @param p2 The proportion in group 2.
#' @param RR0 The relative risk under the null (p2/p1); defaults to 1.
#' @param alpha The significance level (type 1 error rate); defaults to 0.05.
#' @param power The specified level of power.
#' @param sides Either 1 or 2 (default) to specify a one- or two- sided hypothesis test.
#' @param v Either TRUE for verbose output or FALSE (default) to output computed argument only.
#'
#' @return A list of the arguments (including the computed one).
#' @export
#'
#' @examples
#' relrisk(n1 = NULL, n.ratio = 1/6, p1 = 0.1, p2 = 0.1 * 2, power = 0.8)

relrisk <- function (n1 = NULL, n.ratio = 1, p1 = NULL, p2 = NULL, RR0 = 1,
                     alpha = 0.05, power = NULL, sides = 2, v = FALSE) {

  # Check if the arguments are specified correctly
  check.many(list(n1, n.ratio, alpha, power), "oneof")
  check.param(n1, "pos")
  check.param(n.ratio, "pos")
  check.param(p1, "req"); check.param(p1, "unit")
  check.param(p2, "req"); check.param(p2, "unit")
  check.param(RR0, "req"); check.param(RR0, "pos")
  check.param(alpha, "unit")
  check.param(power, "unit")
  check.param(sides, "req"); check.param(sides, "vals", valslist = c(1, 2))
  check.param(v, "req"); check.param(v, "bool")

  # Calculate test statistic
  p.body <- quote({
    RR <- p2 / p1
    d <- abs(log(RR) - log(RR0))
    q1 <- 1 - p1
    q2 <- 1 - p2
    denom <- (1 / n1) * (q1 / p1 + q2 / (n.ratio * p2))
    (stats:: pnorm(qnorm(alpha / sides) + d / sqrt(denom)))
  })

  # Use safe.uniroot function to calculate missing argument
  if (is.null(power)) {
    power <- eval(p.body)
    if (!v) return(power)
  }
  else if (is.null(n1)) {
    n1 <- safe.uniroot(function(n1) eval(p.body) - power, c(2, 1e+07))$root
    if (!v) return(n1)
  }
  else if (is.null(n.ratio)) {
    n.ratio <- safe.uniroot(function(n.ratio) eval(p.body) - power, c(2/n1, 1e+07))$root
    if (!v) return(n.ratio)
  }
  else if (is.null(alpha)) {
    alpha <- safe.uniroot(function(alpha) eval(p.body) - power, c(1e-10, 1 - 1e-10))$root
    if (!v) return(alpha)
  }
  else stop("internal error")

  # Generate output text
  METHOD <- "Relative risk power calculation"
  n <- c(n1, n1 * n.ratio)
  p <- c(p1, p2)
  RR <- c(p2 / p1, RR0)

  # Print output as a power.htest object
  structure(list(`n1, n2` = n, `p1, p2` = p, `RR, RR0` = RR,
                 alpha = alpha, power = power, sides = sides,
                 method = METHOD), class = "power.htest")
}

