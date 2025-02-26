#' Power calculation for rank-sum test
#'
#' @description
#' This function performs power and sample size calculations for the
#' Wilcoxon-Mann-Whitney rank-sum test, also called the Mann-Whitney U test,
#' which is the nonparametric analog of the two independent sample t test.
#' Calculations are based on the approximation given by
#' Noether (1987) Sample size determination for some common nonparametric tests. JASA 82(398):645-647.
#'
#' @details
#' Due to symmetry, the power for p is equal to the power for 1 - p. Therefore,
#' when solving for p, two values, p and 1 - p, are returned.
#'
#' @param n1 The sample size in group 1.
#' @param n.ratio The ratio n2/n1 between the sample sizes of two groups; defaults to 1 (equal group sizes).
#' @param p The probability that an observation in group 2 is greater than an observation in group 1 (P(Y>X).
#' @param alpha The significance level (type 1 error rate); defaults to 0.05.
#' @param power The specified level of power.
#' @param sides Either 1 or 2 (default) to specify a one- or two-sided hypothesis test.
#' @param v Either TRUE for verbose output or FALSE (default) to output computed argument only.
#'
#' @return A list of the arguments (including the computed one).
#' @export
#'
#' @examples
#' ranksum(n1 = 10, n.ratio = 1, p = 0.8, alpha = 0.05, power = NULL, sides = 2)


ranksum <- function (n1 = NULL, n.ratio = 1, p = NULL, alpha = 0.05,
                     power = NULL, sides = 2, v = FALSE) {

  # Check if the arguments are specified correctly
  check.many(list(n1, n.ratio, alpha, power, p), "oneof")
  check.param(n1, "pos"); check.param(n1, "min", min = 2)
  check.param(n.ratio, "pos")
  check.param(alpha, "unit")
  check.param(power, "unit")
  check.param(p, "unit")
  check.param(sides, "req"); check.param(sides, "vals", valslist = c(1, 2))
  check.param(v, "req"); check.param(v, "bool")

  # power equation
  p.body <- quote({
    stats::pnorm(stats::qnorm(alpha / sides) + sqrt(12 * n1 * n.ratio / (1 + n.ratio)) * abs(p - 0.5))
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
  else if (is.null(p)) {
    p <- safe.uniroot(function(p) eval(p.body) - power, c(0.5, 1 - 1e-10))$root
    p <- c(p, 1 - p)
    if (!v) return(p)
  }
  else if (is.null(alpha)) {
    alpha <- safe.uniroot(function(alpha) eval(p.body) - power, c(1e-10, 1 - 1e-10))$root
    if (!v) return(alpha)
  }
  else stop("internal error")

  # Generate output text
  METHOD <- "Rank-sum test power calculation"
  n <- c(n1, n1 * n.ratio)

  # Print output as a power.htest object
  structure(list(n = n, p = p, alpha = alpha,
                 power = power, sides = sides,
                 method = METHOD), class = "power.htest")
}

