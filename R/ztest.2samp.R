#' Power calculation for two-sample z test
#'
#' @description
#' This function performs power and sample size calculations for a two-sample z test. which is
#' analogous to a two-sample t test but with the variances assumed to be known.
#' The function allows for unequal sample sizes and/or variances in the two groups.
#' This function is provided largely for pedagogical purposes; in general,
#' for real studies, the two-sample t test procedure should be used.
#'
#' @details
#' For a noninferiority or superiority by a margin test, the margin should be subtracted
#' as part of the specification of delta and a one-sided test should be used. See Crespi (2025)
#' for specific guidance.
#'
#' @param n1 The sample size for group 1.
#' @param n.ratio The ratio n2/n1 between the sample sizes of two groups; defaults to 1 (equal group sizes).
#' @param delta DeltaA (true difference mu1 - mu2) - Delta0 (difference under the null). For a noninferiority
#' or superiority by a margin test,
#' the margin should be subtracted, that is, delta = DeltaA - Delta0 - margin.
#' @param sd1 The standard deviation for group 1; defaults to 1 (equal standard deviations in the two groups).
#' @param sd.ratio The ratio sd2/sd1 between the standard deviations of the two groups.
#' @param alpha The significance level (type 1 error rate); defaults to 0.05.
#' @param power The specified level of power.
#' @param sides Either 1 or 2 (default) to specify a one- or two- sided hypothesis test.
#' @param v Either TRUE for verbose output or FALSE (default) to output computed argument only.
#'
#' @return A list of the arguments (including the computed one).
#' @export
#'
#' @examples
#' ztest.2samp(n1 = NULL, n.ratio = 1, delta = 0.5, sd1 = 1, power = 0.8, sides = 2)

ztest.2samp <- function (n1 = NULL, n.ratio = 1, delta = NULL,
                         sd1 = 1, sd.ratio = 1,
                         alpha = 0.05, power = NULL, sides = 2,
                         v = FALSE) {

  # Check if the arguments are specified correctly
  check.many(list(n1, n.ratio, delta, sd1, sd.ratio, alpha, power), "oneof")
  check.param(n1, "pos")
  check.param(n.ratio, "pos")
  check.param(sd1, "pos")
  check.param(sd.ratio, "pos")
  check.param(delta, "num")
  check.param(alpha, "unit")
  check.param(power, "unit")
  check.param(sides, "req"); check.param(sides, "vals", valslist = c(1, 2))
  check.param(v, "req"); check.param(v, "bool")

  # Calculate power
  if (sides == 1)
    p.body <- quote({
      d <- abs(delta)
      stats::pnorm(stats::qnorm(alpha) +
                   d / sqrt((sd1 * sd.ratio)^2 / (n1 * n.ratio) + sd1^2 / n1))})
  else if (sides == 2)
    p.body <- quote({
      stats::pnorm(stats::qnorm(alpha / 2) +
                   delta / sqrt((sd1 * sd.ratio)^2 / (n1 * n.ratio) + sd1^2 / n1)) +
      stats::pnorm(stats::qnorm(alpha / 2) -
                   delta / sqrt((sd1 * sd.ratio)^2 / (n1 * n.ratio) + sd1^2 / n1))
    })

  # Use stats::uniroot function to calculate missing argument
  if (is.null(power)) {
    power <- eval(p.body)
    if (!v) return(power)
  }
  else if (is.null(n1)) {
    n1 <- stats::uniroot(function(n1) eval(p.body) - power, c(2, 1e+07))$root
    if (!v) return(n1)
  }
  else if (is.null(n.ratio)) {
    n.ratio <- stats::uniroot(function(n.ratio) eval(p.body) - power,c(2/n1, 1e+07))$root
    if (!v) return(n.ratio)
  }
  else if (is.null(sd1)) {
    sd1 <- stats::uniroot(function(sd1) eval(p.body) - power, delta * c(1e-07, 1e+07))$root
    if (!v) return(sd1)
  }
  else if (is.null(sd.ratio)) {
    sd.ratio <- stats::uniroot(function(sd.ratio) eval(p.body) - power, c(1e-07, 1e+07))$root
    if (!v) return(sd.ratio)
  }
  else if (is.null(delta)) {
    delta <- stats::uniroot(function(delta) eval(p.body) - power,  c(1e-07, 1e+07))$root
    if (!v) return(delta)
  }
  else if (is.null(alpha)) {
    alpha <- stats::uniroot(function(alpha) eval(p.body) - power, c(1e-10, 1 - 1e-10))$root
    if (!v) return(alpha)
  }
  else stop("internal error")

  # Generate output text
  METHOD <- "Two-sample z test power calculation"
  n <- c(n1, n1 * n.ratio)
  sd <- c(sd1, sd1 * sd.ratio)

  # Print output as a power.htest object
  structure(list(`n1, n2` = n, delta = delta, `sd1, sd2` = sd, alpha = alpha,
                 power = power, sides = sides,
                 method = METHOD), class = "power.htest")
}
