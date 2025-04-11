#' Power calculation for two-sample t test
#'
#' @description
#' This function computes power and sample size for a two-sample t test.
#' Unequal sample sizes and/or unequal variances in the two groups are allowed.
#' Can solve for power, n1, n.ratio, sd1, sd.ratio, delta or alpha.
#'
#' @details
#' For a noninferiority or superiority by a margin test, the margin should be subtracted
#' as part of the specification of delta and a one-sided test should be used. See Crespi (2025)
#' for specific guidance. For an equivalence test for two means, see the sampleN.TOST function from the
#' PowerTOST package.
#'
#'
#' @param n1 The sample size for group 1.
#' @param n.ratio The ratio n2/n1 between the sample sizes of two groups; defaults to 1 (equal group sizes).
#' @param delta DeltaA (true difference mu1 - mu2) - Delta0 (difference under the null). For a noninferiority
#' or superiority by a margin test,
#' the margin should be subtracted, that is, delta = DeltaA - Delta0 - margin.
#' @param sd1 The standard deviation for group 1; defaults to 1 (equal standard deviations in the two groups).
#' @param sd.ratio The ratio sd2/sd1 between the standard deviations of the two groups.
#' @param df.method Method for calculating the degrees of freedom: "welch" (default) or "classical", which is n1+n2-2.
#' @param alpha The significance level (type 1 error rate); defaults to 0.05.
#' @param power The specified level of power.
#' @param sides Either 1 or 2 (default) to specify a one- or two- sided hypothesis test.
#' @param v Either TRUE for verbose output or FALSE (default) to output computed argument only.
#'
#' @return A list of the arguments (including the computed one).
#' @export ttest.2samp
#'
#' @examples
#' ttest.2samp(n1 = 50, delta = 2, sd1 = 5, sides = 1)
#' ttest.2samp(n1 = NULL, n.ratio = 2, delta = 0.5, sd1 = 1, power = 0.8, sides = 2)
#' ttest.2samp(n1 = 49, n.ratio = 2, delta = 0.5, sd1 = 1, power = NULL, sides = 2)
#' ttest.2samp(n1 = 25, n.ratio = 3, delta = 3, sd1 = 4, sd.ratio = 1.5, alpha = 0.025, sides = 1)
#' ttest.2samp(n1 = NULL, delta = 0.5, sd1 = 1, power = 0.8, sides = 2)

ttest.2samp <- function (n1 = NULL, n.ratio = 1, delta = NULL,
                         sd1 = 1, sd.ratio = 1,
                         df.method = c("welch", "classical"),
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

  # Assign df method
  df.method <- match.arg(df.method)

  # Calculate df and ncp
  if (sides == 1)
    p.body <- quote({
      d <- abs(delta)
      df <- switch(df.method,
                   welch = (sd1^2 / n1 + (sd1 * sd.ratio)^2 / (n1 * n.ratio))^2 /
                   ((sd1^2 / n1)^2 / (n1 - 1) +
                   ((sd1 * sd.ratio)^2 / (n.ratio * n1))^2 / (n1 * n.ratio - 1)),
                   classical = (1 + n.ratio) * n1 - 2)
      stats::pt(stats::qt(alpha, df, lower.tail = FALSE), df,
                sqrt(n1 / (1 + sd.ratio^2 / n.ratio)) * d / sd1, lower.tail = FALSE)
    })
  else if (sides == 2)
    p.body <- quote({
      d <- abs(delta)
      df2 <- switch(df.method,
                    welch = (sd1^2 / n1 + (sd1 * sd.ratio)^2 / (n1 * n.ratio))^2 /
                    ((sd1^2 / n1)^2 / (n1 - 1) +
                    ((sd1 * sd.ratio)^2 / (n.ratio * n1))^2 / (n1 * n.ratio - 1)),
                    classical = (1 + n.ratio) * n1 - 2)
      q <- stats::qf(alpha, 1, df2, lower.tail = FALSE)
      ncp <- sqrt(n1 / (1 + sd.ratio^2 / n.ratio)) * d / sd1
      stats::pf(q, 1, df2, ncp^2, lower.tail = FALSE)
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
  else if (is.null(sd1)) {
    sd1 <- safe.uniroot(function(sd1) eval(p.body) - power, delta * c(1e-07, 1e+07))$root
    if (!v) return(sd1)
  }
  else if (is.null(sd.ratio)) {
    sd.ratio <- safe.uniroot(function(sd.ratio) eval(p.body) - power, c(1e-07, 1e+07))$root
    if (!v) return(sd.ratio)
  }
  else if (is.null(delta)) {
    delta <- safe.uniroot(function(delta) eval(p.body) - power,  c(1e-07, 1e+07))$root
    if (!v) return(delta)
  }
  else if (is.null(alpha)) {
    alpha <- safe.uniroot(function(alpha) eval(p.body) - power, c(1e-10, 1 - 1e-10))$root
    if (!alpha) return(alpha)
  }
  else stop("internal error")

  # Generate output text
  METHOD <- "Two-sample t test power calculation"
  n <- c(n1, n1 * n.ratio)
  sd <- c(sd1, sd1 * sd.ratio)

  # Print output as a power.htest object
  structure(list(`n1, n2` = n, delta = delta, `sd1, sd2` = sd, alpha = alpha,
                 power = power, sides = sides,
                 method = METHOD), class = "power.htest")
}
