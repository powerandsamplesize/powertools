#' Power for test of heterogeneity of treatment effect in a multisite trial
#'
#' @description
#' Calculates power for a test of heterogeneity of the treatment effect
#' in a multisite trial with a continuous outcome variable. In particular,
#' let sigma.u1^2 be the variance of the site-level treatment effects.
#' The test of heterogeneity tests the null hypothesis that sigma.u1^2 = 0
#' versus the alternative that sigma.u1^2 > 0.
#' Can solve for power, J, or m.
#'
#' @details
#' See Crespi CM (2025) for details.
#'
#'

#' @param m The mean number of subjects per site.
#' @param alloc.ratio The allocation ratio of intervention/control subjects per site; defaults to 1.
#' @param J The total number of sites.
#' @param VR The variance ratio (variance of the site-level treatment effects divided by variance of
#' observations within sites) under the alternative.
#' @param alpha The significance level (type 1 error rate); defaults to 0.05.
#' @param power The specified level of power.
#' @param v Either TRUE for verbose output or FALSE (default) to output computed argument only.
#'
#' @return A list of the arguments (including the computed power).
#' @export
#'
#' @examples
#' multisite.hte(m = 10, J = 30, VR = 8 / 36)

multisite.hte <- function (m = NULL, alloc.ratio = 1, J = NULL, VR = NULL,
                           alpha = 0.05, power = NULL, v = FALSE) {

  # Check if the arguments are specified correctly
  check.many(list(m, J, alpha, power), "oneof")
  check.param(m, "pos")
  check.param(alloc.ratio, "req"); check.param(alloc.ratio, "pos")
  check.param(J, "min", min = 2)
  check.param(VR, "req"); check.param(VR, "pos")
  check.param(alpha, "unit")
  check.param(power, "unit")
  check.param(v, "req"); check.param(v, "bool")

  p.body <- quote({
    omega <- 1 + m * VR / 4
    df1 <- J - 1
    df2 <- J * (m - 2)
    crit <- stats::qf(1 - alpha, df1, df2)
    1 - stats::pf(crit / omega, df1, df2)
  })

  # Use uniroot to calculate missing argument
  if (is.null(alpha)) {
    alpha <- safe.uniroot(function(alpha) eval(p.body) - power, c(1e-10, 1 - 1e-10))$root
    if (!v) return(alpha)
  }
  else if (is.null(power)) {
    power <- eval(p.body)
    if (!v) return(power)
  }
  else if (is.null(J)) {
    J <- safe.uniroot(function(J) eval(p.body) - power, c(2 + 1e-10, 1e+07))$root
    if (!v) return(J)
  }
  else if (is.null(m)) {
    m <- safe.uniroot(function(m) eval(p.body) - power, c(2 + 1e-10, 1e+07))$root
    if (!v) return(m)
  }
  else stop("internal error")

  # Generate output text
  METHOD <-"Power for test of heterogeneity of treatment effect in a multisite trial"
  NOTE <- "m1, m2 are the number of subjects within site in condition 1, condition 2\n      (total of m1 + m2 per site)"
  c <- m / (alloc.ratio + 1)
  t <- alloc.ratio * c
  m <- paste0(t, ", ", c)

  # Print output as a power.htest object
  structure(list(`m1, m2` = m, J = J, VR = VR,
                 alpha = alpha, power = power,
                 method = METHOD, note = NOTE), class = "power.htest")

}
