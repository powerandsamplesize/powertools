#' Power for test of treatment effect in longitudinal cluster randomized trial with baseline measurement
#'
#' @description
#' This function computes power and sample size for a cluster randomized trial in which a continuous
#' outcome variable is measured during both baseline and follow-up periods among
#' the cluster members, and it is planned that the outcome data will be analyzed using
#' a linear mixed model in which the dependent variable vector includes both baseline and follow up
#' measurements and there is a random intercept for cluster.
#' This function can solve for power, J1, J.ratio, m or delta.
#'
#'
#' @param m The number of subjects measured during each cluster-period.
#' @param J1 The number of clusters in arm 1.
#' @param J.ratio The ratio J2/J1 between the number of clusters in the two arms; defaults to 1 (equal clusters per arm).
#' @param delta The difference between the intervention and control means under the alternative minus the difference under the null hypothesis.
#' @param sd The total standard deviation of the outcome variable; defaults to 1.
#' @param icc The within-cluster, within-period intraclass correlation coefficient; defaults to 0.
#' @param cac The cluster autocorrelation; defaults to 0.
#' @param sac The subject autocorrelation; defaults to 0.
#' @param alpha The significance level (type 1 error rate); defaults to 0.05.
#' @param power The specified level of power.
#' @param sides Either 1 or 2 (default) to specify a one- or two- sided hypothesis test.
#' @param v Either TRUE for verbose output or FALSE (default) to output computed argument only.
#' @return A list of the arguments (including the computed one).
#' @export
#'
#' @details
#' The intraclass correlation coefficient (icc) is the correlation between
#' two observations from different subjects in the same cluster and same time period.
#' Denote the correlation between observations from two
#' different subjects in the same cluster but different time periods as iccb.
#' The cluster autocorrelation (cac) is iccb/icc and is interpreted as the proportion of
#' the cluster-level variance that is time-invariant.
#' Denote the correlation between two observations from the same subject
#' in different time periods as rhoa.
#' The subject autocorrelation (sac) is (rhoa - icc)/(iccb - icc) and is interpreted as
#' the proportion of the subject-level variance that is time-invariant.
#' The sac is only relevant for design in which the same subjects are measured at both baseline and follow up.
#' If different subjects are measured during different time periods, sac should be set to zero.
#'
#'
#' @examples
#' crt.long.cont(m = 30, J1 = 8, delta = 0.3, icc = 0.05, cac = 0.4, sac = 0.5)


crt.long.cont <- function (m = NULL, J1 = NULL, J.ratio = 1, delta = NULL, sd = 1,
                           icc = 0, cac = 0, sac = 0,
                           alpha = 0.05, power = NULL, sides = 2,
                           v = FALSE) {

  # Check if the arguments are specified correctly
  check.many(list(m, J1, J.ratio, delta, alpha, power), "oneof")
  check.param(m, "pos")
  check.param(J.ratio, "pos")
  check.param(J1, "min", min = 2)
  if (!is.null(J1) & !is.null(J.ratio))
    check.param(J1 * J.ratio, "min", min = 2)
  check.param(delta, "num")
  check.param(sd, "req"); check.param(sd, "pos")
  check.param(icc, "req"); check.param(icc, "uniti")
  check.param(cac, "req"); check.param(cac, "uniti")
  check.param(sac, "req"); check.param(sac, "uniti")
  check.param(alpha, "unit")
  check.param(power, "unit")
  check.param(sides, "req"); check.param(sides, "vals", valslist = c(1, 2))
  check.param(v, "req"); check.param(v, "bool")

  # Calculate power
  if (sides == 1)
    p.body <- quote({
      J <- J1 + J1 * J.ratio
      N <- m * J
      df <- J - 2
      d <- delta / sd
      de.pa <- 1 + (m - 1) * icc
      r <- m * icc * cac / de.pa + (1 - icc) * sac / de.pa
      de <- 4 * (1 - r^2) * de.pa
      ncp <- d / sqrt(de / N)
      crit <- stats::qt(1 - alpha, df)
      1 - stats::pt(crit, df, ncp)
    })
  else if (sides == 2)
    p.body <- quote({
      J <- J1 + J1 * J.ratio
      N <- m * J
      df2 <- J - 2
      d <- delta / sd
      de.pa <- 1 + (m - 1) * icc
      r <- m * icc * cac / de.pa + (1 - icc) * sac / de.pa
      de <- 4 * (1 - r^2) * de.pa
      ncp <- d / sqrt(de / N)
      crit <- stats::qf(1 - alpha, 1, df2)
      1 - stats::pf(crit, 1, df2, ncp^2)
    })

  # Use uniroot to calculate missing argument
  if (is.null(alpha)) {
    alpha <- stats::uniroot(function(alpha) eval(p.body) - power, c(1e-10, 1 - 1e-10))$root
    if (!v) return(alpha)
  }
  else if (is.null(power)) {
    power <- eval(p.body)
    if (!v) return(power)
  }
  else if (is.null(J1)) {
    J1 <- stats::uniroot(function(J1) eval(p.body) - power, c(2 + 1e-10, 1e+07))$root
    if (!v) return(J1)
  }
  else if (is.null(J.ratio)) {
    J.ratio <- stats::uniroot(function(J.ratio) eval(p.body) - power, c(1e-10, 1e+07))$root
    if (!v) return(J.ratio)
  }
  else if (is.null(m)) {
    m <- stats::uniroot(function(m) eval(p.body) - power, c(2 + 1e-10, 1e+07))$root
    if (!v) return(m)
  }
  else if (is.null(delta)) {
    delta <- stats::uniroot(function(delta) eval(p.body) - power, c(1e-07, 1e+07))$root
    if (!v) return(delta)
  }

  # Generate output text
  METHOD <- "Power for test of treatment effect in a cluster randomized trial with baseline measurement"
  J <- c(J1, J1 * J.ratio)
  de.pa <- 1 + (m - 1) * icc
  r <- m * icc * cac / de.pa + (1 - icc) * sac / de.pa
  ccs <- c(icc, cac, sac)
  out <- list(m = m, `J1, J2` = J, delta = delta, sd = sd, `icc, cac, sac` = ccs, r = r,
              alpha = alpha, power = power, sides = sides, method = METHOD)

  # Print output as a power.htest object
  structure(out, class = "power.htest")

}
