#' Power for test of treatment effect in 2x2 crossover cluster randomized trial
#'
#' @description
#' Power and sample size calculation for a 2x2 crossover cluster randomized trial.
#' Can solve for power, number of clusters per arm (assumes equal number of
#' cluster per arm), m, delta or alpha.
#'
#'
#' @param m The number of subjects measured during each cluster-period.
#' @param J.arm The number of clusters in each arm.
#' @param delta The difference between the intervention and control means under the alternative minus the difference under the null hypothesis.
#' @param sd The total standard deviation of the outcome variable; defaults to 1.
#' @param icc The within-cluster, within-period intraclass correlation coefficient; defaults to 0.
#' @param icca The  within-cluster, within-subject correlation (correlation between two measurements within the same subject); defaults to 0.
#' @param iccb The within-cluster, between-period intraclass correlation coefficient. Either iccb OR cac must be specified.
#' @param cac The cluster autocorrelation. Either iccb OR cac must be specified.
#' @param sac The subject autocorrelation; defaults to 0.
#' @param alpha The significance level or type 1 error rate; defaults to 0.05.
#' @param power The specified level of power.
#' @param sides Either 1 or 2 (default) to specify a one- or two- sided hypothesis test.
#' @param v Either TRUE for verbose output or FALSE to output computed argument only.
#'
#' @return A list of the arguments (including the computed one).
#' @export
#'
#' @examples
#' crt.xo.cont(m = 30, J.arm = 4, delta = 0.3, icc = 0.05, cac = 0.8, sac = 0.4)
#' crt.xo.cont(m = 30, J.arm = 4, delta = 0.3, icc = 0.05, icca = 0.42, iccb = 0.04)
#' crt.xo.cont(m = 30, J.arm = 4, delta = 0.3, icc = 0.05, cac = 0.5)

crt.xo.cont <- function (m = NULL, J.arm = NULL, delta = NULL, sd = 1,
                         icc = 0, icca = 0, iccb = NULL, cac = NULL, sac = 0,
                         alpha = 0.05, power = NULL, sides = 2,
                         v = FALSE) {

  # Check if the arguments are specified correctly
  check.many(list(iccb, cac), "oneof")
  check.many(list(m, J.arm, delta, alpha, power), "oneof")
  check.param(m, "pos")
  check.param(J.arm, "min", min = 2)
  check.param(delta, "num")
  check.param(sd, "req"); check.param(sd, "pos")
  check.param(icc, "req"); check.param(icc, "uniti")
  check.param(icca, "req"); check.param(icca, "uniti")
  check.param(iccb, "uniti")
  check.param(cac, "uniti")
  check.param(sac, "req"); check.param(sac, "uniti")
  check.param(alpha, "unit")
  check.param(power, "unit")
  check.param(sides, "req"); check.param(sides, "vals", valslist = c(1, 2))
  check.param(v, "req"); check.param(v, "bool")

  if (is.null(iccb))
    iccb <- icc * cac
  if (is.null(cac))
    cac <- iccb / icc
  if (!is.null(cac))
    icca <- ifelse(sac == 0, 0, (1 - icc) * sac + icc * cac)

  # Calculate power
  if (sides == 1)
    p.body <- quote({
      J <- 2 * J.arm
      N <- m * J
      df <- 2 * J - 3
      d <- delta / sd

      de <- ifelse(sac == 0, 2 * (1 + (m - 1) * icc - m * iccb) , 2 * (1 - icca + (m - 1) * (icc - iccb)))
      ncp <- d / sqrt(de / N)
      crit <- stats::qt(1 - alpha, df)
      1 - stats::pt(crit, df, ncp)
    })
  else if (sides == 2)
    p.body <- quote({
      J <- 2 * J.arm
      N <- m * J
      df2 <- 2 * J - 3
      d <- delta / sd

      de <- ifelse(sac == 0, 2 * (1 + (m - 1) * icc - m * iccb) , 2 * (1 - icca + (m - 1) * (icc - iccb)))
      ncp <- d / sqrt(de / N)
      crit <- stats::qf(1 - alpha, 1, df2)
      1 - stats::pf(crit, 1, df2, ncp^2)
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
  else if (is.null(J.arm)) {
    J.arm <- safe.uniroot(function(J.arm) eval(p.body) - power, c(2 + 1e-10, 1e+07))$root
    if (!v) return(J.arm)
  }
  else if (is.null(m)) {
    m <- safe.uniroot(function(m) eval(p.body) - power, c(2 + 1e-10, 1e+07))$root
    if (!v) return(m)
  }
  else if (is.null(delta)) {
    delta <- safe.uniroot(function(delta) eval(p.body) - power, c(1e-07, 1e+07))$root
    if (!v) return(delta)
  }

  # Generate output text
  METHOD <- "Power for test of treatment effect in a 2x2 cluster randomized crossover trial"
  m <- m
  J <- c(J.arm, J.arm)
  iccs <- c(icc, icca, iccb)
  acs <- c(cac, sac) # outputs only SAC when CAC is NULL...how to generate CAC


  # Print output as a power.htest object
  out <- list(m = m, `J.arm1, J.arm2` = J, delta = delta, sd = sd,
              `icc, icca, iccb` = iccs, `CAC, SAC` = acs,
              alpha = alpha, power = power, sides = sides, method = METHOD)
  structure(out, class = "power.htest")

}



