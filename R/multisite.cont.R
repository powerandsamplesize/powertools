#' Power for test of average treatment effect in a multisite trial
#'
#' @description
#' Performs power and sample size calculations for a multisite trial with a
#' continuous (normal) outcome variable. Can solve for power, J, m or alpha.
#'
#' @details
#' In a multisite trial design, participants are randomized to conditions
#' within site.
#'
#' @param m The mean cluster/site size (number of participants per site).
#' @param m.sd The standard deviation of cluster/site sizes (provide if unequal number of participants per site);
#' defaults to 0.
#' @param alloc.ratio The allocation ratio of condition 2/condition 1 within site; defaults to 1.
#' @param J The total number of sites.
#' @param delta The difference between the condition 1 and condition 2 means under the alternative
#'  minus the difference under the null hypothesis.
#' @param sd The total standard deviation of the outcome variable; defaults to 1.
#' @param icc0 The proportion of total variance of the outcome attributable to variation in site-level means.
#' @param icc1 The proportion of total variance of the outcome attributable to variation in the treatment effect across sites.
#' @param Rsq The estimated R^2 for regressing the outcome on the covariates; defaults to 0.
#' @param alpha The significance level (type 1 error rate); defaults to 0.05.
#' @param power The specified level of power.
#' @param sides Either 1 or 2 (default) to specify a one- or two- sided hypothesis test.
#' @param v Either TRUE for verbose output or FALSE (default) to output computed argument only.
#'
#' @return A list of the arguments (including the computed one).
#' @export
#'
#' @examples
#' multisite.cont(m = 20, J = 10, delta = 3, sd = sqrt(40), icc0 = 0.1, icc1 = 0)
#' multisite.cont(m = 20, J = 10, delta = 3, sd = sqrt(48), icc0 = 0.095, icc1 = 0.048)
#' multisite.cont(m = 20, alloc.ratio = 1.5, J = 10, delta = 0.43, icc0 = 0.095, icc1 = 0.048)
#' multisite.cont(m = 10, J = NULL, delta = 0.5, sd = 1, icc0 = 0, icc1 = 0.05, power = 0.8)
#' multisite.cont(m = 20, m.sd = 5, J = 10, delta = 3, sd = sqrt(48), icc0 = 0.095, icc1 = 0.048)
#' multisite.cont(m = 20, J = 10, delta = 3, sd = sqrt(48), icc0 = 0.095,
#' icc1 = 0.048, Rsq = 0.5^2)

multisite.cont <- function (m = NULL, m.sd = 0, alloc.ratio = 1, J = NULL,
                            delta = NULL, sd = 1,
                            icc0 = NULL, icc1 = NULL, Rsq = 0,
                            alpha = 0.05, power = NULL, sides = 2,
                            v = FALSE) {

  # Check if the arguments are specified correctly
  check.many(list(m, J, delta, alpha, power), "oneof")
  check.param(m, "pos")
  check.param(m.sd, "req"); check.param(m.sd, "min", min = 0)
  check.param(alloc.ratio, "req"); check.param(alloc.ratio, "pos")
  check.param(J, "min", min = 2)
  check.param(delta, "num")
  check.param(sd, "req"); check.param(sd, "pos")
  check.param(icc0, "req"); check.param(icc0, "uniti")
  check.param(icc1, "req"); check.param(icc1, "uniti")
  check.param(Rsq, "req"); check.param(Rsq, "uniti")
  check.param(alpha, "unit")
  check.param(power, "unit")
  check.param(sides, "req"); check.param(sides, "vals", valslist = c(1, 2))
  check.param(v, "req"); check.param(v, "bool")

  # Calculate power
  if (sides == 1)
    p.body <- quote({
      N <- m * J
      df <- J - 1
      d <- delta / (sd * sqrt((1 - Rsq)))

      RE <- re.clustsize.cont(m = m, m.sd = m.sd, icc = icc1)

      c <- (1 + alloc.ratio)^2 / alloc.ratio
      ncp <- d / sqrt(c * (1 - icc0 + (4 * m / c - 1) * icc1) / N / RE)
      crit <- stats::qt(1 - alpha, df)
      1 - stats::pt(crit, df, ncp)
    })
  else if (sides == 2)
    p.body <- quote({
      N <- m * J
      df2 <- J - 1
      d <- delta / (sd * sqrt((1 - Rsq)))

      RE <- re.clustsize.cont(m = m, m.sd = m.sd, icc = icc1)

      c <- (1 + alloc.ratio)^2 / alloc.ratio
      ncp <- d / sqrt(c * (1 - icc0 + (4 * m / c - 1) * icc1) / N / RE)
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
  else if (is.null(J)) {
    J <- safe.uniroot(function(J) eval(p.body) - power, c(2 + 1e-10, 1e+07))$root
    if (!v) return(J)
  }
  else if (is.null(m)) {
    m <- safe.uniroot(function(m) eval(p.body) - power, c(2 + 1e-10, 1e+07))$root
    if (!v) return(m)
  }
  # else if (is.null(alloc.ratio)) {
  #   alloc.ratio <- safe.uniroot(function(alloc.ratio) eval(p.body) - power + 1e-5, c(1 + 1e-10, 1e+07))$root
  #   alloc.ratio <- c(alloc.ratio, 1 / alloc.ratio)
  #   if (!v) return(alloc.ratio)
  # }
  else if (is.null(delta)) {
    delta <- safe.uniroot(function(delta) eval(p.body) - power, c(1e-07, 1e+07))$root
    if (!v) return(delta)
  }
  else stop("internal error")

  # Generate output text
  METHOD <- "Power for test of average treatment effect in multisite trials"
  NOTE <- "m1, m2 are the number of subjects within site in condition 1, condition 2\n      (total of m1 + m2 per site). m1, m2 OR m2, m1 produce equivalent results."
  icc <- c(icc0, icc1)
  c <- round(m / (alloc.ratio + 1), 3)
  t <- round(alloc.ratio * c, 3)
  m <- ifelse(m.sd == 0, paste0(t, ", ", c), paste0(t, ", ",  c, " (sd = ", m.sd, ")"))

  out <- list(`m1, m2` = m, J = J, delta = delta, sd = sd,
              `icc0, icc1` = icc, Rsq = Rsq,
              alpha = alpha, power = power, sides = sides,
              method = METHOD, note = NOTE)

  # Print output as a power.htest object
  if (Rsq < 0.00000000001) out <- out[!names(out) %in% c("Rsq")]
  structure(out, class = "power.htest")

}

