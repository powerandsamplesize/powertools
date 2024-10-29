#' Power for individually randomized group treatment trial with continuous outcome
#'
#' @description
#' Computes power and sample size for an individually randomized group treatment trial with
#' a continuous outcome, in which after individual randomization, individuals in the
#' intervention/treatment arm are clustered. Can solve for power, J, m, n, delta or alpha.
#'
#' @details
#' Power is solved for using noncentral t or F distribution; other quantities
#' (for example, sample sizes) are solved
#' for using normal approximation.
#'
#'
#' @param m The number of subjects per cluster in the intervention arm.
#' @param J The total number of clusters in the intervention arm.
#' @param n The total number of participants in the control arm.
#' @param delta The difference between the intervention and control means under the alternative minus the difference under the null hypothesis.
#' @param sd The total standard deviation of the outcome variable in the control arm; defaults to 1.
#' @param icc The intraclass correlation coefficient in the treatment arm; defaults to 0.
#' @param Theta The ratio of the total variance in the intervention and control arms; defaults to 1.
#' @param alpha The significance level (type 1 error rate); defaults to 0.05.
#' @param power The specified level of power.
#' @param sides Either 1 or 2 (default) to specify a one- or two-sided hypothesis test.
#' @param tol The desired accuracy (convergence tolerance) for uniroot.
#' @param v Either TRUE for verbose output or FALSE (default) to output computed argument only.
#'
#' @return A list of the arguments (including the computed one).
#' @export
#'
#' @examples
#' irgtt.cont(m = 10, J = 12, n = 120, delta = 0.4, icc = 0.05, Theta = 1, power = NULL)
#' irgtt.cont(m = 10, J = 12, n = NULL, delta = 0.4, icc = 0.05, Theta = 1, power = 0.8)

irgtt.cont <- function (m = NULL, J = NULL, n = NULL, delta = NULL, sd = 1,
                        icc = 0, Theta = 1, alpha = 0.05, power = NULL, sides = 2,
                        tol = .Machine$double.eps^0.25, v = FALSE) {

  # Check if the arguments are specified correctly
  check.many(list(m, J, n, delta, alpha, power), "oneof")
  check.param(m, "pos")
  check.param(J, "min", min = 2)
  check.param(n, "pos")
  check.param(delta, "num")
  check.param(sd, "req"); check.param(sd, "pos")
  check.param(icc, "req"); check.param(icc, "uniti")
  check.param(Theta, "req"); check.param(Theta, "pos")
  check.param(alpha, "unit")
  check.param(power, "unit")
  check.param(sides, "req"); check.param(sides, "vals", valslist = c(1, 2))
  check.param(v, "req"); check.param(v, "bool")

  # Calculate power
  if (sides == 1)
    p.body <- quote({
      de <- 1 + (m - 1) * icc
      Uc <- sd^2 / n
      Ue <- Theta * sd^2 * de / (m * J)
      df <- (Ue^2 * (J + 1)/(J - 1) + Uc^2 * (n+1)/(n-1)) / (Ue^2 * (J+1)/(J-1)^2 + Uc^2 * (n+1)/(n-1)^2)
      d <- delta / sd
      lambda <- d / sqrt(Theta * de / (m * J) + 1 / n)
      crit <- stats::qt(1 - alpha, df = df)
      1 - stats::pt(crit, df = df, ncp = lambda)
    })
  else if (sides == 2)
    p.body <- quote({
      de <- 1 + (m - 1) * icc
      Uc <- sd^2 / n
      Ue <- Theta * sd^2 * de / (m * J)
      df2 <- (Ue^2 * (J + 1)/(J - 1) + Uc^2 * (n+1)/(n-1)) / (Ue^2 * (J+1)/(J-1)^2 + Uc^2 * (n+1)/(n-1)^2)
      d <- delta / sd
      lambda <- d / sqrt(Theta * de / (m * J) + 1 / n)
      crit <- stats::qf(1 - alpha, df1 = 1, df2 = df2)
      1 - stats::pf(crit, df1 = 1, df2 = df2, ncp = lambda^2)
    })

  # use this function (normal approx) if solving for anything other than power
  if (sides == 1)
    p.body2 <- quote({
      de <- 1 + (m - 1) * icc
      df <- 10000
      d <- delta / sd
      lambda <- d / sqrt(Theta * de / (m * J) + 1 / n)
      crit <- stats::qt(1 - alpha, df = df)
      1 - stats::pt(crit, df = df, ncp = lambda)
    })
  else if (sides == 2)
    p.body2 <- quote({
      de <- 1 + (m - 1) * icc
      df2 <- 10000
      d <- delta / sd
      lambda <- d / sqrt(Theta * de / (m * J) + 1 / n)
      crit <- stats::qf(1 - alpha, df1 = 1, df2 = df2)
      1 - stats::pf(crit, df1 = 1, df2 = df2, ncp = lambda^2)
    })

  # Use uniroot to calculate missing argument
  if (is.null(alpha)) {
    alpha <- stats::uniroot(function(alpha) eval(p.body2) - power, interval = c(1e-10, 1 - 1e-10), tol = tol)$root
    if (!v) return(alpha)
  }
  else if (is.null(power)) {
    power <- eval(p.body)
    if (!v) return(power)
  }
  else if (is.null(J)) {
    J <- stats::uniroot(function(J) eval(p.body2) - power, interval = c(2 + 1e-10, 1e+07), tol = tol)$root
    if (!v) return(J)
  }
  else if (is.null(m)) {
    m <- stats::uniroot(function(m) eval(p.body2) - power, interval = c(2 + 1e-10, 1e+07), tol = tol)$root
    if (!v) return(m)
  }
  else if (is.null(n)) {
    n <- stats::uniroot(function(n) eval(p.body2) - power, interval = c(2 + 1e-10, 1e+07), tol = tol)$root
    if (!v) return(n)
  }
  else if (is.null(delta)) {
    delta <- stats::uniroot(function(delta) eval(p.body2) - power, interval = c(1e-07, 1e+07), tol = tol)$root
    if (!v) return(delta)
  }

  mjn <- c(m, J, n)

  # Print output as a power.htest object
  METHOD <- "Power for individually randomized group treatment trial with continuous outcome"
  structure(list(`m, J, n` = mjn, delta = delta, sd = sd, icc = icc,
                 Theta = Theta, alpha = alpha, power = power, sides = sides,
                 method = METHOD), class = "power.htest")

}

