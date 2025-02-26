#' Power for individual randomized group treatment trial with binary outcome
#'
#' @description
#' Computes power and sample size for an individually randomized group treatment trial with
#' a binary outcome, in which after individual randomization, individuals in the
#' intervention/treatment arm are clustered. Can solve for power, J, m, n, or alpha.
#'
#'
#' @param m The number of subjects per cluster in the intervention arm.
#' @param J The total number of clusters in the intervention arm.
#' @param n The total number of participants in the control arm.
#' @param p1 The probability of the outcome in the control arm.
#' @param p2 The probability of the outcome in the intervention arm.
#' @param icc The intraclass correlation coefficient in the intervention arm; defaults to 0.
#' @param alpha The significance level (type 1 error rate); defaults to 0.05.
#' @param power The specified level of power.
#' @param sides Either 1 or 2 (default) to specify a one- or two-sided hypothesis test.
#' @param v Either TRUE for verbose output or FALSE (default) to output computed argument only.
#'
#' @return A list of the arguments (including the computed one).
#' @export
#'
#' @examples
#' irgtt.bin(m = 20, J = 5, n = 100, p1 = 0.8, p2 = 0.6, icc = 0.04, sides = 2)
#' irgtt.bin(m = 20, J = 6, n = 120, p1 = 0.8, p2 = 0.6, icc = 0.04, sides = 2)

irgtt.bin <- function (m = NULL, J = NULL, n = NULL, p1 = NULL, p2 = NULL,
                       icc = 0, alpha = 0.05, power = NULL, sides = 2,
                       v = FALSE) {

  # Check if the arguments are specified correctly
  check.many(list(m, J, n, alpha, power), "oneof")
  check.param(m, "pos")
  check.param(J, "min", min = 2)
  check.param(n, "pos")
  check.param(p1, "req"); check.param(p1, "unit")
  check.param(p2, "req"); check.param(p2, "unit")
  check.param(icc, "req"); check.param(icc, "uniti")
  check.param(alpha, "unit")
  check.param(power, "unit")
  check.param(sides, "req"); check.param(sides, "vals", valslist = c(1, 2))
  check.param(v, "req"); check.param(v, "bool")

  # Calculate power
  p.body <- quote({
    de <- 1 + (m - 1) * icc
    vard <- p2 * (1 - p2) * de / (m * J) + p1 * (1 - p1) / n
    ncp <- abs(p2 - p1) / sqrt(vard)
    crit <- stats::qnorm(alpha / sides, lower.tail = F)
    stats::pnorm(crit, ncp, lower.tail = F)
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
  else if (is.null(n)) {
    n <- safe.uniroot(function(n) eval(p.body) - power, c(4 + 1e-10, 1e+07))$root
    if (!v) return(n)
  }

  mjn <- c(m, J, n)
  p <- c(p1, p2)

  # Print output as a power.htest object
  METHOD <- "Power for an individual randomized group treatment trial with binary outcomes"
  structure(list(`m, J, n` = mjn, `p1, p2` = p, icc = icc,
                 alpha = alpha, power = power, sides = sides,
                 method = METHOD), class = "power.htest")

}

