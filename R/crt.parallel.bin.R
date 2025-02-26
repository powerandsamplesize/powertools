#' Power calculation for cluster randomized trial with binary outcome
#'
#' @description
#' This function performs power and sample size calculations for a two-arm cluster randomized trial
#' with a binary outcome. It assumes the outcome analysis will be conducted using a mixed effect logistic
#' regression model that has a random intercept for cluster. Equal allocation of clusters to arms
#' is assumed. Can solve for power, J, m or alpha.
#'
#' @details
#' For help selecting a reasonable value for sigma.u, consider using the crt.varexplore function.
#'
#'
#' @param m The number of subjects per cluster.
#' @param m.sd The standard deviation of cluster sizes (provide if unequal number of participants per cluster); defaults to 0.
#' @param J The total number of clusters (over both arms).
#' @param pc The probability of the outcome in control clusters.
#' @param pt The probability of the outcome in treatment clusters.
#' @param sigma.u Standard deviation of the cluster random effect (random intercept).
#' @param alpha The significance level (type 1 error rate); defaults to 0.05.
#' @param power The specified level of power.
#' @param sides Either 1 or 2 (default) to specify a one- or two- sided hypothesis test.
#' @param v Either TRUE for verbose output or FALSE (default) to output computed argument only.
#'
#' @return A list of the arguments (including the computed one).
#' @export
#'
#' @examples
#' crt.parallel.bin(m = 60, J = NULL, pc = 0.25, pt = 0.15, sigma.u = 0.3, power = 0.8)
#' crt.parallel.bin(m = 60, m.sd = 1, J = NULL, pc = 0.25, pt = 0.15, sigma.u = 0.3, power = 0.8)

crt.parallel.bin <- function (m = NULL, m.sd = 0, J = NULL,
                              pc = NULL, pt = NULL, sigma.u = NULL,
                              alpha = 0.05, power = NULL, sides = 2,
                              v = FALSE) {

  # Check if the arguments are specified correctly
  check.many(list(m, J, alpha, power), "oneof")
  check.param(m, "pos")
  check.param(J, "pos")
  check.param(alpha, "unit")
  check.param(power, "unit")
  check.param(m.sd, "req"); check.param(m.sd, "min", min = 0)
  check.param(pc, "req"); check.param(pc, "unit")
  check.param(pt, "req"); check.param(pt, "unit")
  check.param(sigma.u, "req"); check.param(sigma.u, "pos")
  check.param(sides, "req"); check.param(sides, "vals", valslist = c(1, 2))
  check.param(v, "req"); check.param(v, "bool")

  # Calculate power
  p.body <- quote({ # Where does RE go?
    or <- (pt / (1 - pt)) / (pc / (1 - pc))
    gammaA <- abs(log(or))
    ssq.e <- (1 / 2) * (1 / (pc * (1 - pc)) + 1 / (pt * (1 - pt)))

    RE <- re.clustsize.bin(m = m, m.sd = m.sd, pc = pc, pt = pt, sigma.u = sigma.u)

    var <- 1.2 * (4 * (ssq.e  + m * sigma.u^2)) / (m *J) / RE
    za <- stats::qnorm(1 - alpha / sides)
    stats::pnorm(gammaA / sqrt(var) - za)
  })

  if (is.null(alpha)) {
    alpha <- safe.uniroot(function(alpha) eval(p.body)  - power, c(1e-10, 1 - 1e-10))$root
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

  # Generate output text
  METHOD <-"Power for a cluster randomized trial with a binary outcome"
  m <- ifelse(m.sd == 0, m, paste0(m, " (sd = ", m.sd, ")"))
  p <- c(pc, pt)

  # Print output as a power.htest object depending on which inputs were given
  structure(list(m = m, J = J, `pc, pt` = p, sigma.u = sigma.u,
                 alpha = alpha, power = power,
                 method = METHOD), class = "power.htest")

}
