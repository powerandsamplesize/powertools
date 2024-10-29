#' Power for a multisite trial with a binary outcome
#'
#' @description
#' Performs power and sample size calculations for a multisite trial with a
#' binary outcome. Can solve for power, J, m or alpha.
#'
#' @details
#' In a multisite trial design, participants are randomized to conditions
#' within site. Consider using ms.varexplore to select plausible values
#' for sigma.u.
#'
#'
#' @param m The total number of subjects in condition 1 + condition 2.
#' @param alloc.ratio The allocation ratio of condition 1/condition 2 within site; defaults to 1.
#' @param J The total number of sites.
#' @param pc The probability of the outcome in the control condition.
#' @param pt The probability of the outcome in the treatment condition.
#' @param sigma.u Standard deviation of the treatment effect across sites.
#' @param alpha The significance level or type 1 error rate; defaults to 0.05.
#' @param power The specified level of power.
#' @param sides Either 1 or 2 (default) to specify a one- or two- sided hypothesis test.
#' @param v Either TRUE for verbose output or FALSE (default) to output computed argument only.
#'
#' @return A list of the arguments (including the computed one).
#' @export
#'
#' @examples
#' multisite.bin(m = 30, J = 25, pc = 0.1, pt = 0.2, sigma.u = 0.4, power = NULL)
#' multisite.bin(m = 30, J = NULL, pc = 0.1, pt = 0.2, sigma.u = 0.4, power = 0.9)

multisite.bin <- function (m = NULL, alloc.ratio = 1, J = NULL,
                           pc = NULL, pt = NULL, sigma.u = NULL,
                           alpha = 0.05, power = NULL, sides = 2,
                           v = FALSE) {

  # Check if the arguments are specified correctly
  check.many(list(m, J, alpha, power), "oneof")
  check.param(m, "pos")
  check.param(alloc.ratio, "req"); check.param(alloc.ratio, "pos")
  check.param(J, "min", min = 2)
  check.param(pc, "req"); check.param(pc, "unit")
  check.param(pt, "req"); check.param(pt, "unit")
  check.param(sigma.u, "req"); check.param(sigma.u, "pos")
  check.param(alpha, "unit")
  check.param(power, "unit")
  check.param(sides, "req"); check.param(sides, "vals", valslist = c(1, 2))
  check.param(v, "req"); check.param(v, "bool")

  # Calculate power
  p.body <- quote({
    or <- (pt / (1 - pt)) / (pc / (1 - pc))
    gammaA <- abs(log(or))
    prop.t <- alloc.ratio / (1 + alloc.ratio)
    ssq.e <- (1 / 4) * (1 / ((1 - prop.t) * pc * (1 - pc)) + 1 /
                          (prop.t * pt * (1 - pt)))
    var <- 1.2 * (4 * (ssq.e / m + 0.25 * sigma.u^2)) / J
    za <- stats::qnorm(1 - alpha / sides)
    stats::pnorm(gammaA / sqrt(var) - za)
  })

  if (is.null(alpha)) {
    alpha <- stats::uniroot(function(alpha) eval(p.body) - power, c(1e-10, 1 - 1e-10))$root
    if (!v) return(alpha)
  }
  else if (is.null(power)) {
    power <- eval(p.body)
    if (!v) return(power)
  }
  else if (is.null(J)) {
    J <- stats::uniroot(function(J) eval(p.body) - power, c(2 + 1e-10, 1e+07))$root
    if (!v) return(J)
  }
  else if (is.null(m)) {
    m <- stats::uniroot(function(m) eval(p.body) - power, c(2 + 1e-10, 1e+07))$root
    if (!v) return(m)
  }
  # else if (is.null(alloc.ratio)) {
  #   alloc.ratio <- stats::uniroot(function(alloc.ratio) eval(p.body) - power + 1e-5, c(1 + 1e-10, 1e+07))$root
  #   if (!v) return(alloc.ratio)
  # }
  else stop("internal error")

  # Generate output text
  METHOD <-"Power for multisite trials with binary outcomes"
  NOTE <- "m1, m2 are the number of subjects within site in condition 1, condition 2\n      (total of m1 + m2 per site). m1, m2 OR m2, m1 produce equivalent results."
  p <- c(pc, pt)
  c <- round(m / (alloc.ratio + 1), 3)
  t <- round(alloc.ratio * c, 3)
  m <- c(t, c)

  # Print output as a power.htest object depending on which inputs were given
  structure(list(`m1, m2` = m, J = J, `pc, pt` = p, sigma.u = sigma.u,
                 alpha = alpha, power = power,
                 method = METHOD, note = NOTE), class = "power.htest")

}
