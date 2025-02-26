#' Power for detecting treatment effect heterogeneity in an individually randomized trial with a continuous outcome
#'
#' @description
#' This function performs power and sample size calculations for detecting a treatment-by-covariate
#' interaction effect in a two-arm randomized trial
#' with a continuous outcome. Can solve for power, beta, n1 or n.ratio.
#'
#' @details
#' Shieh G (2009) Detecting interaction effects in moderated multiple regression with continuous variables:
#' power and sample size considerations. Organizational Research Methods 12(3):510-528.
#'
#' Yang S, Li F, Starks MA, Hernandez AF, Mentz RJ, Choudhury KR (2020) Sample size requirements for detecting
#' treatment effect heterogeneity in cluster randomized trials. Statistics in Medicine 39:4218-4237.
#'
#' @param n1 The sample size for group 1.
#' @param n.ratio The ratio n2/n1 between the sample sizes of two groups; defaults to 1 (equal group sizes).
#' @param beta The regression coefficient for the treatment-by-covariate interaction term.
#' @param sd.x The standard deviation of the covariate.
#' @param sd.yx The standard deviation of the outcome variable adjusting for the covariate.
#' @param alpha The significance level (type 1 error rate); defaults to 0.05.
#' @param power The specified level of power.
#' @param sides Either 1 or 2 (default) to specify a one- or two- sided hypothesis test.
#' @param v Either TRUE for verbose output or FALSE (default) to output computed argument only.
#'
#' @return A list of the arguments (including the computed one).
#' @export
#'
#' @examples
#' irt.hte(n1 = 540, n.ratio = 1, beta = 1, sd.x = 12.7, sd.yx = 71)

irt.hte <- function (n1 = NULL, n.ratio = 1, beta = NULL,
                     sd.x = NULL, sd.yx = NULL,
                     alpha = 0.05, power = NULL, sides = 2,
                     v = FALSE) {

  # Check if the arguments are specified correctly
  check.many(list(n1, n.ratio, beta, alpha, power), "oneof")
  check.param(n1, "pos")
  check.param(n.ratio, "pos")
  check.param(beta, "num")
  check.param(sd.x, "req"); check.param(sd.x, "pos")
  check.param(sd.yx, "req"); check.param(sd.yx, "pos")
  check.param(alpha, "unit")
  check.param(power, "unit")
  check.param(sides, "req"); check.param(sides, "vals", valslist = c(1, 2))
  check.param(v, "req"); check.param(v, "bool")

  # Calculate power
  p.body <- quote({
    sd.w <- sqrt(n.ratio) / (1 + n.ratio)
    N <- n1 + n1 * n.ratio
    A <- abs(beta) * sd.w * sd.x * sqrt(N) / sd.yx
    stats::pnorm(stats::qnorm(alpha / sides) + A)
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
    n.ratio <- safe.uniroot(function(n.ratio) eval(p.body) - power,c(2/n1, 1e+07))$root
    if (!v) return(n.ratio)
  }
  else if (is.null(beta)) {
    beta <- safe.uniroot(function(beta) eval(p.body) - power,  c(1e-07, 1e+07))$root
    if (!v) return(beta)
  }
  else if (is.null(alpha)) {
    alpha <- safe.uniroot(function(alpha) eval(p.body) - power, c(1e-10, 1 - 1e-10))$root
    if (!v) return(alpha)
  }
  else stop("internal error")

  # Generate output text
  METHOD <- "Power for test of HTE in an individually randomized trial"
  n <- c(n1, n1 * n.ratio)
  sd <- c(sd.x, sd.yx)

  # Print output as a power.htest object
  structure(list(`n1, n2` = n, beta = beta, `sd1, sd2` = sd, alpha = alpha,
                 power = power, sides = sides,
                 method = METHOD), class = "power.htest")
}
