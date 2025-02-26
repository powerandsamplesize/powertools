#' Power calculations for ordinal categorical variable under proportional odds assumption
#'
#' @description
#' Performs power and sample size calculation for a comparison of two groups on
#' an ordinal categorical response variable.
#' Assumes that response probabilities follow the proportional
#' odds assumption. Can solve for power, n1, n.ratio and alpha.
#'
#' @details
#' Whitehead J. (1993) Sample size calculations for ordered categorical data.
#' Statistics in Medicine, 12(24):2257–2271
#'
#' @param pC Vector of response probabilities in control group (group 1). Must sum to 1. Categories are ordered from best to worst.
#' @param OR Odds ratio when the alternative is true. Must be greater than 1.
#' @param n1 Sample size for group 1 (control group).
#' @param n.ratio The ratio n2/n1 between the sample sizes of two groups; defaults to 1 (equal group sizes).
#' @param alpha The significance level (type 1 error rate); defaults to 0.05.
#' @param power The specified level of power; defaults to 0.8.
#' @param v Either TRUE for verbose output or FALSE (default) to output computed argument only.
#' @import Hmisc
#' @return A list of the arguments (including the computed one).
#' @export
#'
#' @examples
#' library(Hmisc)
#' pC <- c(0.2, 0.5, 0.2, 0.1)
#' propodds(pC = pC, OR = 2.5, n1 = 65, n.ratio = 1, alpha = 0.05)

propodds <- function(pC, OR, n1, n.ratio = 1, alpha = 0.05,
                     power = NULL, v = FALSE){

  # Check if the arguments are specified correctly
  check.many(list(n1, n.ratio, power, alpha), "oneof")
  check.param(n1, "pos"); check.param(n1, "min", min = 2)
  check.param(n.ratio, "pos")
  check.param(power, "unit")
  check.param(alpha, "unit")
  check.param(pC, "req"); check.param(pC, "sum")
  check.param(OR, "req"); check.param(OR, "mini", min = 1)
  check.param(v, "req"); check.param(v, "bool")

  pC <- pC[!is.na(pC)]
  if (abs(sum(pC) - 1) > 1e-05)
    stop("probabilities in pC do not sum to 1")
  pT <- Hmisc::pomodm(p = pC, odds.ratio = 1 / OR)
  pavg <- (pC + pT) / 2
  ps <- 1 - sum(pavg^3)

  p.body <- quote({
    N <- n1 * (1 + n.ratio)
    n2 <- n1 * n.ratio
    V <- n1 * n2 * N/3/((N + 1)^2) * ps
    za <- stats::qnorm(1 - alpha/2)
    stats::pnorm(abs(log(OR)) * sqrt(V) - za)
  })

  if (is.null(power)){
    power <- eval(p.body)
    if (!v) return(power)
  }
  else if (is.null(n1)){
    n1 <- safe.uniroot(function(n1) eval(p.body) - power, c(2 + 1e-10, 1e+09))$root
    if (!v) return(n1)
  }
  else if (is.null(n.ratio)) {
    n.ratio <- safe.uniroot(function(n.ratio) eval(p.body) - power, c(2/n1, 1e+07))$root
    if (!v) return(n.ratio)
  }
  else if (is.null(alpha)) {
    alpha <- safe.uniroot(function(alpha) eval(p.body) - power, c(1e-10, 1 - 1e-10))$root
    if (!v) return(alpha)
  }
  else stop("internal error")

  METHOD = "Power calculation for ordinal categorical outcome under proportional odds"
  n <- c(round(n1, 3), round(n1 * n.ratio, 3))
  structure(list(n = n, pC = pC, pT = round(pT, digits = 3),
                 OR = OR, alpha = alpha, power = power, method = METHOD),
            class = "power.htest")
}


