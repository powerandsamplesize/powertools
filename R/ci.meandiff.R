#' Power calculation for precision analysis (confidence interval) for a difference between two means
#'
#' @description
#' Calculates the "power" of a confidence interval for a difference between two means, that is, the probability
#' of achieving a 100(1 - alpha) percent confidence interval with halfwidth not greater
#' than a specified value. This function can solve for power, n1, n.ratio or alpha.
#'
#' @details
#' The unconditional probability is the probability of obtaining the desired precision
#' (i.e., that the observed halfwidth does not exceed the desired halfwidth)
#' regardless of whether or not the confidence interval includes the true parameter value.
#' The conditional probability is the probability of both obtaining the desired precision and having
#' the interval include the true parameter value.
#'
#'
#' @param n1 The sample size for group 1.
#' @param n.ratio The ratio n2/n1 between the sample sizes of two groups; defaults to 1 (equal group sizes).
#' @param halfwidth The desired halfwidth for the difference in means.
#' @param sd The estimated standard deviation; defaults to 1. Equal SDs in each group are assumed.
#' @param alpha The significance level (type 1 error rate); defaults to 0.05.
#' @param power The specified level of power.
#' @param cond Specify using unconditional or conditional probability. Defaults to FALSE.
#' @param v Either TRUE for verbose output or FALSE (default) to output computed argument only.
#'
#' @return A list of the arguments (including the computed one).
#' @import PowerTOST
#' @export
#'
#' @examples
#' ci.meandiff(n1 = NULL, halfwidth = 0.25, power = 0.8)
#' ci.meandiff(n1 = 134, halfwidth = 0.25, cond = TRUE)

ci.meandiff <- function (n1 = NULL, n.ratio = 1, halfwidth = NULL, sd = 1,
                         alpha = 0.05, power = NULL, cond = FALSE,
                         v = FALSE) {

  # Check if the arguments are specified correctly
  check.many(list(n1, n.ratio, alpha, power), "oneof")
  check.param(n1, "pos")
  check.param(n.ratio, "pos")
  check.param(alpha, "unit")
  check.param(power, "unit")
  check.param(halfwidth, "req"); check.param(halfwidth, "num")
  check.param(sd, "req"); check.param(sd, "pos")
  check.param(cond, "req"); check.param(cond, "bool")
  check.param(v, "req"); check.param(v, "bool")

  d <- halfwidth / sd

  p.body <- quote({
    df <- n1 * (n.ratio + 1) - 2
    t <- stats::qt(1 - alpha / 2, df)
    b <- d * sqrt(n1 * n.ratio * df) / (t * sqrt(n.ratio + 1))
    if (cond) {
      uQ <- PowerTOST::OwensQ(nu = df, t = t, delta = 0, a = 0, b = b)
      lQ <- PowerTOST::OwensQ(nu = df, t = 0, delta = 0, a = 0, b = b)
      (2 / (1 - alpha)) * (uQ - lQ)
    } else {
      stats::pchisq(b^2, df)
    }
  })

  if (is.null(power)) {
    power <- eval(p.body)
    if (!v) return(power)
  }
  else if (is.null(n1)) {
    n1 <- stats::uniroot(function(n1) eval(p.body) - power, c(2, 1e+07))$root
    if (!v) return(n1)
  }
  else if (is.null(n.ratio)) {
    n.ratio <- stats::uniroot(function(n.ratio) eval(p.body) - power, c(2/n1, 1e+07))$root
    if (!v) return(n.ratio)
  }
  else if (is.null(alpha)) {
    alpha <- stats::uniroot(function(alpha) eval(p.body) - power, c(1e-10, 1 - 1e-10))$root
    if (!v) return(alpha)
  }
  else stop("internal error")

  # Generate output text
  METHOD <- paste0("Precision analysis for difference between means\n     using ",
                   ifelse(cond, "", "un"), "conditional probability")
  n <- c(n1, n1 * n.ratio)

  # Print output as a power.htest object
  structure(list(n = n, halfwidth = halfwidth, sd = sd,
                 alpha = alpha, power = power,
                 method = METHOD), class = "power.htest")

}
