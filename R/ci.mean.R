#' Power calculation for precision analysis (confidence interval) for one mean
#'
#' @description
#' Calculates the "power" of a confidence interval for one mean, that is, the probability
#' of achieving a 100(1 - alpha) percent confidence interval with halfwidth not greater
#' than a specified value. Can solve for power, N or alpha.
#'
#' @details
#' The unconditional probability is the probability of obtaining the desired precision
#' (i.e., that the observed halfwidth does not exceed the desired halfwidth)
#' regardless of whether or not the confidence interval includes the true parameter value.
#' The conditional probability is the probability of both obtaining the desired precision and having
#' the interval include the true parameter value.
#'
#'
#'
#' @param N The sample size.
#' @param halfwidth The desired halfwidth.
#' @param sd The estimated standard deviation; defaults to 1.
#' @param alpha The significance level (type 1 error rate); defaults to 0.05.
#' @param power The specified level of power.
#' @param cond Specify whether to use unconditional or conditional probability. Defaults to FALSE (unconditional).
#' @param v Either TRUE for verbose output or FALSE (default) to output computed argument only.
#'
#' @return A list of the arguments (including the computed one).
#' @import PowerTOST
#' @export
#'
#' @examples
#' ci.mean(N = NULL, halfwidth = 0.25, power = 0.8)
#' ci.mean(N = 62, halfwidth = 0.25, power = NULL)
#' ci.mean(N = 73, halfwidth = 0.25, cond = TRUE)

ci.mean <- function (N = NULL, halfwidth = NULL, sd = 1,
                     alpha = 0.05, power = NULL, cond = FALSE,
                     v = FALSE) {

  # Check if the arguments are specified correctly
  check.many(list(N, alpha, power), "oneof")
  check.param(N, "pos")
  check.param(alpha, "unit")
  check.param(power, "unit")
  check.param(halfwidth, "req"); check.param(halfwidth, "num")
  check.param(sd, "req"); check.param(sd, "pos")
  check.param(cond, "req"); check.param(cond, "bool")
  check.param(v, "req"); check.param(v, "bool")

  d <- halfwidth / sd

  p.body <- quote({
    df <- N - 1
    t <- stats::qt(1 - alpha / 2, df)
    b <- d * sqrt(N * df) / t
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
  else if (is.null(N)) {
    N <- stats::uniroot(function(N) eval(p.body) - power, c(2, 1e+07))$root
    if (!v) return(N)
  }
  else if (is.null(alpha)) {
    alpha <- stats::uniroot(function(alpha) eval(p.body) - power, c(1e-10, 1 - 1e-10))$root
    if (!v) return(alpha)
  }
  else stop("internal error")

  # Generate output text
  METHOD <- paste0("Precision analysis for one mean\n     using ",
                   ifelse(cond, "", "un"), "conditional probability")

  # Print output as a power.htest object
  structure(list(N = N, halfwidth = halfwidth, sd = sd,
                 alpha = alpha, power = power,
                 method = METHOD), class = "power.htest")

}
