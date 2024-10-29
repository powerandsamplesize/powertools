#' Power calculation for a simple linear regression
#'
#' @description
#' Performs sample size and power calculations for a simple linear regression.
#' Can solve for power, N or alpha. Required inputs include the values of the
#' slope regression coefficient under the null hypothesis and the alternative,
#' the variance of the covariate X, and the SD of the error. Power calculations
#' for a simple linear regression can also be conducted using the mlrF.overall function,
#' which requires fewer inputs.
#'
#' @param N The sample size.
#' @param beta10 The slope regression coefficient under the null hypothesis; defaults to 0.
#' @param beta1A The slope regression coefficient under the alternative hypothesis.
#' @param var.x The variance of the covariate X.
#' @param sigma.e The standard deviation of the error terms.
#' @param alpha The significance level (type 1 error rate); defaults to 0.05.
#' @param power The specified level of power.
#' @param sides Either 1 or 2 (default) to specify a one- or two- sided hypothesis test.
#' @param v Either TRUE for verbose output or FALSE (default) to output computed argument only.
#'
#' @return A list of the arguments (including the computed one).
#' @export
#'
#' @examples
#' # Yi = beta0 + beta1 * Xi + ei, i = 1,...,N
#' slr(N = 100, beta10 = 1, beta1A = 1.5, var.x = 25, sigma.e = 10, sides = 1)

slr <- function (N = NULL, beta10 = 0, beta1A = NULL,
                 var.x = NULL, sigma.e = NULL,
                 alpha = 0.05, power = NULL, sides = 2,
                 v = FALSE) {

  # Check if the arguments are specified correctly
  check.many(list(N, alpha, power), "oneof")
  check.param(N, "pos"); check.param(N, "min", min = 4)
  check.param(alpha, "unit")
  check.param(power, "unit")
  check.param(beta10, "req"); check.param(beta10, "num")
  check.param(beta1A, "req"); check.param(beta1A, "num")
  check.param(var.x, "req"); check.param(var.x, "pos")
  check.param(sigma.e, "req"); check.param(sigma.e, "pos")
  check.param(sides, "req"); check.param(sides, "vals", valslist = c(1, 2))
  check.param(v, "req"); check.param(v, "bool")

  # Calculate power
  if (sides == 1)
    p.body <- quote({
      sxx <- (N - 1) * var.x
      ncp <- abs(beta1A - beta10) * sqrt(sxx) / sigma.e
      df <- N - 2
      1 - stats::pt(stats::qt(1 - alpha, df), df, ncp)
    })
  else if (sides == 2)
    p.body <- quote({
      sxx <- (N - 1) * var.x
      ncp <- abs(beta1A - beta10) * sqrt(sxx) / sigma.e
      df2 <- N - 2
      1 - stats::pf(stats::qf(1 - alpha, 1, df2), 1, df2, ncp^2)
    })

  # Use stats::uniroot function to calculate missing argument
  if (is.null(power)) {
    power <- eval(p.body)
    if (!v) return(power)
  }
  else if (is.null(N)) {
    N <- stats::uniroot(function(N) eval(p.body) - power, c(4, 1e+09))$root
    if (!v) return(N)
  }
  else if (is.null(alpha)) {
    alpha <- stats::uniroot(function(alpha) eval(p.body) - power, c(1e-10, 1 - 1e-10))$root
    if (!v) return(alpha)
  }
  else stop("internal error")

  # Generate output text
  METHOD <- "Power calculation for a simple linear regression"

  # Print output as a power.htest object
  structure(list(N = N, beta10 = beta10, beta1A = beta1A,
                 var.x = var.x, sigma.e = sigma.e,
                 alpha = alpha, power = power, sides = sides,
                 method = METHOD), class = "power.htest")

}
