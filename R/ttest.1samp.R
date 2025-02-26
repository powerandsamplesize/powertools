#' Power calculation for one-sample t test
#'
#' @description
#' This function computes power and sample size for a one-sample t test.
#' Can solve for power, N, delta, sd or alpha.
#'
#'
#' @param N The sample size.
#' @param delta muA (the true mean) - mu0 (the mean under the null).
#' @param sd The standard deviation; defaults to 1.
#' @param alpha The significance level (type 1 error rate); defaults to 0.05.
#' @param power The specified level of power.
#' @param sides Either 1 or 2 (default) to specify a one- or two-sided hypothesis test.
#' @param v Either TRUE for verbose output or FALSE (default) to output computed argument only.
#'
#' @return A list of the arguments (including the computed one).
#' @export ttest.1samp
#'
#' @examples
#' ttest.1samp(N = 36, delta = 4.9 - 5.7, sd = 2, sides = 1)
#' ttest.1samp(N = 36, delta = 6.3 - 5.7, sd = 2, sides = 1)
#' ttest.1samp(N = 36, delta = 4.9 - 5.7, sd = 2, sides = 2)
#' ttest.1samp(delta = 0.6, sd = 1, power = 0.8, sides = 1)

ttest.1samp <- function (N = NULL, delta = NULL, sd = 1,
                         alpha = 0.05, power = NULL, sides = 2,
                         v = FALSE) {

  # Check if the arguments are specified correctly
  check.many(list(N, delta, sd, alpha, power), "oneof")
  check.param(N, "pos")
  check.param(delta, "num")
  check.param(sd, "pos")
  check.param(alpha, "unit")
  check.param(power, "unit")
  check.param(sides, "req"); check.param(sides, "vals", valslist = c(1, 2))
  check.param(v, "req"); check.param(v, "bool")

  # Calculate df and ncp
  if (sides == 1)
    p.body <- quote({
      d <- abs(delta) / sd
      df <- N - 1
      stats::pt(stats::qt(alpha, df, lower.tail = FALSE), df,
                sqrt(N) * d, lower.tail = FALSE)
    })
  else if (sides == 2)
    p.body <- quote({
      d <- abs(delta) / sd
      ncp <- sqrt(N) * d
      df2 <- N - 1
      stats::pf(stats::qf(alpha, 1, df2, lower.tail = FALSE),
                1, df2, ncp^2, lower.tail = FALSE)
    })

  # Use safe.uniroot function to calculate missing argument
  if (is.null(power)) {
    power <- eval(p.body)
    if (!v) return(power)
  }
  else if (is.null(N)) {
    N <- safe.uniroot(function(N) eval(p.body) - power, c(2, 1e+07))$root
    if (!v) return(N)
  }
  else if (is.null(sd)) {
    sd <- safe.uniroot(function(sd) eval(p.body) - power, delta * c(1e-07, 1e+07))$root
    if (!v) return(sd)
  }
  else if (is.null(delta)) {
    delta <- safe.uniroot(function(delta) eval(p.body) - power,  c(1e-07, 1e+07))$root
    if (!v) return(delta)
  }
  else if (is.null(alpha)) {
    alpha <- safe.uniroot(function(alpha) eval(p.body) - power, c(1e-10, 1 - 1e-10))$root
    if (!v) return(alpha)
  }

  else stop("internal error")

  # Generate output text
  METHOD <- "One-sample t test power calculation"

  # Print output as a power.htest object
  structure(list(N = N, delta = delta, sd = sd, alpha = alpha,
                 power = power, sides = sides,
                 method = METHOD), class = "power.htest")

}
