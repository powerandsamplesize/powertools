#' Power calculation for paired t test
#'
#' @description
#' Performs power and sample size calculations for a paired t test.
#' Can solve for power, N, delta or alpha.
#'
#'
#' @param N The sample size; if the observations are paired differences, this is the number of pairs.
#' @param delta DeltaA (the true mean difference) - Delta0 (the mean difference under the null).
#' @param sd1 The pre standard deviation; defaults to 1.
#' @param sd2 The post standard deviation; defaults to the value of sd1.
#' @param rho The correlation between pre and post measurements on the same individual.
#' @param alpha The significance level (type 1 error rate); defaults to 0.05.
#' @param power The specified level of power.
#' @param sides Either 1 or 2 (default) to specify a one- or two- sided hypothesis test.
#' @param v Either TRUE for verbose output or FALSE (default) to output computed argument only.
#'
#' @return A list of the arguments (including the computed one).
#' @export ttest.paired
#'
#' @examples
#' ttest.paired(N = NULL, delta = 4, sd1 = 10, sd2 = 10, rho = 0.4, power = 0.8, sides = 2)

ttest.paired <- function (N = NULL, delta = NULL,
                          sd1 = 1, sd2 = NULL, rho = NULL,
                          alpha = 0.05, power = NULL, sides = 2,
                          v = FALSE) {

  # Check if the arguments are specified correctly
  check.many(list(N, delta, alpha, power), "oneof")
  check.param(N, "pos")
  check.param(sd1, "req"); check.param(sd1, "pos")
  check.param(rho, "req"); check.param(rho, "uniti")
  check.param(delta, "num")
  check.param(alpha, "unit")
  check.param(power, "unit")
  check.param(sides, "req"); check.param(sides, "vals", valslist = c(1, 2))
  check.param(v, "req"); check.param(v, "bool")

  # If sd2 is missing, default to sd1
  if (is.null(sd2)) sd2 <- sd1

  # Calculate the standard deviation of differences within pairs
  sigmad <- sqrt(sd1^2 + sd2^2 - 2 * rho * sd1 * sd2)

  # Calculate df and ncp
  if (sides == 1)
    p.body <- quote({
      d <- abs(delta) / sigmad
      df <- N - 1
      stats::pt(stats::qt(alpha, df, lower.tail = FALSE), df,
                sqrt(N) * d, lower.tail = FALSE)
    })
  else if (sides == 2)
    p.body <- quote({
      d <- abs(delta)
      ncp <- sqrt(N) * d / sigmad
      df2 <- N - 1
      stats::pf(stats::qf(alpha, 1, df2, lower.tail = FALSE),
                1, df2, ncp^2, lower.tail = FALSE)
    })

  # Use stats::uniroot function to calculate missing argument
  if (is.null(power)) {
    power <- eval(p.body)
    if (!v) return(power)
  }
  else if (is.null(N)) {
    N <- stats::uniroot(function(N) eval(p.body) - power, c(2, 1e+07))$root
    if (!v) return(N)
  }
  else if (is.null(delta)) {
    delta <- stats::uniroot(function(delta) eval(p.body) - power,  c(1e-07, 1e+07))$root
    if (!v) return(delta)
  }
  else if (is.null(alpha)) {
    alpha <- stats::uniroot(function(alpha) eval(p.body) - power, c(1e-10, 1 - 1e-10))$root
    if (!v) return(alpha)
  }
  else stop("internal error")

  # Generate output text
  METHOD <- "Paired t test power calculation"
  NOTE <- "N is the number of pairs"
  sd <- c(sd1, sd2)

  # Print output as a power.htest object
  structure(list(N = N, delta = delta, `sd1, sd2` = sd, rho = rho,
                 alpha = alpha, power = power, sides = sides,
                 method = METHOD, note = NOTE), class = "power.htest")
}
