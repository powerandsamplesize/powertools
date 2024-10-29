#' Power calculation for test of one correlation coefficient
#'
#' @description
#' Calculates power and sample size for a test that the correlation coefficient in a single population
#' is equal to (or less than or greater than) a specified value. Can solve for power, N or alpha.
#'
#'
#' @param N The sample size.
#' @param rho0 The correlation coefficient under the null hypothesis; defaults to 0.
#' @param rhoA The correlation coefficient under the alternative hypothesis.
#' @param alpha The significance level (type 1 error rate); defaults to 0.05.
#' @param power The specified level of power.
#' @param sides Either 1 or 2 (default) to specify a one- or two- sided hypothesis test.
#' @param v Either TRUE for verbose output or FALSE (default) to output computed argument only.
#'
#' @return A list of the arguments (including the computed one).
#' @export
#'
#' @examples
#' corr.1samp(N = 100, rhoA = 0.2, sides = 1)
#' corr.1samp(N = 100, rho0 = 0.2, rhoA = 0.4, sides = 1)

corr.1samp <- function (N = NULL, rho0 = 0, rhoA = NULL,
                        alpha = 0.05, power = NULL, sides = 2,
                        v = FALSE) {

  # Check if the arguments are specified correctly
  check.many(list(N, alpha, power), "oneof")
  check.param(N, "pos"); check.param(N, "min", min = 4)
  check.param(alpha, "unit")
  check.param(power, "unit")
  check.param(rho0, "req"); check.param(abs(rho0), "uniti")
  check.param(rhoA, "req"); check.param(abs(rhoA), "uniti")
  check.param(sides, "req"); check.param(sides, "vals", valslist = c(1, 2))
  check.param(v, "req"); check.param(v, "bool")

  # Calculate power
  p.body <- quote({
    za <- stats::qnorm(1 - alpha / sides)
    K1 <- (5 + rhoA^2) / (4 * (N - 1))
    K2 <- (11 + 2 * rhoA^2 + 3 * rhoA^4) / (8 * (N - 1)^2)
    K <- 1 + K1 + K2

    om1 <- 0.5 * log((1 + rhoA) / (1 - rhoA))
    om2 <- K * rhoA / (2 * (N - 1))
    om3 <- 0.5 * log((1 + rho0) / (1 - rho0)) + rho0 / (2 * (N - 1))
    omega <- sqrt(N - 3) * (om1 + om2 - om3)

    nu1 <- (4 - rhoA^2) / (2 * (N - 1))
    nu2 <- (22 - 6 * rhoA^2 - 3 * rhoA^4) / (6 * (N - 1)^2)
    nu <- ((N - 3) / (N - 1)) * (1 + nu1 + nu2)

    stats::pnorm((abs(omega) - za) / sqrt(nu))
  })

  # Use stats::uniroot function to calculate missing argument
  if (is.null(power)) {
    power <- eval(p.body)
    if (!v) return(power)
  }
  else if (is.null(N)) {
    N <- stats::uniroot(function(N) eval(p.body) - power, c(4 + 1e-10, 1e+09))$root
    if (!v) return(N)
  }
  else if (is.null(alpha)) {
    alpha <- stats::uniroot(function(alpha) eval(p.body) - power, c(1e-10, 1 - 1e-10))$root
    if (!v) return(alpha)
  }
  else stop("internal error")

  # Generate output text
  METHOD <- "Single correlation coefficient power calculation"

  # Print output as a power.htest object
  structure(list(N = N, rho0 = rho0, rhoA = rhoA, alpha = alpha, power = power,
                 method = METHOD), class = "power.htest")

}
