#' Power calculation for comparing two correlation coefficients
#'
#' @description
#' Calculates power and sample size for a test that the correlation coefficients in two groups/populations
#' are equal. Can solve for power, n1, n.ratio or alpha.
#'
#'
#' @param n1 The sample size for group 1.
#' @param n.ratio The ratio n2/n1 between the sample sizes of two groups; defaults to 1 (equal group sizes).
#' @param rho1 The correlation coefficient in group 1.
#' @param rho2 The correlation coefficient in group 2.
#' @param alpha The significance level (type 1 error rate); defaults to 0.05.
#' @param power The specified level of power.
#' @param sides Either 1 or 2 (default) to specify a one- or two- sided hypothesis test.
#' @param v Either TRUE for verbose output or FALSE (default) to output computed argument only.
#'
#' @return A list of the arguments (including the computed one).
#' @export
#'
#' @examples
#' corr.2samp(n1 = 300, rho1 = 0.3, rho2 = 0.1, sides = 1)

corr.2samp <- function (n1 = NULL, n.ratio = 1, rho1 = NULL, rho2 = NULL,
                        alpha = 0.05, power = NULL, sides = 2,
                        v = FALSE) {

  # Check if the arguments are specified correctly
  check.many(list(n1, n.ratio, alpha, power), "oneof")
  check.param(n1, "pos"); check.param(n1, "min", min = 4)
  check.param(n.ratio, "pos")
  check.param(alpha, "unit")
  check.param(power, "unit")
  check.param(rho1, "req"); check.param(abs(rho1), "uniti")
  check.param(rho2, "req"); check.param(abs(rho2), "uniti")
  check.param(sides, "req"); check.param(sides, "vals", valslist = c(1, 2))
  check.param(v, "req"); check.param(v, "bool")

  fisherz <- function(rho){
    0.5 * log((1 + rho)/(1 - rho))
  }

  # Calculate power
  p.body <- quote({
    za <- stats::qnorm(1 - alpha / sides)
    f1 <- fisherz(rho1) + rho1 / (2 * (n1 - 1))
    f2 <- fisherz(rho2) + rho2 / (2 * (n1 * n.ratio - 1))
    DeltaA <- abs(f1 - f2)
    lambda <- DeltaA / sqrt(1 / (n1 - 3) + 1 / (n1 * n.ratio - 3))

    stats::pnorm(abs(lambda) - za)
  })

  # Use safe.uniroot function to calculate missing argument
  if (is.null(power)) {
    power <- eval(p.body)
    if (!v) return(power)
  }
  else if (is.null(n1)) {
    n1 <- safe.uniroot(function(n1) eval(p.body) - power, c(4 + 1e-10, 1e+09))$root
    if (!v) return(n1)
  }
  else if (is.null(n.ratio)) {
    n.ratio <- safe.uniroot(function(n.ratio) eval(p.body) - power, c(4/n1, 1e+07))$root
    if (!v) return(n.ratio)
  }
  else if (is.null(alpha)) {
    alpha <- safe.uniroot(function(alpha) eval(p.body) - power, c(1e-7, 1 - 1e-7))$root
    if (!v) return(alpha)
  }
  else stop("internal error")

  # Generate output text
  METHOD <- "Power calculation for comparing two correlation coefficients"
  n <- c(n1, n1 * n.ratio)

  # Print output as a power.htest object
  structure(list(n = n, rho1 = rho1, rho2 = rho2, alpha = alpha, power = power,
                 method = METHOD), class = "power.htest")


}
