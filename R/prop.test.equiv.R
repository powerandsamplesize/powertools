#' Power calculation for test of equivalence of two proportions
#'
#' @description
#' Performs power and sample size calculations for a test of equivalence for two
#' proportions. Can solve for power, n1, n.ratio or alpha.
#'
#'
#' @param n1 The sample size for group 1.
#' @param n.ratio The ratio n2/n1 between the sample sizes of two groups; defaults to 1 (equal group sizes).
#' @param p1 The outcome proportion in group 1.
#' @param p2 The outcome proportion in group 2.
#' @param margin The equivalence margin.
#' @param alpha The significance level or type 1 error rate; defaults to 0.05.
#' @param power The specified level of power.
#' @param v Either TRUE for verbose output or FALSE to output computed argument only.
#'
#' @return A list of the arguments (including the computed one).
#' @export
#'
#' @examples
#' prop.test.equiv(n1 = NULL, p1 = 0.5, p2 = 0.5, margin = 0.1, alpha = 0.05, power = 0.8)

prop.test.equiv <- function (n1 = NULL, n.ratio = 1, p1 = NULL, p2 = NULL, margin = NULL,
                             alpha = 0.05, power = NULL, v = FALSE) {

  # Check if the arguments are specified correctly
  check.many(list(n1, n.ratio, alpha, power), "oneof")
  check.param(n1, "pos")
  check.param(n.ratio, "pos")
  check.param(p1, "req"); check.param(p1, "unit")
  check.param(p2, "req"); check.param(p2, "unit")
  check.param(margin, "req"); check.param(margin, "num")
  check.param(alpha, "unit")
  check.param(power, "unit")
  check.param(v, "req"); check.param(v, "bool")

  # Calculate n1
  p.body <- quote({
    d <- abs(p1 - p2)
    beta <- 1 - power
    ((stats::qnorm(1 - alpha) + stats::qnorm(1 - beta / 2))^2 *
    (p1 * (1 - p1) + p2 * (1 - p2) / n.ratio) / (margin - d)^2)
  })

  # Use safe.uniroot function to calculate missing argument
  if (is.null(n1)) {
    n1 <- eval(p.body)
    if (!v) return(n1)
  }
  else if (is.null(power)) {
    power <- safe.uniroot(function(power) eval(p.body) - n1, c(1e-05, 0.99999))$root
    if (!v) return(power)
  }
  else if (is.null(n.ratio)) {
    n.ratio <- safe.uniroot(function(n.ratio) eval(p.body) - n1, c(2/n1, 1e+07))$root
    if (!v) return(n.ratio)
  }
  else if (is.null(alpha)) {
    alpha <- safe.uniroot(function(alpha) eval(p.body) - n1, c(1e-10, 1 - 1e-10))$root
    if (!v) return(alpha)
  }
  else stop("internal error", domain = NA)

  # Generate output text
  METHOD <- "Test for equivalence of two proportions power calculation"
  n <- c(n1, n1 * n.ratio)
  p <- c(p1, p2)

  # Print output as a power.htest object
  structure(list(`n1, n2` = n, `p1, p2` = p,
                 margin = margin, alpha = alpha,
                 power = power,
                 method = METHOD), class = "power.htest")
}

