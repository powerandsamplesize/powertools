#' Power calculation for chi-square test of independence
#'
#' @description
#' Performs power and sample size calculations for a chi-square test of independence.
#' The user inputs a matrix of cell probabilities for a two-way table. The function computes
#' the power (or required total sample size) for a test of no association between the two factors.
#'
#'
#'
#' @param pmatrix The two-way probability table under the alternative hypothesis. The probabilities must sum to 1.
#' @param N The total number of observations.
#' @param alpha The significance level (type 1 error rate); defaults to 0.05.
#' @param power The specified level of power.
#' @param v Either TRUE for verbose output or FALSE (default) to output computed argument only.
#'
#' @return A list of the arguments (including the computed one).
#' @export
#'
#' @examples
#' chisq.indep(pmatrix = matrix(c(0.050, 0.350, 0.100, 0.075, 0.250, 0.175),
#' nrow = 2, byrow = TRUE), N = 230)
#' chisq.indep(pmatrix = matrix(c(0.3, 0.2, 0.4, 0.1), nrow = 2, byrow = TRUE), N = 200)

chisq.indep <- function (pmatrix = NULL, N = NULL, alpha = 0.05, power = NULL,
                         v = FALSE) {

  # Check if the arguments are specified correctly
  check.many(list(N, alpha, power), "oneof")
  check.param(pmatrix, "req"); check.param(pmatrix, "mat"); check.param(pmatrix, "sum")
  check.param(N, "pos"); check.param(N, "min", min = 2)
  check.param(alpha, "unit")
  check.param(power, "unit")
  check.param(v, "req"); check.param(v, "bool")

  # Calculate effect size and df
  pi <- apply(pmatrix, 1, sum)
  pj <- apply(pmatrix, 2, sum)
  p0 <- pi %*% t(pj)
  es <- sqrt(sum((pmatrix - p0)^2/p0))
  df <- (nrow(pmatrix) - 1) * (ncol(pmatrix) - 1)

  # Calculate test statistic
  p.body <- quote({
    stats::pchisq(stats::qchisq(alpha, df, lower = FALSE),
                  df, N * es^2, lower = FALSE)
  })

  # Use safe.uniroot function to calculate missing argument
  if (is.null(power)) {
    power <- eval(p.body)
    if (!v) return(power)
  }
  else if (is.null(N)) {
    N <- safe.uniroot(function(N) eval(p.body) - power, c(1 + 1e-10, 1e+09))$root
    if (!v) return(N)
  }
  else if (is.null(alpha)) {
    alpha <- safe.uniroot(function(sig.level) eval(p.body) - power, c(1e-10, 1 - 1e-10))$root
    if (!v) return(alpha)
  }
  else stop("internal error")

  # Generate output text
  METHOD <- "Chi-square test of independence power calculation"

  # Print output as a power.htest object
  structure(list(pmatrix = matrix.format(pmatrix),
                 `w effect size` = es, N = N, alpha = alpha,
                 power = power, method = METHOD), class = "power.htest")
}
