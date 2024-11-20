#' Power calculation for chi-square goodness-of-fit test
#'
#' @description
#' Performs sample size and power calculations for chi-square goodness-of-fit test,
#' which is used to test whether a sample of data arises from a population with a specific
#' discrete distribution.
#' This function can solve for power, total sample size or alpha.
#'
#'
#' @param p0vec Vector of probabilities for the specified population distribution. Must sum to 1.
#' @param p1vec Vector of expected probabilities for the sample. Must sum to 1.
#' @param N The total number of observations.
#' @param alpha The significance level (type 1 error rate); defaults to 0.05.
#' @param power The specified level of power.
#' @param v Either TRUE for verbose output or FALSE (default) to output computed argument only.
#'
#' @return A list of the arguments (including the computed one).
#' @export
#'
#' @examples
#' chisq.gof(p0vec = c(0.5, 0.3, 0.2), p1vec = c(0.7, 0.15, 0.15), N = 50)

chisq.gof <- function (p0vec = NULL, p1vec = NULL,
                       N = NULL, alpha = 0.05, power = NULL,
                       v = FALSE) {

  # Check if the arguments are specified correctly
  check.many(list(N, alpha, power), "oneof")
  check.param(p0vec, "req"); check.param(p0vec, "vec"); check.param(p0vec, "sum")
  check.param(p1vec, "req"); check.param(p1vec, "vec"); check.param(p1vec, "sum")
  check.param(N, "pos"); check.param(N, "min", min = 2)
  check.param(alpha, "unit")
  check.param(power, "unit")
  check.param(v, "req"); check.param(v, "bool")

  if (length(p0vec) != length(p1vec))
    stop("the two proportion vectors must have the same lengths")

  if (any(p0vec <= 0) | any(p0vec >= 1) | any(p1vec <= 0) | any(p1vec  >= 1))
    stop("all proportions must be between 0 and 1")

  if (sum(p0vec) != 1 | sum(p1vec) != 1)
    stop("proportions must sum to 1 across each vector")

  # Calculate effect size and df
  es <- sqrt(sum((p1vec - p0vec)^2 / p0vec))
  df <- length(p0vec) - 1

  # Calculate test statistic
  p.body <- quote({
    stats::pchisq(stats::qchisq(alpha, df, lower = FALSE),
                  df, N * es^2, lower = FALSE)
  })

  # Use stats::uniroot function to calculate missing argument
  if (is.null(power)) {
    power <- eval(p.body)
    if (!v) return(power)
  }
  else if (is.null(N)) {
    N <- stats::uniroot(function(N) eval(p.body) - power, c(1 + 1e-10, 1e+09))$root
    if (!v) return(N)
  }
  else if (is.null(alpha)) {
    alpha <- stats::uniroot(function(sig.level) eval(p.body) - power, c(1e-10, 1 - 1e-10))$root
    if (!v) return(alpha)
  }
  else stop("internal error")

  # Generate output text
  METHOD <- "Chi-square goodness-of-fit power calculation"
  structure(list(p1vec = p1vec, p0vec = p0vec, `w effect size` = es, N = N, alpha = alpha,
                 power = power, method = METHOD), class = "power.htest")
}
