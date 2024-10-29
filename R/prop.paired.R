#' Power calculation for McNemar test of two correlated proportions
#'
#' @description
#' Performs power and sample size calculation for McNemar test of two correlated
#' proportions using normal approximation. Can solve for power, N or alpha.
#'
#' @details
#' Either p1, p2 and phi, OR paid and dpr must be specified.
#'
#'
#'
#'
#' @param N The sample size; the number of pairs.
#' @param p1 The outcome proportion under condition 1.
#' @param p2 The outcome proportion under condition 2.
#' @param phi The correlation between the two responses from an individual.
#' @param paid The smaller of the two discordant probabilities. Either p1, p2 and phi, OR paid and dpr must be specified.
#' @param dpr The discordant proportion ratio. Must be greater than or equal to 1. Either p1, p2, and phi, OR paid and dpr must be specified.
#' @param alpha The significance level (type 1 error rate); defaults to 0.05.
#' @param power The specified level of power.
#' @param sides Either 1 or 2 (default) to specify a one- or two- sided hypothesis test.
#' @param v Either TRUE for verbose output or FALSE (default) to output computed argument only.
#'
#' @return A list of the arguments (including the computed one).
#' @export
#'
#' @examples
#' prop.paired(N = NULL, p1 = 0.8, p2 = 0.9, phi = 0, power = 0.9, sides = 2)
#' prop.paired(N = NULL, paid = 0.08, dpr = 0.18 / 0.08, power = 0.9, sides = 2)

prop.paired <- function (N = NULL, p1 = NULL, p2 = NULL, phi = NULL,
                         paid = NULL, dpr = NULL, alpha = 0.05,
                         power = NULL, sides = 2, v = FALSE) {

  # Check if the arguments are specified correctly
  if ((is.null(p1) | is.null(p2) | is.null(phi)) & (is.null(dpr) | is.null(paid)))
    stop("p1, p2, and phi OR dpr and paid must be specified")
  check.many(list(N, alpha, power), "oneof")
  check.param(N, "pos")
  check.param(p1, "unit")
  check.param(p2, "unit")
  check.param(phi, "uniti")
  check.param(paid, "unit")
  check.param(dpr, "min", min = 1)
  check.param(alpha, "unit")
  check.param(power, "unit")
  check.param(sides, "req"); check.param(sides, "vals", valslist = c(1, 2))
  check.param(v, "req"); check.param(v, "bool")

  # Calculate paid and dpr if not given
  if (is.null(paid) & is.null(dpr)) {
    p01 <- p1 * (1 - p2) - phi * sqrt((1 - p2) * p1 * (1 - p1) * p2)
    p10 <- p01 + p2 - p1
    paid <- p01
    dpr <- p10 / p01
  }

  # Calculate test statistic
  p.body <- quote(stats::pnorm((sqrt(N * paid) * (dpr - 1) -
                  stats::qnorm(alpha / sides, lower.tail = FALSE) *
                  sqrt(dpr + 1)) / sqrt((dpr + 1) - paid * (dpr - 1)^2)))

  NOTE <- "N is the number of pairs"
  if (!v) print(paste("NOTE:", NOTE))

  # Use stats::uniroot function to calculate missing argument
  if (is.null(power)) {
    power <- eval(p.body)
    if (!v) return(power)
  }
  else if (is.null(N)) {
    N <- stats::uniroot(function(N) eval(p.body) - power, c(ceiling(log(alpha) / log(0.5)), 1e+07))$root
    if (!v) return(N)
  }
  else if (is.null(alpha)) {
    alpha <- stats::uniroot(function(alpha) eval(p.body) - power, c(1e-10, 1 - 1e-10))$root
    if (!v) return(alpha)
  }
  else stop("internal error", domain = NA)

  # Generate output text
  METHOD <-"McNemar paired comparison of proportions\n     approximate power calculation"

  # Print output as a power.htest object depending on which inputs were given
  if (!is.null(p1) & !is.null(p2) & !is.null(phi)) {
    p <- c(p1, p2)
    structure(list(N = N, `p1, p2` = p, phi = phi, alpha = alpha,
                   power = power, sides = sides, note = NOTE,
                   method = METHOD), class = "power.htest")
  }
  else if (!is.null(paid) & !is.null(dpr))
    structure(list(N = N, paid = paid, dpr = dpr, alpha = alpha,
                   power = power, sides = sides, note = NOTE,
                   method = METHOD), class = "power.htest")
}

