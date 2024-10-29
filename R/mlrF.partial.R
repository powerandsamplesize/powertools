#' Power calculation for a multiple linear regression partial F test
#'
#' @description
#' Conducts power and sample size calculations for a partial F test
#' in a multiple linear regression model.
#' This is a test that one or more coefficients are equal to zero
#' after controlling for a set of control predictors.
#' Can solve for power, N or alpha.
#'
#'
#'
#'
#' @param N The sample size.
#' @param p The number of control predictors.
#' @param q The number of test predictors.
#' @param Rsq.red The squared population multiple correlation coefficient for the reduced model. Either both Rsq terms OR pc must be specified.
#' @param Rsq.full The squared population multiple correlation coefficient for the full model. Either both Rsq terms OR pc must be specified.
#' @param pc The partial correlation coefficient. Either both Rsq terms OR pc must be specified.
#' @param alpha The significance level or type 1 error rate; defaults to 0.05.
#' @param power The specified level of power.
#' @param v Either TRUE for verbose output or FALSE to output computed argument only.
#'
#' @return A list of the arguments (including the computed one).
#' @export
#'
#' @examples
#' mlrF.partial(N = 80, p = 3, q = 2, Rsq.red = 0.25, Rsq.full = 0.35)
#' mlrF.partial(N = 150, p = 4, pc = 0.2)

mlrF.partial <- function (N = NULL, p = NULL, q = NULL, pc = NULL,
                          Rsq.red = NULL, Rsq.full = NULL,
                          alpha = 0.05, power = NULL, v = FALSE) {

  # Check if the arguments are specified correctly
  if ((is.null(pc) & (is.null(p) | is.null(q))) | (!is.null(pc) & is.null(p)))
    stop("please specify the number of predictors")
  if ((is.null(Rsq.red) | is.null(Rsq.full)) & is.null(pc))
    stop("please specify Rsq.red and Rsq.full OR pc")

  check.many(list(N, alpha, power), "oneof")
  check.param(N, "pos"); check.param(N, "min", min = 7)
  check.param(alpha, "unit")
  check.param(power, "unit")
  check.param(p, "int")
  check.param(q, "int")
  check.param(Rsq.red, "unit")
  check.param(Rsq.full, "unit")
  check.param(pc, "uniti")
  check.param(v, "req"); check.param(v, "bool")

  # Calculate power
  if (is.null(pc)) {
    p.body <- quote({
      ncp <- N * (Rsq.full - Rsq.red) / (1 - Rsq.full)
      df2 <- N - p - q - 1
      crit <- stats::qf(1 - alpha, q, df2)
      1 - stats::pf(crit, q, df2, ncp)
    })
  } else {
    p.body <- quote({
      ncp <- N * pc^2 / (1 - pc^2)
      df2 <- N - p - 2
      crit <- stats::qf(1 - alpha, 1, df2)
      1 - stats::pf(crit, 1, df2, ncp)
    })
  }

  # Use stats::uniroot function to calculate missing argument
  if (is.null(power)) {
    power <- eval(p.body)
    if (!v) return(power)
  }
  else if (is.null(N)) {
    N <- stats::uniroot(function(n) eval(p.body) - power, c(7, 1e+09))$root
    if (!v) return(N)
  }
  else if (is.null(alpha)) {
    alpha <- stats::uniroot(function(alpha) eval(p.body) - power, c(1e-10, 1 - 1e-10))$root
    if (!v) return(alpha)
  }
  else stop("internal error")

  # Generate output text
  METHOD <- "Power calculation for a multiple linear regression\n     partial F test"

  # Print output as a power.htest object
  if (is.null(pc)) {
    structure(list(N = N, p = p, q = q,
                   Rsq.red = Rsq.red, Rsq.full = Rsq.full,
                   alpha = alpha, power = power,
                   method = METHOD), class = "power.htest")
  } else {
    structure(list(N = N, p = p, pc = pc,
                   alpha = alpha, power = power,
                   method = METHOD), class = "power.htest")
  }

}
