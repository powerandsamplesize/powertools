#' Power calculation for one-way balanced analysis of variance contrast test
#'
#' @description
#' Performs sample size and power calculations for a test of a contrast in a one-way
#' ANOVA with balanced data (that is, equal sized groups). Can be used to solve for
#' power, n (sample size per group), or alpha. For unbalanced data, see
#' anova1way.c.unbal.
#'
#'
#' @param n The sample size per group.
#' @param mvec A vector of group means c(mu1, mu2, ...).
#' @param cvec A vector of contrast coefficients c(c1, c2, ...).
#' @param sd The estimated standard deviation within each group; defaults to 1.
#' @param Rsq The estimated R^2 for regressing the outcome on the covariates; defaults to 0.
#' @param ncov The number of covariates adjusted for in the model; defaults to 0.
#' @param alpha The significance level (type 1 error rate); defaults to 0.05.
#' @param power The specified level of power.
#' @param v Either TRUE for verbose output or FALSE (default) to output computed argument only.
#'
#' @return A list of the arguments (including the computed one).
#' @export
#'
#' @examples
#' anova1way.c.bal(n = 20, mvec = c(5, 10, 12), cvec = c(1, -1, 0), sd = 10, alpha = 0.025)
#' anova1way.c.bal(n = 20, mvec = c(5, 10, 12), cvec = c(1, 0, -1), sd = 10, alpha = 0.025)

anova1way.c.bal <- function (n = NULL, mvec = NULL, cvec = NULL, sd = 1,
                             Rsq = 0, ncov = 0, alpha = 0.05, power = NULL,
                             v = FALSE) {

  # Check if the arguments are specified correctly
  check.many(list(n, alpha, power), "oneof")
  check.param(n, "pos"); check.param(n, "min", min = 2)
  check.param(mvec, "req"); check.param(mvec, "vec")
  check.param(cvec, "req"); check.param(cvec, "vec")
  if (sum(cvec) != 0)
    stop("sum of contrast coefficients must equal 0")
  check.param(sd, "req"); check.param(sd, "pos")
  check.param(Rsq, "req"); check.param(Rsq, "uniti")
  check.param(ncov, "req"); check.param(ncov, "int")
  check.param(alpha, "unit")
  check.param(power, "unit")
  check.param(v, "req"); check.param(v, "bool")

  a <- length(mvec)
  if (a != length(cvec))
    stop("number of contrast coefficients must be equal to the number of groups")

  if (Rsq > 0 & ncov == 0)
    stop("please specify ncov or set Rsq to 0")

  # Calculate df and ncp
  p.body <- quote({
    lambda <- cvec %*% mvec / sd / sqrt(sum(cvec^2 / n)) / sqrt(1 - Rsq)
    df2 <- n * a - a - ncov
    stats::pf(q = stats::qf(alpha, 1, df2, lower.tail = FALSE),
              1, df2, lambda^2, lower.tail = FALSE)
  })

  # Use safe.uniroot function to calculate missing argument
  if (is.null(power)) {
    power <- eval(p.body)
    if (!v) return(power)
  }
  else if (is.null(n)) {
    n <- safe.uniroot(function(n) eval(p.body) - power, c(2 + ncov/a, 1e+05))$root
    if (!v) return(n)
  }
  else if (is.null(alpha)) {
    alpha <- safe.uniroot(function(alpha) eval(p.body) - power, c(1e-10, 1 - 1e-10))$root
    if (!v) return(alpha)
  }
  else stop("internal error", domain = NA)

  # Generate output text
  METHOD <- paste0("Balanced one-way analysis of ", ifelse(ncov < 1, "", "co"),
                   "variance\n     contrast test power calculation")
  out <- list(n = n, mvec = mvec, cvec = cvec, sd = sd,
              ncov = ncov, Rsq = Rsq, alpha = alpha, power = power,
              method = METHOD)

  # Print output as a power.htest object
  if (ncov < 1) out <- out[!names(out) %in% c("ncov", "Rsq")]
  structure(out, class = "power.htest")

}
