#' Power calculation for balanced one-way ANOVA omnibus F test
#'
#' @description
#' This function performs power and sample size calculations for the overall (omnibus) F test
#' in a balanced (equal-sized groups) one-way analysis of variance (ANOVA). Can be used to solve for
#' power, n (sample size per group), or alpha.
#' For an unbalanced one-way
#' ANOVA F test (that is, unequal group sample sizes), use 'anova1way.F.unbal'.
#' For contrast tests in a one-way ANOVA, see 'anova1way.c.bal' and 'anova1way.c.unbal'.
#'
#'
#' @param n The sample size per group.
#' @param mvec A vector of group means c(mu1, mu2, ...).
#' @param sd The estimated standard deviation within each group; defaults to 1.
#' @param Rsq The estimated R^2 for regressing the outcome on the covariates; defaults to 0.
#' @param ncov The number of covariates adjusted for in the model; defaults to 0.
#' @param alpha The significance level or type 1 error rate; defaults to 0.05.
#' @param power The specified level of power.
#' @param v Either TRUE for verbose output or FALSE to output computed argument only.
#'
#' @return A list of the arguments (including the computed one).
#' @export
#'
#' @examples
#' anova1way.F.bal(n = 20, mvec = c(5, 10, 12), sd = 10)
#' anova1way.F.bal(n = NULL, mvec = c(-0.25, 0.25), sd = 1, Rsq = 0.5^2, ncov = 1, power = 0.8)

anova1way.F.bal <- function (n = NULL, mvec = NULL, sd = 1, Rsq = 0, ncov = 0,
                             alpha = 0.05, power = NULL, v = FALSE) {

  # Check if the arguments are specified correctly
  check.many(list(n, alpha, power), "oneof")
  check.param(n, "pos"); check.param(n, "min", min = 2)
  check.param(mvec, "req"); check.param(mvec, "vec")
  check.param(sd, "req"); check.param(sd, "pos")
  check.param(Rsq, "req"); check.param(Rsq, "uniti")
  check.param(ncov, "req"); check.param(ncov, "int")
  check.param(alpha, "unit")
  check.param(power, "unit")
  check.param(v, "req"); check.param(v, "bool")

  a <- length(mvec)

  if (Rsq > 0 & ncov == 0)
    stop("please specify ncov or set Rsq to 0")

  # Get f effect size
  f <- es.anova.f(means = mvec, sd = sd, v = F)

  # Calculate df's and ncp
  p.body <- quote({
    Lambda <- a * n * f^2 / (1 - Rsq)
    df1 <- a - 1
    df2 <- (n - 1) * a - ncov
    stats::pf(stats::qf(alpha, df1, df2, lower.tail = FALSE),
              df1, df2, Lambda, lower.tail = FALSE)
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
                   "variance\n     omnibus F test power calculation")
  out <- list(n = n, mvec = mvec, sd = sd,
              `f effect size` = f, ncov = ncov, Rsq = Rsq,
              alpha = alpha, power = power,
              method = METHOD)

  # Print output as a power.htest object
  if (ncov < 1) out <- out[!names(out) %in% c("ncov", "Rsq")]
  structure(out, class = "power.htest")

}
