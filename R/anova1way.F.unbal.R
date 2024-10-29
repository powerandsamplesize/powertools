#' Power calculation for unbalanced one-way analysis of variance omnibus F test
#'
#' @description
#' Performs power calculation for an unbalanced (unequal group sizes)
#' one-way ANOVA omnibus F test, which tests for any differences among group means.
#' This function solves for power given other parameters. For balanced data
#' (equal-sized groups), anova1way.F.bal can be used and solves for more
#' parameters.
#'
#'
#'
#'
#' @param nvec A vector of group sample sizes c(n1, n2, ...).
#' @param mvec A vector of group mvec c(mu1, mu2, ...).
#' @param sd The estimated standard deviation within each group.
#' @param Rsq The estimated R^2 for regressing the outcome on the covariates; defaults to 0.
#' @param ncov The number of covariates adjusted for in the model; defaults to 0.
#' @param alpha The significance level or type 1 error rate; defaults to 0.05.
#' @param v Either TRUE for verbose output or FALSE to output computed argument only.
#'
#' @return A list of the arguments (including the computed power).
#' @export
#'
#' @examples
#' anova1way.F.unbal(nvec = c(10, 20, 30), mvec = c(5, 10, 12), sd = 10)

anova1way.F.unbal <- function (nvec = NULL, mvec = NULL, sd = NULL,
                                   Rsq = 0, ncov = 0, alpha = 0.05,
                                   v = FALSE) {

  # Check if the arguments are specified correctly
  check.param(nvec, "req"); check.param(nvec, "vec")
  check.param(mvec, "req"); check.param(mvec, "vec")
  check.param(sd, "req"); check.param(sd, "pos")
  check.param(Rsq, "req"); check.param(Rsq, "uniti")
  check.param(ncov, "req"); check.param(ncov, "int")
  check.param(alpha, "req"); check.param(alpha, "unit")
  check.param(v, "req"); check.param(v, "bool")

  a <- length(mvec)
  if (a != length(nvec))
    stop("number of sample sizes must equal to the number of groups")

  if (any(nvec < 2))
    stop("number of observations in each group must be at least 2")

  if (Rsq > 0 & ncov == 0)
    stop("please specify ncov or set Rsq to 0")

  # Get f effect size
  f <- es.anova.f(means = mvec, sd = sd, v = F)

  # Get marginal mean
  mvec <- matrix(mvec)
  mu <- mean(mvec)
  temp1 <- mvec - mu
  mmA <- rowMeans(temp1)

  # Get ncp
  N <- sum(nvec)
  props <- nvec / N
  ws <- props %*% mmA
  temp <- sapply(X = 1:a, FUN = function(i) nvec[i] * ((mmA[i] - ws) / sd)^2)
  Lambda <- sum(temp) / (1 - Rsq)

  # Calculate power
  df1 <- a - 1
  df2 <- N - a - ncov
  power <- stats::pf(stats::qf(alpha, df1, df2, lower.tail = FALSE),
                     df1, df2, Lambda, lower.tail = FALSE)
  if (!v) return(power)

  # Generate output text
  METHOD <- paste0("Unbalanced one-way analysis of ", ifelse(ncov < 1, "", "co"),
                   "variance\n     omnibus F test power calculation")
  out <- list(nvec = nvec, mvec = mvec, sd = sd,
              `f effect size` = f, ncov = ncov, Rsq = Rsq,
              alpha = alpha, power = power,
              method = METHOD)

  # Print output as a power.htest object
  if (ncov < 1) out <- out[!names(out) %in% c("ncov", "Rsq")]
  structure(out, class = "power.htest")

}

