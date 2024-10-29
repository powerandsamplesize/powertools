#' Power calculation for test of simple effect for two-way unbalanced analysis of variance
#'
#' @description
#' Conducts power calculations for a test of a simple effect in a two-way
#' unbalanced (unequal cell sizes) ANOVA.
#' A "simple effect" is a contrast among the cell means.
#' This function does not solve for sample size.
#' For a test of a contrast in a balanced (equal
#' cell sizes) two-way ANOVA, anova2way.se.bal can also be used and can
#' solve for sample size. For a test of contrast among
#' factor levels, see anova2way.c.unbal.
#'
#'
#' @param nmatrix A matrix of sample sizes (see example).
#' @param mmatrix A matrix of group means (see example).
#' @param cmatrix A matrix of contrast coefficients (see example).
#' @param sd The estimated standard deviation within each group.
#' @param Rsq The estimated R^2 for regressing the outcome on the covariates; defaults to 0.
#' @param ncov The number of covariates adjusted for in the model; defaults to 0.
#' @param alpha The significance level (type 1 error rate); defaults to 0.05.
#' @param sides Either 1 or 2 (default) to specify a one- or two- sided hypothesis test.
#' @param v Either TRUE for verbose output or FALSE (default) to output computed argument only.
#'
#' @return A list of the arguments (including the computed power).
#' @export
#'
#' @examples
#' nmatrix <- matrix(c(30, 30, 30, 30, 30, 30), nrow = 2, byrow = TRUE)
#' mmatrix <- matrix(c(9.3, 8.9, 8.5, 8.7, 8.3, 7.3), nrow = 2, byrow = TRUE)
#' cmatrix <- matrix(c(-1, 0, 0, 1, 0, 0), nrow = 2, byrow = TRUE)
#' anova2way.se.unbal(nmatrix = nmatrix, mmatrix = mmatrix, cmatrix = cmatrix,
#' sd = 2, alpha = 0.025)

anova2way.se.unbal <- function (nmatrix = NULL, mmatrix = NULL, cmatrix = NULL,
                                sd = 0, Rsq = 0, ncov = 0, alpha = 0.05,
                                sides = 2, v = FALSE) {

  # Check if the arguments are specified correctly
  check.param(nmatrix, "req"); check.param(nmatrix, "mat")
  check.param(mmatrix, "req"); check.param(mmatrix, "mat")
  check.param(cmatrix, "req"); check.param(cmatrix, "mat")
  if (sum(cmatrix) != 0)
    stop("sum of contrast coefficients must equal 0")
  check.param(sd, "req"); check.param(sd, "pos")
  check.param(Rsq, "req"); check.param(Rsq, "uniti")
  check.param(ncov, "req"); check.param(ncov, "int")
  check.param(alpha, "req"); check.param(alpha, "unit")
  check.param(sides, "req"); check.param(sides, "vals", valslist = c(1, 2))
  check.param(v, "req"); check.param(v, "bool")

  a <- nrow(mmatrix)
  b <- ncol(mmatrix)
  if (a != nrow(cmatrix) | b != ncol(cmatrix))
    stop("number of contrast coefficients must be equal to the number of groups")
  if (a != nrow(nmatrix) | b != ncol(nmatrix))
    stop("number of sample sizes must be equal to the number of groups")

  if (any(nmatrix < 2))
    stop("number of observations in each group must be at least 2")

  if (Rsq > 0 & ncov == 0)
    stop("please specify ncov or set Rsq to 0")

  # See if there is an interaction
  fAB <- es.anova.f(means = mmatrix, sd = sd)$fAB
  intx <- ifelse(fAB == 0, FALSE, TRUE)

  # Get lambda
  lambda <- sum(cmatrix * mmatrix) / sd / sqrt(sum(cmatrix^2 / nmatrix)) /
    sqrt(1 - Rsq)

  # Calculate power
  N <- sum(nmatrix)
  df <- ifelse(intx, N - a * b - ncov, N - a - b + 1 - ncov)

  if (sides == 1)
    power <- stats::pt(q = stats::qt(alpha, df), df, lambda)
  else if (sides == 2)
    power <- stats::pf(q = stats::qf(alpha, 1, df), 1, df, lambda^2)

  if (!v) return(power)

  # Generate output text
  METHOD <- paste0("Unbalanced two-way analysis of ", ifelse(ncov < 1, "", "co"),
                   "variance\n     simple effects power calculation",
                   ifelse(intx, " with interaction", ""))
  out <- list(nmatrix = matrix.format(nmatrix),
              mmatrix = matrix.format(mmatrix),
              cmatrix = matrix.format(cmatrix),
              sd = sd, ncov = ncov, Rsq = Rsq, alpha = alpha, power = power,
              method = METHOD)

  # Print output as a power.htest object
  if (ncov < 1) out <- out[!names(out) %in% c("ncov", "Rsq")]
  structure(out, class = "power.htest")

}

