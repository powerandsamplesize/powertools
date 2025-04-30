#' Power calculation for two-way unbalanced analysis of variance contrast test
#'
#' @description
#' Calculates power for a test of a contrast between levels of a factor in a two-way
#' ANOVA with unbalanced data (that is, unequal cell sizes).
#' This function
#' only solves for power. For a two-way balanced ANOVA,
#' (equal cell sizes), anova2way.c.bal can also be used, and will solve for
#' quantities other than power. For a test of a contrast among cell means, see
#' anova2way.se.unbal.
#'
#'
#' @param nmatrix A matrix of cell sample sizes (see example).
#' @param mmatrix A matrix of cell means (see example).
#' @param cvec A vector of contrast coefficients c(c1, c2, ...).
#' @param factor Either "a" (rows) or "b" (columns) depending on which factor the contrast test is being made on.
#' @param sd The estimated standard deviation within each cell.
#' @param Rsq The estimated R^2 for regressing the outcome on the covariates; defaults to 0.
#' @param ncov The number of covariates adjusted for in the model; defaults to 0.
#' @param alpha The significance level (type 1 error rate); defaults to 0.05.
#' @param v Either TRUE for verbose output or FALSE (default) to output computed argument only.
#'
#' @return A list of the arguments (including the computed one).
#' @export
#'
#' @examples
#' nmatrix <- matrix(c(30, 30, 30, 30, 30, 30), nrow = 2, byrow = TRUE)
#' mmatrix <- matrix(c(9.3, 8.9, 8.5, 8.7, 8.3, 7.3), nrow = 2, byrow = TRUE)
#' anova2way.c.unbal(nmatrix = nmatrix, mmatrix = mmatrix, cvec = c(1, 0, -1),
#' factor = "b", sd = 2, alpha = 0.05)

anova2way.c.unbal <- function (nmatrix = nmatrix, mmatrix = NULL, cvec = NULL,
                               factor = c("a", "b"), sd = NULL, Rsq = 0, ncov = 0,
                               alpha = 0.05, v = FALSE) {

  # Check if the arguments are specified correctly
  check.param(nmatrix, "req"); check.param(nmatrix, "mat")
  check.param(mmatrix, "req"); check.param(mmatrix, "mat")
  check.param(cvec, "req"); check.param(cvec, "vec"); check.param(cvec, "sum0")
  check.param(factor, "req"); check.param(factor, "vals", valslist = c("a", "b"))
  check.param(sd, "req"); check.param(sd, "pos")
  check.param(Rsq, "req"); check.param(Rsq, "uniti")
  check.param(ncov, "req"); check.param(ncov, "int")
  check.param(alpha, "req"); check.param(alpha, "unit")
  check.param(v, "req"); check.param(v, "bool")

  a <- nrow(mmatrix)
  b <- ncol(mmatrix)
  factor <- match.arg(factor)
  if (switch(factor, "a" = a, "b" = b) != length(cvec))
    stop("number of contrast coefficients must be equal to the number of factor levels")

  if (any(nmatrix < 2))
    stop("number of observations in each cell must be at least 2")

  if (Rsq > 0 & ncov == 0)
    stop("please specify ncov or set Rsq to 0")

  # See if there is an interaction
  fAB <- es.anova.f(means = mmatrix, sd = sd)$fAB
  intx <- ifelse(fAB == 0, FALSE, TRUE)

  # Get grand mean and marginal means
  mu <- mean(mmatrix)
  temp1 <- mmatrix - mu
  mmA <- rowMeans(temp1)
  mmB <- colMeans(temp1)

  # Get lambda (Lambda = lambda^2)
  num <- switch(factor, "a" = mmA, "b" = mmB) %*% cvec
  temp <- switch(factor,
                 "a" = sapply(X = 1:a, FUN = function(i)
                       cvec[i]^2 / sum(nmatrix[i, ])),
                 "b" = sapply(X = 1:b, FUN = function(j)
                       cvec[j]^2 / sum(nmatrix[, j])))
  den <- sd * sqrt(sum(temp))
  lambda <- num / den / sqrt(1 - Rsq)

  # Calculate power
  N <- sum(nmatrix)
  df2 <- ifelse(intx, N - a * b - ncov, N - a - b + 1 - ncov)
  power <- stats::pf(q = stats::qf(alpha, 1, df2, lower.tail = FALSE),
                     1, df2, lambda^2, lower.tail = FALSE)
  if (!v) return(power)

  # Generate output text
  METHOD <- paste0("Unalanced two-way analysis of ", ifelse(ncov < 1, "", "co"),
                   "variance\n     contrast test power calculation",
                   ifelse(intx, " with interaction", ""))
  out <- list(nmatrix = matrix.format(nmatrix),
              mmatrix = matrix.format(mmatrix),
              factor = factor,
              cvec = cvec, sd = sd, ncov = ncov, Rsq = Rsq,
              alpha = alpha, power = power,
              method = METHOD)

  # Print output as a power.htest object
  if (ncov < 1) out <- out[!names(out) %in% c("ncov", "Rsq")]
  structure(out, class = "power.htest")

}

