#' Power calculation for two-way unbalanced analysis of variance F tests
#'
#' @description
#' Performs sample size and power calculations for F tests in a two-way
#' ANOVA with unbalanced data (that is, unequal sized cells). For given
#' matrix of cell means and matrix of cell sample sizes, computes power
#' for each factor and for their interaction, if an interaction is present.
#' This function does not solve for cell sizes.
#' For balanced data (equal cell sizes),
#' see anova2way.F.unbal, which can solve for cell size.
#'
#'
#' @param nmatrix A matrix of cell sample sizes (see example).
#' @param mmatrix A matrix of cell means (see example).
#' @param sd The estimated standard deviation within each cell
#' @param Rsq The estimated R^2 for regressing the outcome on the covariates; defaults to 0.
#' @param ncov The number of covariates adjusted for in the model; defaults to 0.
#' @param alpha The significance level (type 1 error rate); defaults to 0.05.
#' @param v Either TRUE for verbose output or FALSE (default) to output computed argument only.
#'
#' @return A list of the arguments (including the computed power).
#' @export
#'
#' @examples
#' nmatrix <- matrix(c(30, 30, 30, 30, 30, 30), nrow = 2, byrow = TRUE)
#' mmatrix <- matrix(c(9.3, 8.9, 8.5, 8.7, 8.3, 7.9), nrow = 2, byrow = TRUE)
#' anova2way.F.unbal(nmatrix = nmatrix, mmatrix = mmatrix, sd = 2, alpha = 0.05)
#' nmatrix <- matrix(c(30, 30, 30, 30, 30, 30), nrow = 2, byrow = TRUE)
#' mmatrix <- matrix(c(9.3, 8.9, 8.5, 8.7, 8.3, 7.3), nrow = 2, byrow = TRUE)
#' anova2way.F.unbal(nmatrix = nmatrix, mmatrix = mmatrix, sd = 2, alpha = 0.05)
#' nmatrix <- matrix(c(30, 30, 30, 30, 30, 30), nrow = 2, byrow = TRUE)
#' mmatrix <- matrix(c(9.3, 8.9, 8.5, 8.7, 8.3, 7.9), nrow = 2, byrow = TRUE)
#' anova2way.F.unbal(nmatrix = nmatrix, mmatrix = mmatrix, sd = 2, Rsq = 0.4^2,
#' ncov = 1, alpha = 0.05)

anova2way.F.unbal <- function (nmatrix = NULL, mmatrix = NULL, sd = NULL,
                               Rsq = 0, ncov = 0, alpha = 0.05, v = FALSE) {

  # Check if the arguments are specified correctly
  check.param(nmatrix, "req"); check.param(nmatrix, "mat")
  check.param(mmatrix, "req"); check.param(mmatrix, "mat")
  check.param(sd, "req"); check.param(sd, "pos")
  check.param(Rsq, "req"); check.param(Rsq, "uniti")
  check.param(ncov, "req"); check.param(ncov, "int")
  check.param(alpha, "req"); check.param(alpha, "unit")
  check.param(v, "req"); check.param(v, "bool")

  a <- nrow(mmatrix)
  b <- ncol(mmatrix)
  if(nrow(nmatrix) != a | ncol(nmatrix) != b)
    stop("number of sample sizes must equal to the number of cells")

  if (any(nmatrix < 2))
    stop("number of observations in each cell must be at least 2")

  if (Rsq > 0 & ncov == 0)
    stop("please specify ncov or set Rsq to 0")

  # Get f effect sizes
  es <- es.anova.f(means = mmatrix, sd = sd)
  fA <- es$fA
  fB <- es$fB
  fAB <- es$fAB
  intx <- ifelse(fAB == 0, FALSE, TRUE)

  # Get marginal means
  mu <- mean(mmatrix)
  temp1 <- mmatrix - mu
  mmA <- rowMeans(temp1)
  mmB <- colMeans(temp1)
  temp2 <- sweep(x = temp1, MARGIN = 2, STATS = mmB, FUN = "-")
  ints <- sweep(x = temp2, MARGIN = 1, STATS = mmA, FUN = "-")

  # Get Lambdas
  LambdaA <- 0
  for (i in 1:a) {
    for (j in 1:b) {
      temp <- mmA[i] / (sd / sqrt(nmatrix[i, j]))
      LambdaA <- LambdaA + temp^2
    }
  }
  LambdaA <- LambdaA / (1 - Rsq)
  LambdaB <- 0
  for (j in 1:b) {
    for (i in 1:a) {
      temp <- mmB[j] / (sd / sqrt(nmatrix[i, j]))
      LambdaB <- LambdaB + temp^2
    }
  }
  LambdaB <- LambdaB / (1 - Rsq)
  LambdaAB <- 0
  for (i in 1:a) {
    for (j in 1:b) {
      temp <- ints[i, j] / (sd / sqrt(nmatrix[i, j]))
      LambdaAB <- LambdaAB + temp^2
    }
  }
  LambdaAB <- LambdaAB / (1 - Rsq)

  NOTE <- "The 3rd value for f and power or n is for the interaction"
  if(!v & intx) print(paste("NOTE:", NOTE))

  # Calculate power
  N <- sum(nmatrix)
  df1A <- a - 1; df1B <- b - 1; df1AB <- (a - 1) * (b - 1)
  df2 <- ifelse(intx, N - a * b - ncov, N - a - b + 1 - ncov)
  powerA <- round(stats::pf(stats::qf(alpha, df1A, df2, lower.tail = FALSE),
                  df1A, df2, LambdaA, lower.tail = FALSE), 4)
  powerB <- round(stats::pf(stats::qf(alpha, df1B, df2, lower.tail = FALSE),
                  df1B, df2, LambdaB, lower.tail = FALSE), 4)
  if (intx) {
    powerAB <- round(stats::pf(stats::qf(alpha, df1AB, df2, lower.tail = FALSE),
                     df1AB, df2, LambdaAB, lower.tail = FALSE), 4)
    if (!v) return(c(powerA = powerA, powerB = powerB, powerAB = powerAB))
  }
  else {
    powerAB <- 0
    if (!v) return(c(powerA = powerA, powerB = powerB))
  }


  # Generate output text
  if (intx) power <- c(powerA, powerB, powerAB) else power <- c(powerA, powerB)
  if (intx) f <- c(round(fA, 4), round(fB, 4), round(fAB, 4))
  else f <- c(round(fA, 4), round(fB, 4))
  METHOD <- paste0("Unalanced two-way analysis of ", ifelse(ncov < 1, "", "co"),
                   "variance\n     omnibus F test power calculation")
  out <- list(nmatrix = matrix.format(nmatrix),
              mmatrix = matrix.format(mmatrix),
              sd = sd, `f effect size` = f, ncov = ncov, Rsq = Rsq,
              alpha = alpha, power = power,
              method = METHOD, note = NOTE)

  # Print output as a power.htest object
  if (ncov < 1) out <- out[!names(out) %in% c("ncov", "Rsq")]
  if (!intx) out <- out[!names(out) == "note"]
  structure(out, class = "power.htest")

}

