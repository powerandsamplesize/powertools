#' Cohen f effect size calculation for one- or two- way analysis of variance
#'
#' @description
#' Calculates teh Cohen f effect size for a one- or two-way ANOVA.
#' Takes as input the cell or group means for a one- or two-way ANOVA and the common
#' standard deviation and outputs the f effect size, as defined by Cohen (1988). Note that
#' this effect size calculation is only valid when cell/group sizes are equal.
#'
#' @details
#' Cohen J (1988) Statistical Power Analysis for the Behavioral Sciences, 2nd edition.
#' Lawrence Erlbaum Associates, Hillsdale, New Jersey
#'
#'
#' @param means A vector or matrix of group means.
#' @param sd The estimated standard deviation within each group.
#' @param v Either TRUE for verbose output or FALSE to output computed argument only.
#'
#' @return Various calculated f effect sizes.
#' @export
#'
#' @examples
#' es.anova.f(means = c(5, 10, 12), sd = 10)
#' mmatrix <- matrix(c(9.3, 8.9, 8.5, 8.7, 8.3, 7.9), nrow = 2, byrow = TRUE)
#' es.anova.f(means = mmatrix, sd = 2)

es.anova.f <- function (means = NULL, sd = NULL, v = TRUE) {

  # Check if the arguments are specified correctly
  check.param(means, "req")
  check.param(sd, "req"); check.param(sd, "pos")
  check.param(v, "req"); check.param(v, "bool")

  if (is.vector(means)) {
    a <- length(means)
    b <- 1
    check.param(means, "vec")
    mmatrix <- matrix(means)
  }
  else if (is.matrix(means)) {
    a <- nrow(means)
    b <- ncol(means)
    check.param(means, "mat")
    mmatrix <- means
  }

  # Get grand mean, marginal means, and interaction effects
  mu <- mean(mmatrix)
  temp1 <- mmatrix - mu
  mmA <- rowMeans(temp1)
  mmB <- colMeans(temp1)
  temp2 <- sweep(x = temp1, MARGIN = 2, STATS = mmB, FUN = "-")
  ints <- sweep(x = temp2, MARGIN = 1, STATS = mmA, FUN = "-")

  # Get sds and effect sizes
  sdA <- sqrt(sum(mmA^2) / a)
  sdB <- sqrt(sum(mmB^2) / b)
  sdAB <- sqrt(sum(ints^2) / (a * b))
  fA <- sdA / sd
  fB <- sdB / sd
  fAB <- sdAB / sd
  if (fAB < 0.000001) fAB <- 0

  # Print output as a power.htest object
  if (is.vector(means)) {
    if (!v) return(fA)
    METHOD <- "Cohen's f effect size calculation for\n     one-way analysis of variance"
    structure(list(fA = fA, method = METHOD), class = "power.htest")

  } else if (is.matrix(means)) {
    METHOD <- "Cohen's f effect size calculation for\n     two-way analysis of variance"
    mrows <- c()
    for (i in 1:a) mrows <- c(mrows, paste(means[i,], collapse = ', '))
    means <- paste(mrows, collapse = "\n                  ")
    structure(list(fA = fA, fB = fB, fAB = fAB,
                   method = METHOD), class = "power.htest")
  }

}
