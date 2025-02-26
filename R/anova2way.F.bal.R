#' Power calculation for two-way balanced analysis of variance F tests
#'
#' @description
#' Performs sample size and power calculations for F tests in a two-way
#' ANOVA with balanced data (that is, equal cell sizes). For a given
#' matrix of cell means, computes power or required cell size
#' for each factor and for their interaction, if an interaction is present.
#' For unbalanced data (unequal cell sizes),
#' see anova2way.F.unbal.
#'
#'
#' @param n The sample size per cell
#' @param mmatrix A matrix of cell means (see example).
#' @param sd The estimated standard deviation within each cell; defaults to 1.
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
#' mmatrix <- matrix(c(9.3, 8.9, 8.5, 8.7, 8.3, 7.9), nrow = 2, byrow = TRUE)
#' anova2way.F.bal(n = 30, mmatrix = mmatrix, sd = 2, alpha = 0.05)
#' mmatrix <- matrix(c(9.3, 8.9, 8.5, 8.7, 8.3, 7.3), nrow = 2, byrow = TRUE)
#' anova2way.F.bal(n = 30, mmatrix = mmatrix, sd = 2, alpha = 0.05)
#' mmatrix <- matrix(c(9.3, 8.9, 8.5, 8.7, 8.3, 7.9), nrow = 2, byrow = TRUE)
#' anova2way.F.bal(n = 30, mmatrix = mmatrix, sd = 2, Rsq = 0.4, ncov = 1, alpha = 0.05)

anova2way.F.bal <- function (n = NULL, mmatrix = NULL, sd = 1,
                             Rsq = 0, ncov = 0, alpha = 0.05, power = NULL,
                             v = FALSE) {

  # Check if the arguments are specified correctly
  check.many(list(n, alpha, power), "oneof")
  check.param(n, "pos"); check.param(n, "min", min = 2)
  check.param(mmatrix, "req"); check.param(mmatrix, "mat")
  check.param(sd, "req"); check.param(sd, "pos")
  check.param(Rsq, "req"); check.param(Rsq, "uniti")
  check.param(ncov, "req"); check.param(ncov, "int")
  check.param(alpha, "unit")
  check.param(power, "unit")
  check.param(v, "req"); check.param(v, "bool")

  a <- nrow(mmatrix)
  b <- ncol(mmatrix)

  if (Rsq > 0 & ncov == 0)
    stop("please specify ncov or set Rsq to 0")

  # Set default values if given
  nA <- n; nB <- n; nAB <- n
  powerA <- power; powerB <- power; powerAB <- power

  # Get f effect sizes
  es <- es.anova.f(means = mmatrix, sd = sd)
  fA <- es$fA
  fB <- es$fB
  fAB <- es$fAB
  intx <- ifelse(fAB == 0, FALSE, TRUE)

  # Calculate df's and ncp's
  p.body.A <- quote({
    N <- n * a * b
    LambdaA <- N * fA^2 / (1 - Rsq)
    df1 <- a - 1
    df2 <- ifelse(intx, N - a * b - ncov, N - a - b + 1 - ncov)
    stats::pf(q = stats::qf(alpha, df1, df2, lower.tail = FALSE),
              df1, df2, LambdaA, lower.tail = FALSE)
  })
  p.body.B <- quote({
    N <- n * a * b
    LambdaB <- N * fB^2 / (1 - Rsq)
    df1 <- b - 1
    df2 <- ifelse(intx, N - a * b - ncov, N - a - b + 1 - ncov)
    stats::pf(q = stats::qf(alpha, df1, df2, lower.tail = FALSE),
              df1, df2, LambdaB, lower.tail = FALSE)
  })
  p.body.AB <- quote({
    N <- n * a * b
    LambdaAB <- N * fAB^2 / (1 - Rsq)
    df1 <- (a - 1) * (b - 1)
    df2 <- N - a * b - ncov
    stats::pf(q = stats::qf(alpha, df1, df2, lower.tail = FALSE),
              df1, df2, LambdaAB, lower.tail = FALSE)
  })

  NOTE <- "The 3rd value for f and power or n is for the interaction"
  if(!v & intx) print(paste("NOTE:", NOTE))

  # Use safe.uniroot function to calculate missing argument
  if (is.null(power)) {
    powerA <- round(eval(p.body.A), 4)
    powerB <- round(eval(p.body.B), 4)
    if (intx) powerAB <- round(eval(p.body.AB), 4)
    if (!v) return(c(powerA = powerA, powerB = powerB, powerAB = powerAB))
  }
  else if (is.null(n)){
    nA <- round(safe.uniroot(function(n) eval(p.body.A) - power, c(2 + ncov/(a * b), 1e+05))$root, 4)
    nB <- round(safe.uniroot(function(n) eval(p.body.B) - power, c(2 + ncov/(a * b), 1e+05))$root, 4)
    if (intx)
      nAB <- round(safe.uniroot(function(n) eval(p.body.AB) - power, c(2 + ncov/(a * b), 1e+05))$root, 4)
    if (!v) return(c(nA = nA, nB = nB, nAB = nAB))
  }
  else if (is.null(alpha)) {
    alpha <- safe.uniroot(function(alpha) eval(p.body.A) - power, c(1e-10, 1 - 1e-10))$root
    if (!v) return(alpha)
  }
  else stop("internal error", domain = NA)

  # Generate output text
  if (is.null(power))
    if (intx) power <- c(powerA, powerB, powerAB) else power <- c(powerA, powerB)
  else if (is.null(n))
    if (intx) n <- c(nA, nB, nAB) else n <- c(nA, nB)
  if (intx) f <- c(round(fA, 4), round(fB, 4), round(fAB, 4))
  else f <- c(round(fA, 4), round(fB, 4))
  METHOD <- paste0("Balanced two-way analysis of ", ifelse(ncov < 1, "", "co"),
                   "variance\n     omnibus F test power calculation",
                   ifelse(intx, " with interaction", ""))
  out <- list(n = n, mmatrix = matrix.format(mmatrix),
              sd = sd, `f effect size` = f, ncov = ncov, Rsq = Rsq,
              alpha = alpha, power = power,
              method = METHOD, note = NOTE)

  # Print output as a power.htest object
  if (ncov < 1) out <- out[!names(out) %in% c("ncov", "Rsq")]
  if (!intx) out <- out[!names(out) == "note"]
  structure(out, class = "power.htest")

}

