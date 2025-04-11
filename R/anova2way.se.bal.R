#' Power calculation for test of simple effect for two-way balanced analysis of variance
#'
#' @description
#' Conducts power and sample size calculations for a test of a simple effect in a two-way
#' balanced (equal cell sizes) ANOVA. A "simple effect" is a contrast among the cell means.
#' For a test of a contrast in an unbalanced (unequal
#' cell sizes) two-way ANOVA, see anova2way.se.unbal. For a test of contrast among
#' factor levels, see anova2way.c.bal.
#'
#'
#' @param n The sample size per cell.
#' @param mmatrix A matrix of cell means (see example).
#' @param cmatrix A matrix of contrast coefficients (see example).
#' @param sd The estimated standard deviation within each cell; defaults to 1.
#' @param Rsq The estimated R^2 for regressing the outcome on the covariates; defaults to 0.
#' @param ncov The number of covariates adjusted for in the model; defaults to 0.
#' @param alpha The significance level (type 1 error rate); defaults to 0.05.
#' @param power The specified level of power.
#' @param sides Either 1 or 2 (default) to specify a one- or two- sided hypothesis test.
#' @param v Either TRUE for verbose output or FALSE (default) to output computed argument only.
#'
#' @return A list of the arguments (including the computed one).
#' @export
#'
#' @examples
#' mmatrix <- matrix(c(9.3, 8.9, 8.5, 8.7, 8.3, 7.3), nrow = 2, byrow = TRUE)
#' cmatrix <- matrix(c(-1, 0, 0, 1, 0, 0), nrow = 2, byrow = TRUE)
#' anova2way.se.bal(n = 30, mmatrix = mmatrix, cmatrix = cmatrix, sd = 2, alpha = 0.025, sides = 1)

anova2way.se.bal <- function (n = NULL, mmatrix = NULL, cmatrix = NULL,
                              sd = 1, Rsq = 0, ncov = 0,
                              alpha = 0.05, power = NULL, sides = 2,
                              v = FALSE) {

  # Check if the arguments are specified correctly
  check.many(list(n, alpha, power), "oneof")
  check.param(n, "pos"); check.param(n, "min", min = 2)
  check.param(mmatrix, "req"); check.param(mmatrix, "mat")
  check.param(cmatrix, "req"); check.param(cmatrix, "mat")
  if (sum(cmatrix) != 0)
    stop("sum of contrast coefficients must equal 0")
  check.param(sd, "req"); check.param(sd, "pos")
  check.param(Rsq, "req"); check.param(Rsq, "uniti")
  check.param(ncov, "req"); check.param(ncov, "int")
  check.param(alpha, "unit")
  check.param(power, "unit")
  check.param(sides, "req"); check.param(sides, "vals", valslist = c(1, 2))
  check.param(v, "req"); check.param(v, "bool")

  a <- nrow(mmatrix)
  b <- ncol(mmatrix)
  if (a != nrow(cmatrix) | b != ncol(cmatrix))
    stop("number of contrast coefficients must be equal to the number of cells")

   if (Rsq > 0 & ncov == 0)
    stop("please specify ncov or set Rsq to 0")

  # See if there is an interaction
  fAB <- es.anova.f(means = mmatrix, sd = sd)$fAB
  intx <- ifelse(fAB == 0, FALSE, TRUE)

  # Get test statistic
  if (sides == 1)
    p.body <- quote({
      lambda <- sum(cmatrix * mmatrix) / sd / sqrt(sum(cmatrix^2 / n)) /
        sqrt(1 - Rsq)
      N <- a * b * n
      df <- ifelse(intx, N - a * b - ncov, N - a - b + 1 - ncov)
      stats::pt(q = stats::qt(alpha, df, lower.tail = FALSE), df, abs(lambda), lower.tail = FALSE)
    })
  else if (sides == 2)
    p.body <- quote({
      lambda <- sum(cmatrix * mmatrix) / sd / sqrt(sum(cmatrix^2 / n)) /
        sqrt(1 - Rsq)
      N <- a * b * n
      df2 <- ifelse(intx, N - a * b - ncov, N - a - b + 1 - ncov)
      stats::pf(q = stats::qf(alpha, 1, df2, lower.tail = FALSE), 1, df2, lambda^2, lower.tail = FALSE)
    })

  # Use safe.uniroot function to calculate missing argument
  if (is.null(power)) {
    power <- eval(p.body)
    if (!v) return(power)
  }
  else if (is.null(n)) {
    n <- safe.uniroot(function(n) eval(p.body) - power, c(2 + ncov/(a * b), 1e+05))$root
    if (!v) return(n)
  }
  else if (is.null(alpha)) {
    alpha <- safe.uniroot(function(alpha) eval(p.body) - power, c(1e-10, 1 - 1e-10))$root
    if (!v) return(alpha)
  }
  else stop("internal error", domain = NA)

  # Generate output text
  METHOD <- paste0("Balanced two-way analysis of ", ifelse(ncov < 1, "", "co"),
                   "variance\n     simple effects power calculation",
                   ifelse(intx, " with interaction", ""))
  out <- list(n = n, mmatrix = matrix.format(mmatrix),
              cmatrix = matrix.format(cmatrix),
              sd = sd, ncov = ncov, Rsq = Rsq, alpha = alpha, power = power,
              method = METHOD)

  # Print output as a power.htest object
  if (ncov < 1) out <- out[!names(out) %in% c("ncov", "Rsq")]
  structure(out, class = "power.htest")

}

