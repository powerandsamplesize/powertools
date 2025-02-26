#' Power for detecting treatment effect heterogeneity in a cluster randomized trial with a continuous outcome
#'
#' @description
#' This function performs power and sample size calculations for detecting a treatment-by-covariate
#' interaction effect in a two-arm cluster randomized trial
#' with a continuous outcome when the data will be analyzed using a linear mixed effects model
#' (random intercept for cluster and fixed effect for the treatment-by-covariate interaction).
#' Can solve for power, beta, J1, J.ratio or m.
#'
#'
#' @details
#' This function is based on Yang et al (2020). If the covariate is a cluster-level covariate,
#' then icc.x should be set to 1 (the covariate does not vary within cluster).
#'
#' Yang S, Li F, Starks MA, Hernandez AF, Mentz RJ, Choudhury KR (2020) Sample size requirements for detecting
#' treatment effect heterogeneity in cluster randomized trials. Statistics in Medicine 39:4218-4237.
#'
#'
#'
#' @param m The number of subjects per cluster.
#' @param J1 The number of clusters in arm 1.
#' @param J.ratio The ratio J2/J1 between the number of clusters in the two arms; defaults to 1 (equal clusters per arm).
#' @param beta The regression coefficient for the treatment-by-covariate interaction term.
#' @param sd.x The standard deviation of the covariate.
#' @param sd.yx The standard deviation of the outcome variable adjusting for the covariate.
#' @param icc.x The intraclass correlation coefficient for the covariate; defaults to 0.
#' @param icc.yx The intraclass correlation coefficient for the outcome adjusting for the covariate; defaults to 0.
#' @param alpha The significance level (type 1 error rate); defaults to 0.05.
#' @param power The specified level of power.
#' @param sides Either 1 or 2 (default) to specify a one- or two- sided hypothesis test.
#' @param v Either TRUE for verbose output or FALSE (default) to output computed argument only.
#'
#' @return A list of the arguments (including the computed one).
#' @export
#'
#' @examples
#' crt.parallel.hte(beta = 1, m = 27, J1 = 20, sd.x = 12.7, sd.yx = 71, icc.x = 0.08, icc.yx = 0.04)


crt.parallel.hte <- function (m = NULL, J1 = NULL, J.ratio = 1, beta = NULL,
                              sd.x = NULL, sd.yx = NULL,
                              icc.x = 0, icc.yx = 0,
                              alpha = 0.05, power = NULL, sides = 2,
                              v = FALSE) {

  # Check if the arguments are specified correctly
  check.many(list(m, J1, J.ratio, beta, alpha, power), "oneof")
  check.param(m, "pos")
  check.param(J.ratio, "pos")
  check.param(J1, "min", min = 2)
  if (!is.null(J1) & !is.null(J.ratio))
    check.param(J1 * J.ratio, "min", min = 2)
  check.param(beta, "num")
  check.param(sd.x, "req"); check.param(sd.x, "pos")
  check.param(sd.yx, "req"); check.param(sd.yx, "pos")
  check.param(icc.x, "req"); check.param(icc.x, "unitii")
  check.param(icc.yx, "req"); check.param(icc.yx, "uniti")
  check.param(alpha, "unit")
  check.param(power, "unit")
  check.param(sides, "req"); check.param(sides, "vals", valslist = c(1, 2))
  check.param(v, "req"); check.param(v, "bool")

  # Calculate power
    p.body <- quote({
      sd.w <- sqrt(J.ratio) / (1 + J.ratio)
      J <- J1 + J1 * J.ratio
      A <- (abs(beta) * sd.w * sd.x)
      B <- sd.yx * sqrt(((1 - icc.yx) * (1 + (m-1) * icc.yx)) / (J * m * (1 + (m-2) * icc.yx - (m-1) * icc.x * icc.yx)))
      stats::pnorm(stats::qnorm(alpha / sides) + A / B)
    })

  # Use uniroot to calculate missing argument
  if (is.null(alpha)) {
    alpha <- safe.uniroot(function(alpha) eval(p.body) - power, c(1e-10, 1 - 1e-10))$root
    if (!v) return(alpha)
  }
  else if (is.null(power)) {
    power <- eval(p.body)
    if (!v) return(power)
  }
  else if (is.null(J1)) {
    J1 <- safe.uniroot(function(J1) eval(p.body) - power, c(2 + 1e-10, 1e+07))$root
    if (!v) return(J1)
  }
  else if (is.null(J.ratio)) {
    J.ratio <- safe.uniroot(function(J.ratio) eval(p.body) - power, c(1e-10, 1e+07))$root
    if (!v) return(J.ratio)
  }
  else if (is.null(m)) {
    m <- safe.uniroot(function(m) eval(p.body) - power, c(2, 1e+07))$root
    if (!v) return(m)
  }
  else if (is.null(beta)) {
    beta <- safe.uniroot(function(beta) eval(p.body) - power, c(1e-07, 1e+07))$root
    if (!v) return(beta)
  }

  # Generate output text
  METHOD <- "Power for treatment-by-covariate interaction in a cluster randomized trial with a continuous outcome"
  J <- c(J1, J1 * J.ratio)
  sd <- c(sd.x, sd.yx)
  icc <- c(icc.x, icc.yx)
  out <- list(m = m, `J1, J2` = J, beta = beta, `sd.x, sd.yx` = sd, `icc.x, icc.yx` = icc,
              alpha = alpha, power = power, sides = sides, method = METHOD)

  # Print output as a power.htest object

  structure(out, class = "power.htest")

}

