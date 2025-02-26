#' Power calculation for multiple alternative (at least one) primary continuous endpoints
#' assuming known covariance matrix
#'
#' @description
#' Calculates power and sample size for the case of comparing two groups on the means of K
#' continuous endpoints and concluding that the trial is a 'success' if the
#' null hypothesis is rejected for at least one of the K endpoints.
#' All mean differences must be specified as positive; the scale for some outcomes may need to be reversed
#' to meet this condition. All tests are assumed to be upper-tailed, one-sided tests. Can
#' solve for power, n1, n.ratio or alpha.
#'
#' To use a Bonferroni correction for multiple comparisons, specify alpha as
#' the desired familywise error rate (FWER) divided by K. For example, for one-sided FWER of 0.025
#' and K = 2 endpoints, specify alpha as 0.0125.
#'
#' Either sd and rho or Sigma must be specified.
#'
#' A known covariance matrix is assumed, which can result in a slight overestimate of power and
#' underestimate of required sample size.
#'
#'
#' @details
#' Sozu T, Sugimoto T, Hamasaki T, Evans SR (2015)
#' Sample Size Determination in Clinical Trials with Multiple Endpoints.
#' Springer International Publishing, Switzerland.
#'
#'
#' @param K The number of endpoints.
#' @param n1 The sample size for group 1.
#' @param n.ratio The ratio n2/n1 between the sample sizes of two groups; defaults to 1 (equal group sizes).
#' @param delta A vector of length K of the true mean differences mu1k - mu2k; must all be positive.
#' @param Sigma The covariance matrix of the K outcomes, of dimension K x K.
#' @param sd A vector of length K of the standard deviations of the K outcomes.
#' @param rho A vector of length 0.5*K*(K-1) of the correlations among the K outcomes.
#' @param alpha The significance level (type 1 error rate) for each test; defaults to 0.025.
#' A one-sided test is assumed.
#' @param power The specified level of power.
#' @param v Either TRUE for verbose output or FALSE (default) to output computed argument only.
#'
#' @return A list of the arguments (including the computed one).
#' @import mvtnorm
#' @export
#'
#' @examples
#' altprimary(K = 2, n1 = 100, delta = c(0.4, 0.5), sd = c(1, 1), rho = 0.3,
#' alpha = 0.025 / 2, power = NULL)
#'
#' Sigma <- matrix(c(1, 0.3, 0.3, 0.3, 1, 0.3, 0.3, 0.3, 1), nrow = 3, ncol = 3)
#' altprimary(K = 3, n1 = 100, delta = c(0.2, 0.2, 0.4), Sigma = Sigma,
#'    alpha = 0.025 / 3, power = NULL)

altprimary <- function(K, n1 = NULL, n.ratio = 1, delta = NULL, Sigma, sd, rho,
                       alpha = 0.025, power = NULL, v = FALSE){

  # Check if the arguments are specified correctly
  check.many(list(n1, n.ratio, power, alpha), "oneof")
  check.param(n1, "pos")
  check.param(power, "unit")
  check.param(alpha, "unit")
  check.param(K, "req"); check.param(K, "min", min = 1); check.param(K, "int")
  check.param(n.ratio, "pos")
  check.param(v, "req"); check.param(v, "bool")
  check.param(delta, "req"); check.param(delta, "vec")
  if(length(delta) != K)
    stop("length of 'delta' must be equal to 'K'")
  if(!all(delta > 0))
    stop("all effect sizes need to be positive")

  if(!missing(Sigma)){ # Sigma is given
    check.param(Sigma, "mat")
    if(nrow(Sigma) != ncol(Sigma))
      stop("covariance matrix 'Sigma' must be quadratic")
    if(nrow(Sigma) != K)
      stop("covariance matrix must have dimension 'K' x 'K'")
    if(max(abs(Sigma - t(Sigma))) > 1e-10)
      stop("matrix 'Sigma' must be symmetric")
  }
  if(missing(Sigma)){ # make Sigma with sd & rho
    check.param(sd, "vec")
    if(missing(sd) || missing(rho))
      stop("if 'Sigma' is missing 'sd' and 'rho' must be given.")
    if(length(sd) != K)
      stop("length of 'sd' must be equal to 'K'")
    if(!all(sd > 0))
      stop("all standard deviations need to be positive")
    if(length(rho) != 0.5*K*(K-1))
      stop("length of 'rho' must be equal to '0.5*K*(K-1)'")
    if(!all(rho >= -1 & rho < 1))
      stop("all correlations need to be between 0 and 1")
    Sigma <- matrix(0, nrow = K, ncol = K)
    iter <- 0
    for(i in 1:(K-1)){
      for(j in (i+1):K){
        iter <- iter + 1
        Sigma[i,j] <- rho[iter]*sd[i]*sd[j]
      }
    }
    Sigma <- Sigma + t(Sigma)
    diag(Sigma) <- sd^2
  }
  if(!all(eigen(Sigma)$values > 0))
    stop("matrix 'Sigma' must be positive definite")
  Sigma.cor <- stats::cov2cor(Sigma)

  ## calculations

  p.body <- quote({
    std.effect <- delta/sqrt(diag(Sigma))
    z.alpha <- stats::qnorm(1-alpha)
    crit.vals <- z.alpha - sqrt(n1*(n.ratio/(1+n.ratio)))*std.effect
    out <- 1 - mvtnorm::pmvnorm(lower = -crit.vals, sigma = Sigma.cor)
    out[1]
  })

  if(is.null(power)){
    power <- eval(p.body)
    if (!v) return(power)
  }
  else if (is.null(n1)) {
    n1 <- safe.uniroot(function(n1) eval(p.body) - power, c(2, 1e+07))$root
    if (!v) return(n1)
  }
  else if (is.null(n.ratio)) {
    n.ratio <- safe.uniroot(function(n.ratio) eval(p.body) - power,c(2/n1, 1e+07))$root
    if (!v) return(n.ratio)
  }
  else if (is.null(alpha)) {
    alpha <- safe.uniroot(function(alpha) eval(p.body) - power, c(1e-10, 1 - 1e-10))$root
    if (!v) return(alpha)
  }
  else stop("internal error")

  n <- c(n1, n1*n.ratio)
  METHOD <- "Power calculation for alternative (at least one) primary endpoints"
  structure(list(n = n, delta = delta, Sigma = matrix.format(Sigma),
                 alpha = alpha, power = power, method = METHOD),
            class = "power.htest")
}

