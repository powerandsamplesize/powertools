#' Power calculations for multiple co-primary continuous endpoints assuming known covariance matrix
#'
#' @description
#' Computes power for test involving multiple co-primary continuous endpoints, assuming that the
#' covariance matrix (variances and covariances between endpoints) is known and therefore z-based test
#' statistics will be used. Studies with co-primary endpoints use “all-or-none” testing procedures and
#' only declare the trial to be a “success” if all endpoints are affirmed. All true mean differences
#' must be positive (the scale for some outcomes may need to be reversed
#' to meet this condition) and upper-tailed one-sided tests are assumed.
#' For the more realistic case that the covariance matrix is not known, see coprimary.t.
#'
#' Either sd and rho or Sigma must be specified.
#'
#' @details
#' See Crespi et al. (2025) for more details.
#' This function is based on the power.known.var function from the mpe R package and material from
#' Sozu T, Sugimoto T, Hamasaki T, Evans SR. (2015) Sample Size Determination in Clinical Trials
#' with Multiple Endpoints. Springer International Publishing, Switzerland.
#'
#'
#'
#' @param K The number of endpoints.
#' @param n1 The sample size for group 1.
#' @param n.ratio The ratio n2/n1 between the sample sizes of two groups; defaults to 1 (equal group sizes).
#' @param delta A vector of length K of the true mean differences mu1k - mu2k; must all be positive.
#' @param Sigma The covariance matrix of the K outcomes, of dimension K x K.
#' @param sd A vector of length K of the standard deviations of the K outcomes.
#' @param rho A vector of length 0.5*K*(K-1) of the correlations among the K outcomes.
#' @param alpha The significance level or type 1 error rate; defaults to 0.025. A one-sided test is assumed.
#' @param power The specified level of power.
#' @param v Either TRUE for verbose output or FALSE to output computed argument only.
#'
#' @return A list of the arguments (including the computed one).
#' @import mvtnorm
#' @export
#'
#' @examples
#' coprimary.z(K = 2, n1 = 100, delta = c(0.4, 0.5), sd = c(1, 1), rho = 0.3,
#' alpha = 0.025, power = NULL)
#'
#' Sigma <- matrix(c(1, 0.3, 0.3, 0.3, 1, 0.3, 0.3, 0.3, 1), nrow = 3, ncol = 3)
#' coprimary.z(K = 3, n1 = NULL, delta = c(0.2, 0.3, 0.4), Sigma = Sigma, alpha = 0.025, power = 0.8)
#'

coprimary.z <- function(K, n1 = NULL, n.ratio = 1, delta = NULL, Sigma, sd, rho,
                        alpha = 0.025, power = NULL, v = FALSE) {

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

  if(!missing(Sigma)){
    check.param(Sigma, "mat")
    if(nrow(Sigma) != ncol(Sigma))
      stop("covariance matrix 'Sigma' must be square")
    if(nrow(Sigma) != K)
      stop("covariance matrix must have dimension 'K' x 'K'")
    if(max(abs(Sigma - t(Sigma))) > 1e-10)
      stop("matrix 'Sigma' must be symmetric")
  }
  if(missing(Sigma)){
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

  p.body <- quote({
    std.effect <- delta/sqrt(diag(Sigma))
    z.alpha <- stats::qnorm(1-alpha)
    crit.vals <- z.alpha - sqrt(n1*(n.ratio/(1+n.ratio)))*std.effect
    out <- mvtnorm::pmvnorm(upper = -crit.vals, sigma = Sigma.cor)
    out[1]
  })

  ## calculations
  if(is.null(power)){
    power <- eval(p.body)
    if (!v) return(power)
  }
  else if (is.null(n1)) {
    n1 <- stats::uniroot(function(n1) eval(p.body) - power, c(2, 1e+07))$root
    if (!v) return(n1)
  }
  else if (is.null(n.ratio)) {
    n.ratio <- stats::uniroot(function(n.ratio) eval(p.body) - power,c(2/n1, 1e+07))$root
    if (!v) return(n.ratio)
  }
  else if (is.null(alpha)) {
    alpha <- stats::uniroot(function(alpha) eval(p.body) - power, c(1e-10, 1 - 1e-10))$root
    if (!v) return(alpha)
  }
  else stop("internal error")

  n <- c(n1, n1*n.ratio)
  METHOD <- "Power calculation for multiple co-primary endpoints (covariance matrix known)"
  structure(list(n = n, delta = delta, Sigma = matrix.format(Sigma),
                 alpha = alpha, power = power, method = METHOD),
            class = "power.htest")

}
