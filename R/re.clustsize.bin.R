#' Relative efficiency of a cluster randomized trial with binary outcome
#' with varying cluster sizes
#'
#' @description
#' For a binary outcome, computes the relative efficiency (ratio of the variances) of a cluster randomized trial
#' with varying cluster sizes to that of a cluster randomized trial with constant cluster sizes,
#' assuming equal total number of subjects.
#'
#' @details
#' Candel MJJM and van Breukelen GJP (2010) Sample size adjustments for varying cluster sizes in
#' cluster randomized trials with binary outcomes analyzed with second-order PQL mixed logistic regression.
#' Statistics in Medicine 29(14):1488-1501.
#'
#'
#'
#' @param m The number of subjects per cluster or the mean cluster size (if unequal number of participants per cluster).
#' @param m.sd The standard deviation of cluster sizes (in case of unequal number of participants per cluster).
#' @param pc The probability of the outcome in control clusters.
#' @param pt The probability of the outcome in treatment clusters.
#' @param sigma.u Standard deviation of the cluster random effect.
#' @return The computed RE.
#' @export
#'
#' @examples
#' re.clustsize.bin(m = 60, m.sd = 45, pc = 0.25, pt = 0.15, sigma.u = 0.3)

re.clustsize.bin <- function(m, m.sd, pc, pt, sigma.u){
  check.param(m, "req"); check.param(m, "pos")
  check.param(m.sd, "req"); check.param(m.sd, "min", min = 0)
  check.param(pc, "req"); check.param(pc, "unit")
  check.param(pt, "req"); check.param(pt, "unit")
  check.param(sigma.u, "req"); check.param(sigma.u, "pos")

  logoddsc <- log(pc / (1 - pc))
  logoddst <- log(pt / (1 - pt))
  gam0 <- (logoddsc + logoddst) / 2
  gam1 <- logoddst - logoddsc

  sigmasqet <- 2 + exp(-gam0 - 2 * gam1) + exp(gam0 + 2 * gam1)
  sigmasqec <- 2 + exp(-gam0 + 2 * gam1) + exp(gam0 - 2 * gam1)

  alpha_t <- sigmasqet / sigma.u^2
  alpha_c <- sigmasqec / sigma.u^2

  cv <- m.sd / m

  lamt <- m / (m + alpha_t)
  lamc <- m / (m + alpha_c)

  re.num <- (1 - cv^2 * lamt * (1 - lamt)) * (1 - cv^2 * lamc * (1 - lamc)) * (lamt + lamc)
  re.denom <- lamt + lamc - cv^2 * (lamt^2 * (1 - lamt) + lamc^2 * (1 - lamt))
  re <- re.num / re.denom

  return(re)
}


