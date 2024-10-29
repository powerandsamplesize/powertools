#' Correlation between a cluster mean at baseline and follow up
#'
#' @description
#' For a cluster randomized trial with a continuous outcome, this function calculates the correlation
#' between a cluster's mean at baseline and at follow up based on various inputs.
#' For cross-sectional sampling of subjects,
#' that is, different subjects are measured at baseline and follow up, specify sac = 0.
#'
#'
#' @param m The number of measurements in each cluster at baseline and follow up.
#' @param icc The intraclass correlation coefficient.
#' @param cac The cluster autocorrelation.
#' @param sac The subject autocorrelation.
#' @return The computed correlation.
#' @export
#'
#' @examples
#' crt.means.r(m = 30, icc = 0.05, cac = 0.4, sac = 0.5)


crt.means.r <- function (m, icc, cac, sac) {
  check.param(m, "req"); check.param(m, "pos")
  check.param(icc, "req"); check.param(icc, "uniti")
  check.param(cac, "req"); check.param(cac, "uniti")
  check.param(sac, "req"); check.param(sac, "uniti")

  denom <- 1 + (m - 1) * icc
  num <- m * icc * cac + (1 - icc) * sac
  r <- num / denom

  return(r)
}
