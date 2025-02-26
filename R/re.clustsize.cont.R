#' Relative efficiency of a cluster randomized or multisite trial with continuous outcome
#' with varying cluster sizes
#'
#' @description
#' For a continuous outcome, computes the relative efficiency (ratio of the variances) of a cluster randomized
#' or multisite trial with varying cluster sizes to that of a trial with constant cluster sizes,
#' assuming equal total number of subjects. This approximation may not be valid when m.sd / m > 2.
#'
#'
#' @param m The number of subjects per cluster or the mean cluster size (if unequal number of participants per cluster).
#' @param m.sd The standard deviation of cluster sizes (in case of unequal number of participants per cluster).
#' @param icc The intraclass correlation coefficient. For a multisite trial this is icc1. For a CRT this is the average of the 2 icc's.
#' @return The computed RE.
#' @export
#'
#' @examples
#' re.clustsize.cont(m = 25, m.sd = 15, icc = 0.05)


re.clustsize.cont <- function (m, m.sd, icc) {
  check.param(m, "req"); check.param(m, "pos")
  check.param(m.sd, "req"); check.param(m.sd, "min", min = 0)
  check.param(icc, "req"); check.param(icc, "uniti")

  cv <- m.sd / m
  K <- (m * icc) / (1 + (m - 1) * icc)
  RE <- 1 - cv^2 * K * (1 - K)

  NOTE <- "A negative RE is invalid. The approximation is inaccurate when m.sd / m > 2."
  if (RE < 0) print(paste("NOTE:", NOTE))

  return(RE)
}
