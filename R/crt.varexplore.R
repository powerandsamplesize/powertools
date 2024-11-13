#' Variance exploration for cluster randomized trials with binary outcomes
#'
#' @description
#' This function can be used to help select a plausible value for the variance/SD of the
#' random intercept for cluster in a cluster randomized trial with a binary outcome.
#' Based on user-supplied values of the outcome proportions in the two arms,
#' the function outputs, for a range of possible values of the SD of the
#' random intercept, the intervals within which we expect about 95% of the cluster-level
#' proportions to lie in each arm.
#'
#' @details
#' The use of this function is illustrated in Crespi CM (2025) Power and Sample Size in R.
#'
#'
#'
#' @param pc The probability of the outcome in control clusters.
#' @param pt The probability of the outcome in treatment clusters.
#' @param print Whether or not to print the results in a table; defaults to TRUE. To retrieve the results even when print = FALSE, assign the function output to an object.
#'
#' @return A list of the arguments and a dataframe of outputs.
#' @import knitr
#' @export
#'
#' @examples crt.varexplore(pc = 0.25, pt = 0.15)
#' output <- crt.varexplore(pc = 0.25, pt = 0.15, print = FALSE)
#' output$pc

crt.varexplore <- function(pc, pt, print = TRUE){
  logoddsc <- log(pc / (1 - pc))
  logoddst <- log(pt / (1 - pt))
  gam0 <- (logoddsc + logoddst) / 2
  gam1 <- logoddst - logoddsc

  sigma.u <- seq(0.1, 1, by = 0.1)
  pc.lo <- exp(gam0 + gam1 * -0.5 - 1.96 * sigma.u)/(1 + exp(gam0 + gam1 * -0.5 - 1.96 * sigma.u))
  pc.lo <- round(pc.lo, 2)
  pc.up <- exp(gam0 + gam1 * -0.5 + 1.96 * sigma.u)/(1 + exp(gam0 + gam1 * -0.5 + 1.96 * sigma.u))
  pc.up <- round(pc.up, 2)

  pt.lo <- exp(gam0 + gam1 * 0.5 - 1.96 * sigma.u)/(1 + exp(gam0 + gam1 * 0.5 - 1.96 * sigma.u))
  pt.lo <- round(pt.lo, 2)
  pt.up <- exp(gam0 + gam1 * 0.5 + 1.96 * sigma.u)/(1 + exp(gam0 + gam1 * 0.5 + 1.96 * sigma.u))
  pt.up <- round(pt.up, 2)

  table <- data.frame("sigma.u" = sigma.u, "pc.lower" = pc.lo, "pc.upper" = pc.up,
                      "pt.lower" = pt.lo, "pt.upper" = pt.up)

  if (print) {
    out <- kable(table, caption = paste("pc:", pc, "; pt:", pt), "simple")
    print(gsub("^Table:", "", out))
  }

  return(invisible(list(pc = pc, pt = pt, table = table)))
}

