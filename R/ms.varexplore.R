#' Variance exploration for a multisite trial with a binary outcome
#'
#' @description
#' This function can be used to help select a plausible value for the variance/SD of the
#' random intercept for site in a multisite trial with a binary outcome.
#' Based on user-supplied values of the outcome proportions in the two conditions,
#' this function outputs, for a range of possible values of the SD, the intervals
#' within which we expect about 95% of the site-specific
#' odds ratios to lie.
#'
#' @details
#' In a multisite trial design, participants are randomized to conditions within site.
#' The use of this function is illustrated in Crespi CM (2025) Power and Sample Size in R.
#'
#'
#'
#' @param pc The probability of the outcome in control clusters.
#' @param pt The probability of the outcome in treatment clusters.
#' @param print Whether or not to print the results in a table; defaults to TRUE. To retrieve the results even when print = FALSE, assign the function output to an object.
#'
#' @return A list containing the calculated OR and a dataframe of outputs.
#' @import knitr
#' @export
#'
#' @examples ms.varexplore(pc = 0.1, pt = 0.2)
#' output <- ms.varexplore(pc = 0.1, pt = 0.2, print = FALSE)
#' output$OR

ms.varexplore <- function(pc, pt, print = TRUE) {

  or <- (pt / (1 - pt)) / (pc / (1 - pc))

  logoddsc <- log(pc / (1 - pc))
  logoddst <- log(pt / (1 - pt))
  gam1 <- logoddst - logoddsc

  sigma_u1 <- seq(0.1, 1, by = 0.1)
  or_lo <- exp(gam1 - 1.96 * sigma_u1)
  or_lo <- round(or_lo, 2)
  or_up <- exp(gam1 + 1.96 * sigma_u1)
  or_up <- round(or_up, 2)

  table <- data.frame("sigma.u1" = sigma_u1, "OR.lower" = or_lo, "OR.upper" = or_up)

  if (print) {
    out <- kable(table, caption = paste("OR:", or), "simple")
    print(gsub("^Table:", "", out))
  }

  return(invisible(list(OR = or, table = table)))
}
