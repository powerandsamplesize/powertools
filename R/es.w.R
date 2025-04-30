#' Cohen w effect size calculation for chi-square tests
#'
#' @description
#' Calculates the w effect size for chi-square tests. For chi-square goodness of fit tests, specify
#' p0vec and p1vec. For chi-square tests of independence, provide a matrix of cell probabilities.
#' Based on Cohen (1988).
#'
#' @details
#' Cohen J (1988) Statistical Power Analysis for the Behavioral Sciences, 2nd edition.
#' Lawrence Erlbaum Associates, Hillsdale, New Jersey
#'
#'
#'
#' @param p0vec The first vector of probabilities. Both p0vec and p1vec, or pmatrix must be specified.
#' @param p1vec The second vector of probabilities. Both p0vec and p1vec, or pmatrix must be specified.
#' @param pmatrix The two-way probability table. Both p0vec and p1vec, or pmatrix must be specified.
#'
#' @return A list of the arguments and the w effect size.
#' @export
#'
#' @examples
#' es.w(p0vec = c(0.5, 0.3, 0.2), p1vec = c(0.7, 0.2, 0.1))
#' es.w(pmatrix = matrix(c(0.050, 0.350, 0.100, 0.075, 0.250, 0.175), nrow = 2, byrow = TRUE))

es.w <- function (p0vec = NULL, p1vec = NULL, pmatrix = NULL) {

  # Check if the arguments are specified correctly
  if(is.null(pmatrix) & (is.null(p0vec) | is.null(p1vec)))
    stop("p0vec and p1vec, or pmatrix must be specified")
  if(!is.null(pmatrix) & (!is.null(p0vec) | !is.null(p1vec)))
    stop("p0vec and p1vec, or pmatrix must be specified")
  check.param(p0vec, "sum1")
  check.param(p1vec, "sum1")
  check.param(pmatrix, "sum1")

  # Calculate w
  if (is.null(pmatrix))
    w <- sqrt(sum((p1vec - p0vec)^2 / p0vec))
  else {
    pi <- apply(pmatrix, 1, sum)
    pj <- apply(pmatrix, 2, sum)
    p0 <- pi %*% t(pj)
    w <- sqrt(sum((pmatrix - p0)^2/p0))
  }

  # Print output as a power.htest object
  METHOD <- "Cohen w effect size calculation for chi-square tests"

  if (is.null(pmatrix))
    structure(list(p0vec = p0vec, p1vec = p1vec, w = w,
                 method = METHOD), class = "power.htest")
  else
    structure(list(pmatrix = matrix.format(pmatrix), w = w,
                   method = METHOD), class = "power.htest")

}
