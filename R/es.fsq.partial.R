#' Cohen f-squared effect size for partial F test in multiple linear regression
#'
#' @description
#' Computes the f-squared (r^2) effect size for a partial F test in a multiple linear regression model
#' based on the model R^2 (Rsq). Based on Cohen (1988). EDIT THIS
#'
#' @details
#' Cohen J (1988) Statistical Power Analysis for the Behavioral Sciences, 2nd edition.
#' Lawrence Erlbaum Associates, Hillsdale, New Jersey
#'
#' @param Rsq.red The squared population multiple correlation coefficient for the reduced model. Either 2 out of 3 Rsq terms OR pc must be specified.
#' @param Rsq.full The squared population multiple correlation coefficient for the full model. Either 2 out of 3 Rsq terms OR pc must be specified.
#' @param Rsq.diff The difference between the squared population multiple correlation coefficient for the full model and the reduced model. Either 2 out of 3 Rsq terms OR pc must be specified.
#' @param pc The partial correlation coefficient. Either 2 out of 3 Rsq terms OR pc must be specified.
#'
#' @return A list of the arguments and the f^2 effect size.
#' @export
#'
#' @examples
#' es.fsq.partial(pc = 0.2)
#' es.fsq.partial(Rsq.red = 0.25, Rsq.full = 0.35)
#' es.fsq.partial(Rsq.red = 0.25, Rsq.diff = 0.1)
#' es.fsq.partial(Rsq.full = 0.35, Rsq.diff = 0.1)

es.fsq.partial <- function (Rsq.red = NULL, Rsq.full = NULL, Rsq.diff = NULL, pc = NULL) {

  # Check if the arguments are specified correctly
  if ((sum(sapply(list(Rsq.red, Rsq.full, Rsq.diff), is.null)) > 1) & is.null(pc))
    stop("please specify 2 out of 3 Rsq parameters OR pc")

  if(is.null(pc)) {
    # Calculate the missing Rsq term
    if (is.null(Rsq.diff))
      Rsq.diff <- Rsq.full - Rsq.red
    else if (is.null(Rsq.red))
      Rsq.red <- Rsq.full - Rsq.diff
    else if (is.null(Rsq.full))
      Rsq.full <- Rsq.diff + Rsq.red

    # Calculate f^2
    fsq <- Rsq.diff / (1 - Rsq.full)

  } else {
    # Calculate f^2
    fsq <- pc^2 / (1 - pc^2)
  }

  # Print output as a power.htest object
  METHOD <- "Cohen f^2 effect size calculation for partial F test"

  if (is.null(pc)) {
    structure(list(Rsq.red = Rsq.red, Rsq.full = Rsq.full,
                   Rsq.diff = Rsq.diff, fsq = fsq,
                   method = METHOD), class = "power.htest")
  } else {
    structure(list(pc = pc, fsq = fsq,
                   method = METHOD), class = "power.htest")
  }
}
