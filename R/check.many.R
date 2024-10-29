#' Internal-use function for quality checking sets of parameters
#'
#' @param paramlist The list of parameters being checked.
#' @param type The expected type of parameter: currently only supports "oneof".
#'
#' @return If the check passes, returns nothing. If the check does not pass, throw an error.
#' @keywords internal
#' @export
#'
#' @examples
#' N <- 10
#' power <- NULL
#' check.many(list(N, power), "oneof")

check.many <- function(paramlist = NULL,
                       type = c("oneof")) {

  names <- lapply(substitute(paramlist), deparse)
  names <- paste(names[-1], collapse = ", ")

  # variables where only one must be missing
  if (type == "oneof") {
    if (sum(sapply(paramlist, is.null)) != 1)
      stop(paste("exactly one of the following arguments must be NULL:", names))
  }
}
