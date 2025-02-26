#' Internal-use wrapper for catching errors induced by stats::uniroot
#'
#' @param f The function for which the root is sought
#' @param interval A vector containing the end-points of the interval to be searched for the root
#'
#' @return If no error occurs, returns the solved root. If an error occurs, output a custom message.
#' @keywords internal
#' @export
#'
#' @examples
#' f <- function(x) x - 3
#' safe.uniroot(f, c(-10, 10))

safe.uniroot <- function(f = NULL, interval = NULL) {
  tryCatch(
    {
      stats::uniroot(f, interval)$root
    },
    error = function(e) {
      stop(simpleError(paste0(e$message, ".
The uniroot function is used to solve the power equation for unknowns.
This error occurs when it cannot find a solution or has found multiple
solutions. Please try again with different parameters.")))
    }
  )
}
