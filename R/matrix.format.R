#' Internal-use function for displaying matrices in function outputs
#'
#' @param matrix The matrix to be displayed.
#' @param nspaces The number of white spaces before the matrix; defaults to 18.
#'
#' @return The matrix printed as a string with rows separated by newline.
#' @keywords internal
#' @export
#'
#' @examples
#' matrix <- matrix(c(1, 2, 3, 4), nrow = 2)
#' matrix.format(matrix = matrix)

matrix.format <- function (matrix = NULL, nspaces = 18) {
  mrows <- c()
  for (i in 1:nrow(matrix))
    mrows <- c(mrows, paste(matrix[i,], collapse = ', '))
  div <- paste0("\n", strrep(" ", nspaces))
  paste(mrows, collapse = div)
}


