#' Internal-use function for quality checking individual parameters
#'
#' @param param The parameter being checked.
#' @param type The expected type of parameter: one of "req", "num", "unit(i)(i)", "pos", "int", "bool", "vals", "min(i)", "vec", "mat", "sum".
#' @param valslist If the type is "vals", a list of allowed values.
#' @param min If the type is "min" or "mini", a lower bound.
#'
#' @return If the check passes, returns nothing. If the check does not pass, throw an error.
#' @keywords internal
#' @export
#'
#' @examples
#' power <- 0.8
#' check.param(power, "unit")

check.param <- function(param = NULL,
                  type = c("req",   # required
                           "num",   # numerical
                           "unit",  # (0, 1)
                           "uniti", # [0, 1)
                           "unitii",# [0, 1]
                           "pos",   # positive values
                           "int",   # positive integers
                           "bool",  # T or F
                           "vals",  # certain values (specify valslist)
                           "min",   # can't be < minimum (specify min)
                           "mini",  # can't be <= minimum (specify min)
                           "vec",   # numerical vector
                           "mat",   # numerical matrix
                           "sum"    # sum to one
                           ),
                  valslist = NULL, min = NULL) {

  name <- deparse(substitute(param))

  # required variable
  if (type == "req") {
    if (is.null(param))
      stop(paste(name, "must be specified"))
  }

  # the following checks are done assuming the variable isn't null
  else {
    if (!is.null(param)){

      # numeric variables
      if (type == "num") {
        if (length(param) != 1)
          stop(paste(name, "should have single value (length 1)"))
        if (!is.numeric(param))
          stop(paste(name, "should be a numeric value"))
      }

      # variables between 0 and 1
      else if (type == "unit") {
        if (length(param) != 1)
          stop(paste(name, "should have single value (length 1)"))
        if (!is.numeric(param))
          stop(paste(name, "should be a numeric value"))
        if (param >= 1 || param <= 0)
          stop(paste(name, "should be between 0 and 1"))
      }

      # variables between 0 and 1 including 0
      else if (type == "uniti") {
        if (length(param) != 1)
          stop(paste(name, "should have single value (length 1)"))
        if (!is.numeric(param))
          stop(paste(name, "should be a numeric value"))
        if (param >= 1 || param < 0)
          stop(paste(name, "should be between 0 and 1 (including 0)"))
      }

      # variables between 0 and 1 including 0 and 1
      else if (type == "unitii") {
        if (length(param) != 1)
          stop(paste(name, "should have single value (length 1)"))
        if (!is.numeric(param))
          stop(paste(name, "should be a numeric value"))
        if (param > 1 || param < 0)
          stop(paste(name, "should be between 0 and 1 (including 0 and 1)"))
      }

      # positive values only
      else if (type == "pos") {
        if (length(param) != 1)
          stop(paste(name, "should have single value (length 1)"))
        if (!is.numeric(param))
          stop(paste(name, "should be a numeric value"))
        if (param <= 0)
          stop(paste(name, "should be a positive numeric value"))
      }

      # integers only
      else if (type == "int") {
        if (length(param) != 1)
          stop(paste(name, "should have single value (length 1)"))
        if (!is.numeric(param))
          stop(paste(name, "should be a numeric value"))
        if (param %% 1 != 0)
          stop(paste(name, "should be a positive integer"))
        if (param < 0)
          stop(paste(name, "should be positive"))
      }

      # true / false only
      else if (type == "bool") {
        if (length(param) != 1)
          stop(paste(name, "should have single value (length 1)"))
        if (!isTRUE(param) & !isFALSE(param))
          stop(paste(name, "should be either TRUE or FALSE"))
      }

      # only certain values allowed
      else if (type == "vals") {
        if (length(param) != 1)
          stop(paste(name, "should have single value (length 1)"))
        if (!(param %in% valslist))
          stop(paste(name, "should be one of:", paste(valslist, collapse = ", ")))
      }

      # value has a minimum
      else if (type == "min") {
        if (length(param) != 1)
          stop(paste(name, "should have single value (length 1)"))
        if (!is.numeric(param))
          stop(paste(name, "should be a numeric value"))
        if (param < min)
          stop(paste(name, "should be greater than or equal to", min))
      }

      # value has a minimum but includes it
      else if (type == "mini") {
        if (length(param) != 1)
          stop(paste(name, "should have single value (length 1)"))
        if (!is.numeric(param))
          stop(paste(name, "should be a numeric value"))
        if (param <= min)
          stop(paste(name, "should be greater than", min))
      }

      # value is numerical vector
      else if (type == "vec") {
        if (!is.numeric(param))
          stop(paste(name, "should be a numeric vector"))
        if (length(param) < 2)
          stop(paste(name, "should have at least 2 groups"))
      }

      # value is numerical matrix
      else if (type == "mat") {
        if (!is.numeric(c(param)))
          stop(paste(name, "should be a numeric matrix"))
        if (nrow(param) < 2 | ncol(param) < 2)
          stop(paste(name, "each factor should have at least 2 groups"))
      }


      # value is numerical vector/matrix that sums to one
      else if (type == "sum") {
        if (any(param <= 0) | any(param >= 1))
          stop(paste("all values of", name, "should be between 0 and 1"))
        if (sum(param) != 1)
          stop(paste("values of", name, "should sum to one"))
      }
    }
  }
}
