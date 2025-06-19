#' @title Matrix or Data Frame Union
#' @description A union class that allows for both matrix and data frame types.
#' @export
methods::setClassUnion(
  "matrixOrDataFrame",
  c("matrix", "data.frame")
)

#' @title corrResult Class
#' @description A class to hold correlation and p-values matrices.
#'
#' The `corrResult` class is designed to store the results of a correlation analysis, including
#' the correlation matrix, p-values matrix, and any parameters used during the analysis.
#'
#' @slot correlation_matrix A matrix or data frame containing correlation values.
#' @slot p_values_matrix A matrix or data frame containing p-values associated with the correlation values.
#' @slot parameters A list of parameters used in the analysis of correlation.
#' @export
methods::setClass("corrResult",
  slots = list(
    correlation_matrix = "matrixOrDataFrame",
    p_values_matrix = "matrixOrDataFrame",
    parameters = "list"
  )
)

#' @title Show Method for corrResult
#' @description Displays the correlation matrix, p-values matrix, and parameters used in the analysis.
#'
#' This method provides a concise summary of the `corrResult` object, showing the first few rows and columns
#' of the correlation matrix and p-values matrix, as well as the parameters used for the analysis.
#'
#' @param object An instance of `corrResult`.
#' @importFrom methods show
#' @export
methods::setMethod(
  "show",
  "corrResult",
  function(object) {
    n_col <- min(ncol(object@correlation_matrix), 3)
    n_row <- min(nrow(object@correlation_matrix), 3)

    cat("@correlation_matrix:\ndim: ")
    cat(dim(object@correlation_matrix), " [", n_row, " x ", n_col, "]\n")
    print(object@correlation_matrix[1:n_row, 1:n_col])

    cat("\n@p_values_matrix:\ndim: ")
    cat(dim(object@p_values_matrix), " [", n_row, " x ", n_col, "]\n")
    print(object@p_values_matrix[1:n_row, 1:n_col])

    cat("\n@parameters:\n")
    print(object@parameters)
  }
)
