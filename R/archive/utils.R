
#' Melt a matrix into a long data.table
#'
#' This function converts a matrix into a long format data.table. It first transforms the matrix into a data.table,
#' then melts it into long format, cleaning up and converting factors and removing non-finite or non-positive values.
#'
#' @param x A matrix that will be converted into a long format data.table.
#'
#' @return A long format data.table containing the melted matrix with 'x', 'y', and 'value' columns.
#' @export
#'
#' @examples
#' # Assuming 'mat' is a matrix:
#' long_dt <- melt_matrix(mat)
melt_matrix <- function(x, remove_zero = TRUE, remove_non_finite = TRUE) {
  has_names <- is.character(colnames(x))
  # Convert the matrix to a data.table
  dt <- data.table::as.data.table(x, na.rm = FALSE)
  dt[, x := 1:.N]
  # Melt the data.table
  long_dt <- melt(dt, id.vars = "x", variable.name = "y", value.name = "value",
                  variable.factor = FALSE)
  if (isFALSE(has_names)) long_dt[, y := as.integer(gsub("V", "", y))]

  if (isTRUE(remove_zero))  long_dt <- long_dt[value > 0]
  if (isTRUE(remove_non_finite))  long_dt <- long_dt[is.finite(value)]

  return(long_dt[])
}


#' Sum Columns of Matrices with Matching Column Names
#'
#' Takes a list of matrices and returns a single matrix with the same number
#' of rows and unique column names across all input matrices.
#' For each column name present in any of the matrices, it creates a column
#' in the resulting matrix, where each entry is the sum of entries from all
#' matrices with the same column name.
#'
#' @param list_of_matrices A list where each element is a matrix of numeric values.
#' Each matrix must have the same number of rows.
#'
#' @return A matrix with the same number of rows as the input matrices and as many columns as there are unique column names across all matrices. The values are the sums of corresponding entries for each column name.
#'
#' @examples
#' sample_comb <- list(
#'   matrix(c(1, 2, 3, 4, 5, 6), nrow = 2),
#'   matrix(c(7, 8, 9, 10, 11, 12), nrow = 2, dimnames = list(NULL, c("V1", "V2", "V3")))
#' )
#' sumMatrices(sample_comb)
#'
#' @export

sum_matrices <- function(list_of_matrices) {
  # Ensure all matrices have the same number of rows
  n_rows <- unique(sapply(list_of_matrices, nrow))
  if(length(n_rows) != 1) {
    stop("All matrices must have the same number of rows.")
  }

  # Create an empty list to store matrices with standardized column names
  standardized_matrices <- vector("list", length(list_of_matrices))

  # Standardize the matrices to have all unique column names
  all_col_names <- unique(unlist(lapply(list_of_matrices, colnames)))

  for (i in seq_along(list_of_matrices)) {
    mat <- list_of_matrices[[i]]
    standardized_mat <- matrix(0, nrow = n_rows, ncol = length(all_col_names))
    colnames(standardized_mat) <- all_col_names

    for (col in colnames(mat)) {
      standardized_mat[, col] <- mat[, col]
    }
    standardized_matrices[[i]] <- standardized_mat
  }

  # Sum all standardized matrices
  summed_matrix <- Reduce(`+`, standardized_matrices)
  return(summed_matrix)
}
