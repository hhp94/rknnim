#' Inject NA Values into a Matrix
#'
#' This function randomly selects positions in a matrix to inject NA values, ensuring that
#' the injection does not exceed specified missingness thresholds for rows and columns.
#' It attempts to find a valid set of positions within a maximum number of iterations.
#'
#' @param obj A matrix or a character string specifying the path to a file to be read
#'   using `safe_read()`. The matrix should contain data where NAs can be injected.
#' @param hash Optional hash value for validation when reading from a file using `safe_read()`.
#' @param prop Proportion of non-NA elements to convert to NA. Must be between 0 and 1 (exclusive of 0).
#' @param min_na Minimum number of NAs to inject. Must be a positive integer.
#' @param max_na Maximum number of NAs to inject. Must be a positive integer and at least `min_na`.
#' @param c_thresh Maximum allowed proportion of missing values per column after injection (default: 0.9).
#' @param r_thresh Maximum allowed proportion of missing values per row after injection (default: 0.9).
#' @param max_iter Maximum number of iterations to attempt finding valid NA positions (default: 1000).
#'
#' @return A vector of integer indices indicating the positions in the matrix
#' where NAs should be injected.
#'
#' @details
#' The function first computes the number of NAs to inject based on `prop`, clamped between `min_na` and `max_na`.
#' It then repeatedly samples random positions from existing non-NA elements and checks if injecting NAs
#' at those positions would exceed the missingness thresholds for any row or column (accounting for existing NAs).
#' If no valid set is found within `max_iter` attempts, an error is thrown.
#'
#' @examples
#' mat <- matrix(1:100, nrow = 10, ncol = 10)
#' na_positions <- inject_na(mat, prop = 0.1, min_na = 5, max_na = 15)
#' mat[na_positions] <- NA
#'
#' @export
inject_na <- function(
    obj,
    hash = NULL,
    prop,
    min_na,
    max_na,
    c_thresh = 0.9,
    r_thresh = 0.9,
    max_iter = 1000
  ) {
  if (is.character(obj)) {
    obj <- safe_read(obj, hash = hash)
  }
  checkmate::assert_number(c_thresh, lower = 0, upper = 1)
  checkmate::assert_number(r_thresh, lower = 0, upper = 1)
  checkmate::assert_count(max_iter, positive = TRUE)
  checkmate::assert_number(prop, lower = 0, upper = 1)
  checkmate::assert_true(prop > 0, .var.name = "prop > 0")
  checkmate::assert_count(min_na, positive = TRUE)
  checkmate::assert_count(max_na, positive = TRUE)
  checkmate::assert_true(max_na >= min_na, .var.name = "max_na >= min_na")

  # subset the matrix to the specified features and samples
  na_mat <- !is.na(obj)
  not_na <- which(na_mat)
  # calculate `na_size` and make sure falls within [min_na, max_na]
  if (min_na > length(not_na)) {
    stop(
      sprintf(
        "'min_na' (%d) exceeds the number of available non-NA elements (%d).
        Adjust 'min_na' or increase feature/sample group size.",
        min_na, length(not_na)
      )
    )
  }
  na_size <- ceiling(length(not_na) * prop)
  na_size <- max(min(na_size, max_na), min_na)
  if (na_size == 0 || na_size > length(not_na)) {
    stop(
      sprintf(
        "Invalid number of NAs to inject: calculated na_size = %d, available non-NA elements = %d.
        Adjust 'prop', 'min_na', 'max_na', or increase feature/sample group size.",
        na_size, length(not_na)
      )
    )
  }
  # max allowed missing counts per column and row
  max_col_miss <- floor(nrow(na_mat) * c_thresh)
  max_row_miss <- floor(ncol(na_mat) * r_thresh)
  # initialize variables for the while loop
  c_miss <- TRUE
  r_miss <- TRUE
  na_loc <- NULL
  iter <- 0
  # inject NAs while ensuring missingness thresholds and iter are not exceeded
  while (c_miss || r_miss) {
    iter <- iter + 1
    if (iter > max_iter) {
      stop("NA injection failed. Lower 'prop' or increase 'c_thresh' and 'r_thresh'.")
    }
    na_mat_test <- na_mat
    na_loc <- sample(not_na, size = na_size)
    na_mat_test[na_loc] <- FALSE
    # calculate the counts of missing values in columns and rows
    col_miss_count <- nrow(na_mat_test) - colSums(na_mat_test)
    row_miss_count <- ncol(na_mat_test) - rowSums(na_mat_test)
    # check if any column or row exceeds the missingness thresholds
    c_miss <- any(col_miss_count > max_col_miss)
    r_miss <- any(row_miss_count > max_row_miss)
  }
  return(na_loc)
}

safe_read <- function(obj, hash = NULL) {
  obj <- qs2::qs_read(obj, validate_checksum = TRUE)
  if (!is.null(hash)) {
    stopifnot("'hash' mismatch" = rlang::hash(obj) == hash)
  }
  return(obj)
}
