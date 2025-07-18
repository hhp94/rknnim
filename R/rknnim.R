#' @title Rolling KNN Imputation
#'
#' @description
#' Performs rolling window KNN imputation on a matrix with missing values.
#' The function divides the input matrix into overlapping windows along the rows,
#' imputes each window using KNN imputation, and combines the results by averaging
#' values in overlapping regions. Any remaining NA values are imputed with row means.
#'
#' @param obj A numeric matrix where rows represent features and columns represent samples.
#' @param n.feat Integer specifying the number of features (rows) per imputation window.
#' @param n.overlap Integer specifying the number of features of overlap
#'   between consecutive windows (default: 10).
#' @param coord Optional coordinate matrix (default: NULL). Currently unused in the function.
#' @inheritParams impute::impute.knn
#' @inheritParams purrr::map
#' @param .parallel Logical indicating whether to use parallel processing for imputation
#'   (default: FALSE). Requires the 'mirai' package if enabled.
#'
#' @return A numeric matrix of the same dimensions as \code{obj} with missing values imputed.
#'
#' @export
#'
#' @examples
#' \dontrun{
#' # Example matrix with missing values
#' set.seed(123)
#' mat <- matrix(rnorm(100), nrow = 20, ncol = 5)
#' mat[sample(100, 10)] <- NA
#' imputed <- rknnim(mat, n.feat = 10, overlap = 0.2, k = 5)
#' }
rknnim <- function(
    obj,
    n.feat,
    n.overlap = 10,
    coord = NULL,
    k = 10,
    rowmax = 0.5,
    colmax = 0.8,
    rng.seed = 362436069,
    .progress = FALSE,
    .parallel = FALSE) {
  # pre-conditioning
  checkmate::assert_matrix(obj, mode = "numeric", null.ok = FALSE)
  if (nrow(obj) < ncol(obj)) {
    warning("Input object is wide, this is unexpected. Rows should be features and samples should be columns")
  }
  checkmate::assert_true(sum(is.infinite(obj)) == 0)
  checkmate::assert_integerish(n.feat, lower = 2, upper = nrow(obj), len = 1, null.ok = FALSE)
  checkmate::assert_integerish(n.overlap, lower = 0, upper = n.feat - 1, len = 1, null.ok = FALSE)

  # To-dos: weighted, write intermediary to files, C++
  # C++ version would be idx 0
  idx <- 1

  # max_i > (nrow(obj) - idx)/(n.feat - n.overlap)
  max_step <- ceiling((nrow(obj) - idx) / (n.feat - n.overlap))
  step <- 0:max_step
  start <- idx + (step * n.feat) - (step * n.overlap)
  end <- start + n.feat - 1
  # Edge case, the end might overshoot the number of rows of the obj
  n_overshoot <- sum(end > nrow(obj))

  # In which case trim off the runs that overshoot
  corrected_length <- length(end) - n_overshoot
  start <- start[idx:corrected_length]
  end <- end[idx:corrected_length]

  # And make the last window extra wide to cover the full end
  end[corrected_length] <- nrow(obj)
  width <- end - start + 1

  # Rolling Imputation
  if (.parallel) {
    mirai::require_daemons()
    fn <- purrr::in_parallel
  } else {
    fn <- function(x, ...) {
      x
    }
  }
  imputed <- purrr::map(
    seq_along(start),
    fn(
      function(i) {
        impute_knn(
          obj = obj[start[i]:end[i], ],
          k = k,
          rowmax = rowmax,
          colmax = colmax,
          rng.seed = rng.seed
        )
      },
      impute_knn = impute_knn,
      obj = obj,
      start = start,
      end = end,
      k = k,
      rowmax = rowmax,
      colmax = colmax,
      rng.seed = rng.seed
    ),
    .progress = .progress
  )

  # To get the final matrix. We add a zero matrix with the imputed values, rows
  # that are outside of the overlap will be added once, while outside of the overlap
  # will be added twice, then we need the count matrix to calculate the average

  # Initialize a matrix to hold the final results and a matrix to count the contributions
  final_imputed <- matrix(0, nrow = nrow(obj), ncol = ncol(obj))
  counts <- matrix(0, nrow = nrow(obj), ncol = ncol(obj))

  # Loop through each imputed window.
  for (i in seq_along(imputed)) {
    # Row range for the current window
    rows_to_update <- start[i]:end[i]

    # Add the imputed values from the window to the accumulator matrix
    final_imputed[rows_to_update, ] <- final_imputed[rows_to_update, ] + imputed[[i]]

    # Increment the count for these rows, noting that they've been imputed
    counts[rows_to_update, ] <- counts[rows_to_update, ] + 1
  }

  # Average the values by dividing the sums by the counts.
  final_imputed <- final_imputed / counts

  # Chunking windows can create NA cols, this final step impute the NA cols with
  # the mean of the rows
  if (anyNA(final_imputed)) {
    final_imputed <- mean_impute_row(final_imputed)
  }

  colnames(final_imputed) <- colnames(obj)
  row.names(final_imputed) <- row.names(obj)

  return(final_imputed)
}

#' @title KNN Imputation Wrapper
#'
#' @description
#' A wrapper around \code{impute::impute.knn} that filters out columns with excessive
#' missing values before performing KNN imputation. If all columns are filtered out,
#' the original matrix is returned unchanged.
#'
#' @inheritParams rknnim
#'
#' @return The imputed matrix, with good columns imputed and bad columns unchanged.
#'
#' @importFrom impute impute.knn
#'
#' @examples
#' \dontrun{
#' mat <- matrix(c(1, NA, 3, 4, 5, NA), nrow = 2)
#' imputed <- impute_knn(mat, k = 1, rowmax = 0.5, colmax = 0.8, rng.seed = 123)
#' }
#' @keywords internal
#' @noRd
impute_knn <- function(obj, k, rowmax, colmax, rng.seed) {
  na_mat <- is.na(obj)
  good_cols <- !(colSums(na_mat) / nrow(na_mat) >= colmax)
  if (sum(good_cols) == 0) {
    return(obj)
  }
  imputed_good <- impute.knn(
    data = obj[, good_cols, drop = FALSE],
    k = k,
    maxp = nrow(obj),
    rowmax = rowmax,
    colmax = colmax,
    rng.seed = rng.seed
  )$data

  result <- obj
  result[, good_cols] <- imputed_good
  return(result)
}

#' @title Row Mean Imputation
#'
#' @description
#' Imputes missing values (NA) in a matrix by replacing them with the mean of their
#' respective rows, computed excluding NA values.
#'
#' @inheritParams rknnim
#'
#' @return The matrix with NA values replaced by row means.
#'
#' @export
#'
#' @examples
#' \dontrun{
#' mat <- matrix(c(1, NA, 3, 4, 5, NA), nrow = 2)
#' imputed <- mean_impute_row(mat)
#' # Result: row 1 mean ( (1+3)/2 = 2 ) for NA; row 2 mean ( (4+5)/2 = 4.5 ) for NA
#' }
mean_impute_row <- function(obj) {
  na_indices <- which(is.na(obj), arr.ind = TRUE)
  row_means <- rowMeans(obj, na.rm = TRUE)
  obj[na_indices] <- row_means[na_indices[, 1]]
  return(obj)
}

dummy <- function() {
  Rcpp::evalCpp()
  carrier::crate()
}
