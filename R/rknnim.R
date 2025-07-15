rknnim.internal <- function(
    obj,
    n.feat,
    overlap = 0.1,
    coord = NULL,
    k = 10,
    rowmax = 0.5,
    colmax = 0.8,
    rng.seed = 362436069) {
  checkmate::assert_matrix(obj, mode = "numeric")
  checkmate::assert_true(sum(is.infinite(khanmiss1)) == 0)

  n.overlap <- round(n.feat * overlap)

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
  ## Initialize
  imputed <- lapply(
    width,
    \(x) {
      matrix(NA, nrow = x, ncol = ncol(obj))
    }
  )

  ## Impute
  for (i in seq_along(start)) {
    imputed[[i]][, ] <- impute_knn(
      obj = obj[start[i]:end[i], ],
      k = k,
      rowmax = rowmax,
      colmax = colmax,
      rng.seed = rng.seed
    )
  }

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
  if(anyNA(final_imputed)) {
    final_imputed <- mean_impute_row(final_imputed)
  }

  colnames(final_imputed) <- colnames(obj)
  row.names(final_imputed) <- row.names(obj)

  return(final_imputed)
}

impute_knn <- function(obj, k, rowmax, colmax, rng.seed) {
  na_mat <- is.na(obj)
  good_cols <- !(colSums(na_mat) / nrow(na_mat) >= colmax)
  if (sum(good_cols) == 0) {
    return(obj)
  }
  imputed_good <- impute::impute.knn(
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

mean_impute_row <- function(obj) {
  na_indices <- which(is.na(obj), arr.ind = TRUE)
  row_means <- rowMeans(obj, na.rm = TRUE)
  obj[na_indices] <- row_means[na_indices[, 1]]
  return(obj)
}
