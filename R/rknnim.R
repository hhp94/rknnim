rknnim.internal <- function(
    data,
    n.feat,
    overlap,
    coord = NULL,
    k,
    rowmax = NULL,
    colmax = NULL,
    rng.seed = NULL) {
  n.overlap <- round(n.feat * overlap)

  # To-dos: weighted, handle rowmax colmax, write intermediary to files, C++
  # C++ version would be idx 0
  idx <- 1

  # max_i > (nrow(data) - idx)/(n.feat - n.overlap)
  max_step <- ceiling((nrow(data) - idx) / (n.feat - n.overlap))
  step <- 0:max_step
  start <- idx + (step * n.feat) - (step * n.overlap)
  end <- start + n.feat - 1
  # Edge case, the end might overshoot the number of rows of the data
  n_overshoot <- sum(end > nrow(data))

  # In which case trim off the runs that overshoot
  corrected_length <- length(end) - n_overshoot
  start <- start[idx:corrected_length]
  end <- end[idx:corrected_length]

  # And make the last window extra wide to cover the full end
  end[corrected_length] <- nrow(data)
  width <- end - start + 1

  # Rolling Imputation
  ## Initialize
  imputed <- lapply(
    width,
    \(x) {
      matrix(NA, nrow = x, ncol = ncol(data))
    }
  )

  ## Impute
  for (i in seq_along(start)) {
    imputed[[i]][, ] <- impute::impute.knn(
      data = data[start[i]:end[i], ],
      k = k,
      maxp = width[i]
    )$data
  }

  # To get the final matrix. We add a zero matrix with the imputed values, rows
  # that are outside of the overlap will be added once, while outside of the overlap
  # will be added twice, then we need the count matrix to calculate the average

  # Initialize a matrix to hold the final results and a matrix to count the contributions
  final_imputed <- matrix(0, nrow = nrow(data), ncol = ncol(data))
  counts <- matrix(0, nrow = nrow(data), ncol = ncol(data))

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
  colnames(final_imputed) <- colnames(data)
  row.names(final_imputed) <- row.names(data)

  return(final_imputed)
}
