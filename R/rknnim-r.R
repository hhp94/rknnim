r_impute <- function(obj, k, colmax = 1, rowmax = 0.5, row_imp = TRUE) {
  # Pre-conditioning
  stopifnot(
    "Can't have greater or equal number of neighbors to number of rows" = k < nrow(obj)
  )
  # Pre-processing
  n <- nrow(obj)
  p <- ncol(obj)
  miss <- is.na(obj)
  if (sum(miss) == 0) {
    message("No missing data")
    return(obj)
  }
  cmiss <- colSums(miss) / n
  if (any(cmiss >= colmax)) {
    stop("Column(s) missing exceeded colmax. Consider removing column(s) with too high NA %")
  }
  rmiss <- rowSums(miss) / p
  # Partition
  mean_imp <- rmiss >= rowmax
  mean_imp_rows <- obj[mean_imp, ] # Mean Imputed Rows
  pre_imp_rows <- obj[!mean_imp, ] # KNN Imputed Rows
  pre_imp_miss <- miss[!mean_imp, ] # KNN Imputed Miss Object
  pre_imp_rows[pre_imp_miss] <- 0 # Pre fill NA with 0 to sum with imputed values after
  if (any(colSums(pre_imp_miss) / nrow(pre_imp_miss) >= colmax)) {
    stop("Column(s) missing exceeded colmax. Consider removing column(s) with too high NA %")
  }
  # Impute
  post_imp <- r_impute_knn(pre_imp_rows, pre_imp_miss, k = k)
  # Final
  ## Reconstruct and impute values still missing
  final <- rbind(mean_imp_rows, post_imp)[row.names(obj), colnames(obj)]
  na_indices <- which(is.na(final), arr.ind = TRUE)
  final[na_indices] <- if (row_imp) {
    rowMeans(obj, na.rm = TRUE)[na_indices[, 1]]
  } else {
    colMeans(obj, na.rm = TRUE)[na_indices[, 2]]
  }
  return(final)
}

# Rows to calculate distance from excluding self
row_index_r <- function(n) {
  v_mat <- matrix(1:(n - 1), nrow = n, ncol = n - 1, byrow = TRUE)
  v_mat <- v_mat + (col(v_mat) >= row(v_mat))
  v_mat
}

# Euclidian distant. If implement other distant then implement here
calc_distant <- function(v1, v2, m1, m2) {
  valid <- (m1 + m2) == 0
  total_valid <- sum(valid)
  if (total_valid == 0) {
    return(Inf)
  }
  dist <- sum((v1[valid] - v2[valid])^2)
  return(dist / total_valid)
}

# Just placeholder for the implementation in C++ with partial sort, return index
# of rows that are closest neighbors
knn_partial_sort <- function(distance, k) {
  sorted <- sort(distance)[1:k]
  as.numeric(names(sorted))
}

r_impute_knn <- function(obj, miss, k) {
  n <- nrow(obj)
  p <- ncol(obj)

  results <- matrix(0, nrow = n, ncol = p)
  miss_rows <- rowSums(miss)
  rows_to_impute <- which(miss_rows > 0)
  ni <- row_index_r(nrow(miss)) # Rows to calculate distance from
  dist_matrix <- matrix(NA, nrow = nrow(ni), ncol = ncol(ni)) # Initialize distant matrix

  for (i in rows_to_impute) {
    # Vector that hold distant, exclude self. More complex but avoid edge case where
    # distant = 0 is not self
    working_dist <- dist_matrix[i, ]
    # neighbor row index, skip self row.
    neigh_index <- ni[i, ]
    # Edge case, if i = nrow(obj), then bumping the index will create an error because
    names(working_dist) <- neigh_index
    ## KNN step 1: calc distant between a row and all rows (self included)
    for (j in seq_along(neigh_index)) {
      working_dist[j] <- calc_distant(
        obj[i, ],
        obj[neigh_index[j], ],
        miss[i, ],
        miss[neigh_index[j], ]
      )
    }

    ## KNN step 2: get the nearest neighbor. nn is the index of rows of nn of obj
    nn <- knn_partial_sort(working_dist, k = k)

    ## KNN step 3: mean impute
    ## Travel along the columns of the rows being currently imputed
    for (z in seq_along(miss[i, ])) {
      # If a value is missing for row i
      if (miss[i, z]) {
        # Calculate total number of valid neighbor. This is the denominator
        valid <- miss[nn, z] == 0
        n_valid <- sum(valid)
        if (sum(valid) == 0) {
          # Failure to impute will be stored as NA. Use with RCpp
          results[i, z] <- NA
        }
        # store imputed results in the results matrix to avoid affecting other imputation
        results[i, z] <- sum(obj[nn, z][valid]) / n_valid
      }
    }
  }
  # obj is pre-imputed data with missing values = 0. results are all zeros, with -1 being
  # failure and non zero being imputed. adding obj and results will give a matrix of
  # imputed values being non zero entry in results while -1 being failure to impute
  # for post processing
  return(obj + results)
}
