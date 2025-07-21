draft <- function() {
  set.seed(1234)
  m <- matrix(sample(1:100, size = 10 * 5, replace = TRUE), ncol = 5)
  miss <- sample(1:length(m), size = 10)
  m[miss] <- NA
  m[1, 3:5] <- NA
  colnames(m) <- paste0("s", 1:5)
  row.names(m) <- paste0("r", 1:10)

  m

  r_impute <- function(obj, colmax = 1, rowmax = 0.5) {
    # Pre-processing
    pn <- dim(obj)
    miss <- is.na(obj)
    cmiss <- colSums(miss) / pn[1]
    if (any(cmiss >= colmax)) {
      stop("Column(s) missing exceeded colmax")
    }
    rmiss <- rowSums(miss) / pn[2]
    mean_imp_rows <- obj[rmiss >= rowmax, ]
    pre_imp_rows <- obj[rmiss < rowmax, ]
    pre_imp_miss <- is.na(pre_imp_rows)
    if (any(colSums(pre_imp_miss) / nrow(pre_imp_miss) >= colmax)) {
      stop("Column(s) missing exceeded colmax after filtering. Consider removing column(s) with too high %NA")
    }
    # Impute
    pre_imp_miss[pre_imp_miss] <- 0
    post_imp <- r_impute_knn(pre_imp_rows, pre_imp_miss)
    # Post-processing
    return()
  }

  # Place holder
  r_impute_knn <- function(obj, miss) {

  }

  # Internals of r_impute_knn
  ## Initializations
  results <- matrix(0, nrow = nrow(obj), ncol(obj))
  miss_rows <- rowSums(miss)
  rows_to_impute <- which(miss_rows > 0)
  ## Rows to be knn imputed
  obj <- m[3:nrow(m), ]
  miss <- is.na(obj)
  results <- matrix(0, nrow = nrow(obj), ncol(obj))
  miss_rows <- rowSums(miss)
  rows_to_impute <- which(miss_rows > 0)

  # For i in rows_to_impute. Example slice i = 3
  i <- 3
  working_dist <- rep(NA, length = nrow(obj))
  for (j in 1:nrow(obj)) {
    working_dist[j] <- calc_distant(obj[i, ], obj[j, ], miss[i, ], miss[j, ])
  }

  ## KNN step 1: calc distant between a row and all rows (self included)
  calc_distant <- function(v1, v2, m1, m2) {
    valid <- (m1 + m2) == 0
    total_valid <- sum(valid)
    if (total_valid == 0) {
      return(NA)
    }
    dist <- sum((v1[valid] - v2[valid])^2)
    return(dist / total_valid)
  }

  ## KNN step 2: get the nearest neighbor. Just placeholder for the implementation in C++ with partial sort,
  ## return index of rows that are closest neighbors
  knn_partial_sort <- function(distance, k) {
    names(distance) <- 1:length(distance)
    sorted <- sort(distance)[1:(k + 1)]
    as.numeric(names(sorted))
  }

  ## KNN step 3: mean impute
  nn <- knn_partial_sort(working_dist, k = 3)
  neighbor <- obj[nn[2:length(nn)], ]
  m_neighbor <- miss[nn[2:length(nn)], ]

  for (z in seq_along(miss[i, ])) {
    if (miss[i, ][z]) {
      valid <- m_neighbor[, z] == TRUE
      n_valid <- sum(valid)
      if (sum(valid) == 0) {
        # Failure to impute will be stored as 2 in the missing matrix.
        miss[i, ][z] <- 2
      }
      # store imputed results in the results matrix
      results[i, z] <- sum(neighbor[, z]) / n_valid
    }
  }
}
