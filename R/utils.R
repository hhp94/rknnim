sim_mat <- function(
    n = 100, m = 100, nchr = 2, ngrp = 1, perc_NA = 0.5, perc_col_NA = 0.5, transpose = TRUE) {
  stopifnot(
    n > 1, m > 1, nchr >= 1, nchr <= 22, nchr < n,
    perc_NA >= 0, perc_NA <= 1, perc_col_NA >= 0, perc_col_NA <= 1
  )

  # create and scale the matrix to between 0 and 1
  d_length <- n * m
  d <- matrix(stats::rnorm(d_length), nrow = n, ncol = m)
  mins <- matrixStats::colMins(d)
  maxs <- matrixStats::colMaxs(d)
  ranges <- maxs - mins
  d <- sweep(d, 2, mins, "-")
  d <- sweep(d, 2, ranges, "/")

  # generate realistic feature and sample names
  feature <- seq_len(n)
  chr <- sample(paste0("chr", seq_len(nchr)), size = n, replace = TRUE)
  group_feature <- data.frame(feature_id = paste0("feat", feature), group = chr)
  colnames(d) <- paste0("s", seq_len(m))
  grp <- sample(paste0("grp", seq_len(ngrp)), size = m, replace = TRUE)
  group_sample <- data.frame(sample_id = colnames(d), group = grp)
  rownames(d) <- group_feature$feature_id
  d <- t(d)

  # d is t()ed, so n and m are swapped
  col_miss_size <- max(floor(perc_col_NA * n), 0)
  NA_size <- max(floor(perc_NA * m), 0)

  if (col_miss_size > 0 && NA_size > 0) {
    col_miss <- sample.int(n, size = col_miss_size)
    row_idx <- c(replicate(col_miss_size, sample.int(m, NA_size)))
    col_idx <- rep(col_miss, each = NA_size)
    d[cbind(row_idx, col_idx)] <- NA
  }

  if (transpose) {
    d <- t(d)
  }
  return(list(input = d, group_feature = group_feature, group_sample = group_sample))
}
