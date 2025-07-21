sim_mat <- function(n = 100, m = 100, nchr = 2, ngrp = 1, perc_NA = 0.5, perc_col_NA = 0.5) {
  stopifnot(
    n > 1, m > 1, nchr >= 1, nchr <= 22,
    perc_NA >= 0, perc_NA <= 1, perc_col_NA >= 0, perc_col_NA <= 1
  )

  # Create and scale the matrix to between 0 and 1 per column
  d_length <- n * m
  d <- matrix(stats::rnorm(d_length), nrow = n, ncol = m)
  mins <- matrixStats::colMins(d)
  maxs <- matrixStats::colMaxs(d)
  ranges <- maxs - mins
  d <- sweep(d, 2, mins, "-")
  d <- sweep(d, 2, ranges, "/")

  # Generate realistic feature and sample names
  feature <- seq_len(n)
  chr <- sample(paste0("chr", seq_len(nchr)), size = n, replace = TRUE)
  group_feature <- data.frame(feature_id = paste0("feat", feature), group = chr)
  colnames(d) <- paste0("s", seq_len(m))
  grp <- sample(paste0("grp", seq_len(ngrp)), size = m, replace = TRUE)
  group_sample <- data.frame(sample_id = colnames(d), group = grp)
  rownames(d) <- group_feature$feature_id

  # Introduce missing values in selected features (rows)
  feature_miss_size <- floor(perc_col_NA * n)
  na_size <- floor(perc_NA * m)

  if (feature_miss_size > 0 && na_size > 0) {
    feature_miss <- sample.int(n, size = feature_miss_size)
    col_idx <- c(replicate(feature_miss_size, sample.int(m, na_size)))
    row_idx <- rep(feature_miss, each = na_size)
    d[cbind(row_idx, col_idx)] <- NA
  }

  return(list(input = d, group_feature = group_feature, group_sample = group_sample))
}
