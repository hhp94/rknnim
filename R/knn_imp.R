#' K-Nearest Neighbor (k-NN) Imputation
#'
#' @description
#' This function imputes missing values in a numeric matrix using the k-Nearest
#' Neighbors algorithm. It follows a two-stage process. First, it imputes
#' columns with a proportion of missing values below `colmax` using k-NN.
#' Second, any remaining missing values (including those in columns that
#' exceeded `colmax`) are imputed using the column mean.
#'
#' @details
#' The distance calculation between columns for identifying nearest neighbors is
#' scaled based on the number of non-missing value pairs. Specifically, the
#' raw distance is penalized by scaling it up for columns that have fewer
#' overlapping observations. This penalizes distances for columns with very few
#' shared observations used for distance calculations. The
#' `impute.knn` method averages the distances over the number of matching positions,
#' so a column with only one matching value to calculate distance from might have a lower
#' raw distance than a column with many matched values. See also [stats::dist()]
#'
#' @param obj A numeric matrix with samples in rows and features in columns.
#' @param k An integer specifying the number of nearest neighbors to use for
#'   imputation. Must be between 1 and the number of columns.
#' @param colmax A numeric value between 0 and 1. This is the threshold for the
#'   proportion of missing values in a column. Columns exceeding this
#'   threshold will be imputed using the column mean instead of k-NN.
#' @param rowmax A numeric value between 0 and 1. This is the maximum
#'   allowable proportion of missing values in any single row. If a row
#'   exceeds this threshold, the function will stop with an error.
#' @param method A character string specifying the distance metric for k-NN.
#'   Acceptable values are `"euclidean"`, `"manhattan"`, or `"impute.knn"`.
#'   Defaults to `"euclidean"`. See details
#'
#' @return A numeric matrix of the same dimensions as `obj` with missing values
#'   imputed.
#'
#' @export
#'
#' @examples
#' # See ?khanmiss1
#' data(khanmiss1)
#' sum(is.na(khanmiss1))
#'
#' # Perform k-NN imputation. `khanmiss1` stores genes in the row so we have to t().
#' # set method to "impute.knn" to mimic how distant is scaled in impute::impute.knn.
#' imputed <- knn_imp(obj = t(khanmiss1), k = 3, colmax = 0.5, rowmax = 0.8, method = "euclidean")
#' imputed[1:5, 1:20]
#' sum(is.na(imputed))
knn_imp <- function(obj, k, colmax = 0.5, rowmax = 0.8, method = c("euclidean", "manhattan", "impute.knn")) {
  # Pre-conditioning
  method <- match.arg(method)
  checkmate::assert_matrix(obj, mode = "numeric", min.rows = 1, min.cols = 2, col.names = "unique")
  checkmate::assert_integerish(k, lower = 1, upper = ncol(obj), len = 1)
  checkmate::assert_numeric(colmax, lower = 0, upper = 1, len = 1)
  checkmate::assert_numeric(rowmax, lower = 0, upper = 1, len = 1)

  if (ncol(obj) < nrow(obj)) {
    warning("we expect features to be on the columns")
  }

  miss <- is.na(obj)
  if (sum(miss) == 0) {
    message("No missing data")
    return(obj)
  }

  rmiss <- rowSums(miss) / ncol(obj)
  if (any(rmiss >= rowmax)) {
    stop("Row(s) missing exceeded rowmax. Consider removing rows(s) with too high NA %")
  }

  cmiss <- colSums(miss)
  mean_imp <- cmiss / nrow(obj) >= colmax
  if (any(mean_imp == 1)) {
    stop("Col(s) with all missing detected. Consider removing before proceed")
  }

  # Partition
  mean_imp_cols <- obj[, mean_imp, drop = FALSE]
  pre_imp_cols <- obj[, !mean_imp, drop = FALSE]
  pre_imp_miss <- miss[, !mean_imp, drop = FALSE]
  p_pre <- ncol(pre_imp_cols)

  if (p_pre > 0 && any(rowSums(pre_imp_miss) / p_pre >= rowmax)) {
    stop("Row(s) missing exceeded rowmax. Consider removing rows(s) with too high NA %")
  }

  # Impute
  post_imp <- impute_knn_naive(
    obj = pre_imp_cols,
    miss = pre_imp_miss,
    k = k,
    n_col_miss = cmiss,
    method = switch(method, "euclidean" = 0L, "manhattan" = 1L, "impute.knn" = 2L)
  )
  colnames(post_imp) <- colnames(pre_imp_cols)

  # Final
  ## Reconstruct and impute values still missing
  final <- obj
  final[, !mean_imp] <- post_imp
  if (anyNA(final)) {
    na_indices <- which(is.na(final), arr.ind = TRUE)
    col_means <- colMeans(obj, na.rm = TRUE)
    final[na_indices] <- col_means[na_indices[, 2]]
  }
  return(final)
}
