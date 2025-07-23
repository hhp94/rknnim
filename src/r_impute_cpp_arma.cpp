// [[Rcpp::depends(RcppArmadillo)]]

#include <RcppArmadillo.h>
#include <vector>
#include <algorithm>

using namespace Rcpp;

arma::mat impute_knn_arma(const arma::mat& obj, const arma::umat& miss, const arma::uword k);

// Helper function to calculate the distance matrix between columns.
// The distance is the mean squared difference, calculated only on non-missing pairs.
// [[Rcpp::export]]
arma::mat distant_matrix(const arma::mat& obj, const arma::umat& miss) {
  arma::uword n_cols = obj.n_cols;
  arma::uword n_rows = obj.n_rows;

  // Initialize the final distance matrix with Inf. This serves two purposes:
  // 1. We can easily find valid (finite) distances later.
  // 2. A column will not be selected as its own nearest neighbor.
  arma::mat dist_mat(n_cols, n_cols);
  dist_mat.fill(arma::datum::inf);

  for (arma::uword i = 1; i < n_cols; ++i) {
    for (arma::uword j = 0; j < i; ++j) {
      // Use column pointers for efficient access.
      const double* obj_col_i_ptr = obj.colptr(i);
      const double* obj_col_j_ptr = obj.colptr(j);
      const arma::uword* miss_col_i_ptr = miss.colptr(i);
      const arma::uword* miss_col_j_ptr = miss.colptr(j);

      double sum_sq_diff = 0.0;
      int valid_pairs = 0;

      // Iterate over all rows to find common non-missing values.
      for (arma::uword r = 0; r < n_rows; ++r) {
        // Check if the value is present in BOTH columns.
        if (miss_col_i_ptr[r] == 0 && miss_col_j_ptr[r] == 0) {
          double diff = obj_col_i_ptr[r] - obj_col_j_ptr[r];
          sum_sq_diff += diff * diff;
          valid_pairs++;
        }
      }

      // If there are any valid pairs, calculate the mean squared difference.
      // Otherwise, the distance remains Inf.
      if (valid_pairs > 0) {
        double distance = sum_sq_diff / valid_pairs;
        dist_mat(i, j) = distance;
        dist_mat(j, i) = distance; // The distance matrix is symmetric.
      }
    }
  }

  return dist_mat;
}

// Helper function to find the indices of the k smallest values in a vector.
// [[Rcpp::export]]
arma::uvec find_knn_indices_arma(const arma::vec& distances, arma::uword k) {
  arma::uvec indices = arma::regspace<arma::uvec>(0, distances.n_elem - 1);

  // Partition the indices vector such that the k-th smallest element is in its
  // head k position. Elements before it are smaller or equal.
  std::nth_element(
    indices.begin(),
    // nth_element is 0-indexed, so we use k-1 to get the k-th element
    indices.begin() + (k - 1),
    indices.end(),
    // Custom comparator lambda function that compares the distances at the given indices.
    [&](arma::uword a, arma::uword b) {
      return distances(a) < distances(b);
    }
  );

  // Return the first k indices, which now correspond to the k smallest distances.
  return indices.head(k);
}

//' @export
// [[Rcpp::export]]
arma::mat impute_knn_arma(const arma::mat& obj, const arma::umat& miss, const arma::uword k) {
  const arma::uword n_rows = obj.n_rows;
  const arma::uword n_cols = obj.n_cols;
  arma::mat imputed = obj; // Start with a copy of the original matrix.

  // 1. Calculate the distance matrix between all columns ONCE.
  arma::mat dist_mat = distant_matrix(obj, miss);

  // 2. Iterate through each column to check for and impute missing values.
  for (arma::uword i = 0; i < n_cols; i++) {
    // If this column has no missing values, skip it
    if (arma::sum(miss.col(i)) == 0) {
      continue;
    }
    // 3. Find the nearest neighbors for the current column 'i'.
    // Get the vector of distances from column 'i' to all other columns
    arma::vec d_vec = dist_mat.col(i);
    // Get indices for sorting
    arma::uvec candidate_indices = arma::find_finite(d_vec);

    // If no neighbors with finite distance exist, we cannot impute this column.
    if (candidate_indices.n_elem == 0) {
        // As a fallback, explicitly set any missing values in this column to NaN.
        for(arma::uword r = 0; r < n_rows; ++r) {
            if(miss(r, i) == 1) {
                imputed(r, i) = arma::datum::nan;
            }
        }
        continue;
    }
    // Create a sub-vector containing only the finite distances.
    arma::vec finite_dists = d_vec.elem(candidate_indices);
    // Determine the actual number of neighbors to use (can't be more than available).
    arma::uword n_neighbors = std::min(k, finite_dists.n_elem);

    // Find the indices of the `n_neighbors` smallest distances within `finite_dists`.
    // Note: These are indices *relative* to the `finite_dists` vector.
    arma::uvec knn_relative_indices = find_knn_indices_arma(finite_dists, n_neighbors);

    // Map these relative indices back to the original matrix column indices.
    arma::uvec nn_col_indices = candidate_indices.elem(knn_relative_indices);

    // 4. Impute missing values in column 'i' using the identified neighbors.
    for (arma::uword r = 0; r < n_rows; r++) {
      // Only impute if the value at this position is marked as missing.
      if (miss(r, i) == 1) {
        double imputation_sum = 0.0;
        arma::uword valid_neighbor_count = 0;

        // Iterate through the neighbor columns.
        for (arma::uword neighbor_idx : nn_col_indices) {
          // Check if the neighbor has a non-missing value at this specific row.
          if (miss(r, neighbor_idx) == 0) {
            imputation_sum += obj(r, neighbor_idx);
            valid_neighbor_count++;
          }
        }

        // If we found at least one valid neighbor value, calculate the mean.
        if (valid_neighbor_count > 0) {
          imputed(r, i) = imputation_sum / valid_neighbor_count;
        } else {
          // Otherwise, no imputation ispossible for this specific value.
          // Set it to NaN. The user can decide how to handle NaNs later.
          imputed(r, i) = arma::datum::nan;
        }
      }
    }
  }
  return imputed;
}
