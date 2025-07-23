// [[Rcpp::depends(RcppArmadillo)]]

#include <RcppArmadillo.h>
#include <vector>
#include <algorithm>
#include <cmath>

using namespace Rcpp;

// [[Rcpp::export]]
double calc_distance(
    const arma::mat &obj,
    const arma::umat &miss,
    const arma::uword idx1,
    const arma::uword idx2,
    const double total_rows)
{
  // TRUE only if both values are not missing
  arma::uvec mask = (miss.col(idx1) == 0) % (miss.col(idx2) == 0);
  arma::uvec valid_indices = arma::find(mask);
  const double n_valid = static_cast<double>(valid_indices.n_elem);
  if (n_valid > 0)
  {
    const double *col1_ptr = obj.colptr(idx1);
    const double *col2_ptr = obj.colptr(idx2);
    double dist = 0.0;
    for (arma::uword i = 0; i < n_valid; ++i)
    {
      arma::uword row = valid_indices(i);
      double diff = col1_ptr[row] - col2_ptr[row];
      dist += diff * diff;
    }
    dist = std::sqrt(dist * (total_rows / n_valid));
    return dist;
  }

  return arma::datum::nan;
}

// [[Rcpp::export]]
arma::mat distance_matrix(
    const arma::mat &obj,
    const arma::umat &miss,
    const arma::uvec &index_miss,
    const arma::uvec &index_not_miss)
{
  arma::mat dist_mat(obj.n_cols, index_miss.n_elem);
  dist_mat.fill(arma::datum::nan);
  const double total_rows = static_cast<double>(obj.n_rows);

  // Compute the lower triangular part of the distance matrix. To optimize computation
  // and avoid rearranging the matrix, columns are ordered with missing values first,
  // using index_miss. Thus, the matrix columns correspond to index_miss (columns with
  // missing values), while the rows are ordered as index_miss followed by columns not in
  // index_miss (i.e., index_not_miss). This ordering is critical because the concatenated
  // rows of the distance matrix represent the neighbor indices used for calculations.
  // For example, if there are 4 columns and columns 1 and 3 have missing values, the
  // distance matrix columns will be [1, 3], and the rows will be [1, 3, 0, 2]. This
  // structure differs from the (i, j) indices, which track the specific row and column
  // positions in the matrix. Just write it out using pencil and paper and you'll get
  // it.
  for (arma::uword i = 1; i < index_miss.n_elem; i++)
  {
    for (arma::uword j = 0; j < i; j++)
    {
      double dist = calc_distance(obj, miss, index_miss(i), index_miss(j), total_rows);
      dist_mat(i, j) = dist;
      dist_mat(j, i) = dist;  // Symmetric assignment for the upper triangular part
    }
  }

  // Compute the remaining part of the matrix (distances between non-missing rows and
  // missing columns). The diagonal elements in the top-left submatrix are automatically
  // set to NaN. If all columns have missing values (i.e., the matrix is square with
  // obj.n_cols == index_miss.n_elem), this loop is skipped.
  for (arma::uword i = index_miss.n_elem; i < obj.n_cols; i++)
  {
    for (arma::uword j = 0; j < index_miss.n_elem; j++)
    {
      // Offset i to index into index_not_miss from 0
      double dist = calc_distance(obj, miss, index_not_miss(i - index_miss.n_elem), index_miss(j), total_rows);
      dist_mat(i, j) = dist;
    }
  }

  return dist_mat;
}

// [[Rcpp::export]]
arma::uvec find_knn_indices_arma(const arma::vec &distances, arma::uword k)
{
  arma::uvec indices = arma::regspace<arma::uvec>(0, distances.n_elem - 1);

  // Partition the indices vector such that the k-th smallest element is in its
  // head k position. Elements before it are smaller or equal.
  std::nth_element(
      indices.begin(),
      // nth_element is 0-indexed, so we use k-1 to get the k-th element
      indices.begin() + (k - 1),
      indices.end(),
      // Custom comparator lambda function that compares the distances at the given indices.
      [&](arma::uword a, arma::uword b)
      {
        return distances(a) < distances(b);
      });

  // Return the first k indices, which now correspond to the k smallest distances.
  return indices.head(k);
}

// [[Rcpp::export]]
arma::mat impute_knn_arma(const arma::mat &obj, const arma::umat &miss, const arma::uword k)
{
  const arma::uword n_rows = obj.n_rows;
  const arma::uword n_cols = obj.n_cols;
  // CCCCrucial. find just columns with missing.
  arma::urowvec any_missing_r = arma::any(miss, 0);
  arma::uvec any_missing = arma::conv_to<arma::uvec>::from(any_missing_r);
  arma::uvec col_index_miss = arma::sort(arma::find(any_missing));
  arma::uvec col_index_non_miss = arma::sort(arma::find(any_missing == 0));
  if (col_index_miss.n_elem == 0)
  {
    return obj;
  }
  // Impute just columns with missing. Most expensive part
  arma::mat imputed = obj;
  const arma::mat dist_mat = distance_matrix(obj, miss, col_index_miss, col_index_non_miss);
  return dist_mat;
}
//   arma::mat dist_mat = distance_matrix(obj, miss);
//
//   // 2. Iterate through each column to check for and impute missing values.
//   for (arma::uword i = 0; i < n_cols; i++) {
//     // If this column has no missing values, skip it
//     if (arma::sum(miss.col(i)) == 0) {
//       continue;
//     }
//     // 3. Find the nearest neighbors for the current column 'i'.
//     // Get the vector of distances from column 'i' to all other columns
//     arma::vec d_vec = dist_mat.col(i);
//     // Get indices for sorting
//     arma::uvec candidate_indices = arma::find_finite(d_vec);
//
//     // If no neighbors with finite distance exist, we cannot impute this column.
//     if (candidate_indices.n_elem == 0) {
//         // As a fallback, explicitly set any missing values in this column to NaN.
//         for(arma::uword r = 0; r < n_rows; ++r) {
//             if(miss(r, i) == 1) {
//                 imputed(r, i) = arma::datum::nan;
//             }
//         }
//         continue;
//     }
//     // Create a sub-vector containing only the finite distances.
//     arma::vec finite_dists = d_vec.elem(candidate_indices);
//     // Determine the actual number of neighbors to use (can't be more than available).
//     arma::uword n_neighbors = std::min(k, finite_dists.n_elem);
//
//     // Find the indices of the `n_neighbors` smallest distances within `finite_dists`.
//     // Note: These are indices *relative* to the `finite_dists` vector.
//     arma::uvec knn_relative_indices = find_knn_indices_arma(finite_dists, n_neighbors);
//
//     // Map these relative indices back to the original matrix column indices.
//     arma::uvec nn_col_indices = candidate_indices.elem(knn_relative_indices);
//
//     // 4. Impute missing values in column 'i' using the identified neighbors.
//     for (arma::uword r = 0; r < n_rows; r++) {
//       // Only impute if the value at this position is marked as missing.
//       if (miss(r, i) == 1) {
//         double imputation_sum = 0.0;
//         arma::uword valid_neighbor_count = 0;
//
//         // Iterate through the neighbor columns.
//         for (arma::uword neighbor_idx : nn_col_indices) {
//           // Check if the neighbor has a non-missing value at this specific row.
//           if (miss(r, neighbor_idx) == 0) {
//             imputation_sum += obj(r, neighbor_idx);
//             valid_neighbor_count++;
//           }
//         }
//
//         // If we found at least one valid neighbor value, calculate the mean.
//         if (valid_neighbor_count > 0) {
//           imputed(r, i) = imputation_sum / valid_neighbor_count;
//         } else {
//           // Otherwise, no imputation ispossible for this specific value.
//           // Set it to NaN. The user can decide how to handle NaNs later.
//           imputed(r, i) = arma::datum::nan;
//         }
//       }
//     }
//   }
//   return imputed;
// }
