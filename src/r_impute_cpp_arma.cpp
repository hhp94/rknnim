// [[Rcpp::depends(RcppArmadillo)]]

#include <RcppArmadillo.h>
#include <vector>
#include <algorithm>
#include <stdexcept>
#include <cmath>

using namespace Rcpp;

// Helper function for Squared Euclidean (method 0)
// [[Rcpp::export]]
double calc_distance_sqeuclid(
    const arma::mat &obj,
    const arma::umat &miss,
    const arma::uword idx1,
    const arma::uword idx2,
    const double total_rows)
{
  const double *col1_ptr = obj.colptr(idx1);
  const double *col2_ptr = obj.colptr(idx2);
  const arma::uword *miss1_ptr = miss.colptr(idx1);
  const arma::uword *miss2_ptr = miss.colptr(idx2);

  double dist = 0.0;
  arma::uword n_valid = 0;

  for (arma::uword r = 0; r < obj.n_rows; ++r) {
    if (miss1_ptr[r] == 0 && miss2_ptr[r] == 0) {
      double diff = col1_ptr[r] - col2_ptr[r];
      dist += diff * diff;
      ++n_valid;
    }
  }

  if (n_valid == 0) {
    return arma::datum::inf;
  }

  const double nn = static_cast<double>(n_valid);
  dist *= (total_rows / nn);

  return dist;
}

// Helper function for Manhattan (method 1)
double calc_distance_manhattan(
    const arma::mat &obj,
    const arma::umat &miss,
    const arma::uword idx1,
    const arma::uword idx2,
    const double total_rows)
{
  const double *col1_ptr = obj.colptr(idx1);
  const double *col2_ptr = obj.colptr(idx2);
  const arma::uword *miss1_ptr = miss.colptr(idx1);
  const arma::uword *miss2_ptr = miss.colptr(idx2);

  double dist = 0.0;
  arma::uword n_valid = 0;

  for (arma::uword r = 0; r < obj.n_rows; ++r) {
    if (miss1_ptr[r] == 0 && miss2_ptr[r] == 0) {
      double diff = col1_ptr[r] - col2_ptr[r];
      dist += std::abs(diff);
      ++n_valid;
    }
  }

  if (n_valid == 0) {
    return arma::datum::inf;
  }

  const double nn = static_cast<double>(n_valid);
  dist *= (total_rows / nn);

  return dist;
}

// Helper function for impute.knn (method 2)
double calc_distance_knn(
    const arma::mat &obj,
    const arma::umat &miss,
    const arma::uword idx1,
    const arma::uword idx2,
    const double total_rows)
{
  const double *col1_ptr = obj.colptr(idx1);
  const double *col2_ptr = obj.colptr(idx2);
  const arma::uword *miss1_ptr = miss.colptr(idx1);
  const arma::uword *miss2_ptr = miss.colptr(idx2);

  double dist = 0.0;
  arma::uword n_valid = 0;

  for (arma::uword r = 0; r < obj.n_rows; ++r) {
    if (miss1_ptr[r] == 0 && miss2_ptr[r] == 0) {
      double diff = col1_ptr[r] - col2_ptr[r];
      dist += diff * diff;
      ++n_valid;
    }
  }

  if (n_valid == 0) {
    return arma::datum::inf;
  }

  const double nn = static_cast<double>(n_valid);
  dist /= nn;

  return dist;
}

// [[Rcpp::export]]
arma::mat distance_matrix(
    const arma::mat &obj,
    const arma::umat &miss,
    const arma::uvec &index_miss,
    const arma::uvec &index_not_miss,
    const int method)
{
  using dist_func_t = double (*)(
      const arma::mat &,
      const arma::umat &,
      const arma::uword,
      const arma::uword,
      const double);

  dist_func_t calc_dist = nullptr;

  switch (method) {
  case 0:
    calc_dist = calc_distance_sqeuclid;
    break;
  case 1:
    calc_dist = calc_distance_manhattan;
    break;
  case 2:
    calc_dist = calc_distance_knn;
    break;
  default:
    throw std::invalid_argument("0 = Euclid; 1 = Manhattan; 2 = impute.knn.");
  }

  arma::mat dist_mat(obj.n_cols, index_miss.n_elem);
  dist_mat.fill(arma::datum::inf);
  const double total_rows = static_cast<double>(obj.n_rows);

  // Compute the lower triangular part of the distance matrix. To optimize computation
  // and avoid rearranging the matrix, columns are ordered with missing values first,
  // using index_miss. Thus, the matrix columns correspond to index_miss (columns with
  // missing values), while the rows are ordered as index_miss followed by columns not in
  // index_miss (i.e., index_not_miss). This ordering is critical because the concatenated
  // rows of the distance matrix represent the neighbor indices used for finding neighbors.
  // For example, if there are 4 columns and columns 1 and 3 have missing values, the
  // distance matrix columns will be [1, 3], and the rows will be [1, 3, 0, 2]. This
  // structure differs from the (i, j) indices, which track the specific row and column
  // positions in the matrix. Just write it out using pencil and paper and you'll get
  // it. The diagonal elements in the top triangular submatrix are automatically
  // set to Inf.
  for (arma::uword i = 1; i < index_miss.n_elem; i++)
  {
    for (arma::uword j = 0; j < i; j++)
    {
      double dist = calc_dist(obj, miss, index_miss(i), index_miss(j), total_rows);
      dist_mat(i, j) = dist;
      dist_mat(j, i) = dist; // Symmetric assignment for the upper triangular part
    }
  }

  // Compute the remaining part of the matrix (distances between non-missing columns and
  // missing columns). If all columns have missing values (i.e., the matrix is square with
  // obj.n_cols == index_miss.n_elem), this loop is automatically skipped.
  for (arma::uword i = index_miss.n_elem; i < obj.n_cols; i++)
  {
    for (arma::uword j = 0; j < index_miss.n_elem; j++)
    {
      // Offset i to index into index_not_miss from 0
      double dist = calc_dist(obj, miss, index_not_miss(i - index_miss.n_elem), index_miss(j), total_rows);
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

//' @export
// [[Rcpp::export]]
arma::mat impute_knn_arma(const arma::mat &obj, const arma::umat &miss, const arma::uword k, const int method)
{
  // Only compute distance of rows with any missing and do it once
  arma::urowvec any_missing_r = arma::any(miss, 0);
  arma::uvec col_index_miss = arma::find(any_missing_r);
  if (col_index_miss.n_elem == 0)
  {
    return obj;
  }
  arma::mat imputed = obj;
  arma::uvec col_index_non_miss = arma::find(any_missing_r == 0);
  // Longest step. dist_mat holds distance from missing column to all other columns
  const arma::mat dist_mat = distance_matrix(obj, miss, col_index_miss, col_index_non_miss, method);
  // dist_mat row index - neighbors - are arranged as col_index_miss then col_index_non_miss
  // while columns are col_index_miss
  // return dist_mat;
  arma::uvec neighbor_index = arma::join_vert(col_index_miss, col_index_non_miss);
  // Loop through only missing columns
  for (arma::uword i = 0; i < col_index_miss.n_elem; i++)
  {
    arma::vec dist_vec = dist_mat.col(i);
    arma::uvec candidate_indices = arma::find_finite(dist_vec);
    if (candidate_indices.n_elem == 0)
    {
      continue;
    }
    // Caps `k` to the number of available neighbors.
    arma::uword n_neighbors = std::min(k, candidate_indices.n_elem);
    // index of nn relative to dist_vec rows
    arma::uvec nn_dist_indices = find_knn_indices_arma(dist_vec, n_neighbors);
    // actual index of nn columns relative to obj columns
    arma::uvec nn_columns = neighbor_index.elem(nn_dist_indices);

    // now we move on to nn imputation. Find the rows to impute for this column
    arma::uword target_col_idx = col_index_miss(i);
    arma::uvec rows_to_impute = arma::find(miss.col(target_col_idx));

    // man, range-for-loop is trippy. kind of not used to this but it works
    for (arma::uword row_idx : rows_to_impute)
    {
      double sum = 0.0;
      arma::uword count = 0;
      for (arma::uword neighbor_col_idx : nn_columns)
      {
        if (miss(row_idx, neighbor_col_idx) == 0)
        {
          sum += obj(row_idx, neighbor_col_idx);
          count++;
        }
      }
      if (count > 0)
      {
        imputed(row_idx, target_col_idx) = sum / count;
      }
    }
  }

  return imputed;
}
