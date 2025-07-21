#include <Rcpp.h>
#include <numeric>
#include <algorithm>
#include <vector>

using namespace Rcpp;

// [[Rcpp::export]]
IntegerMatrix row_index(int n) {
  // Create matrix of size n x (n-1) because self is removed
  IntegerMatrix v_mat(n, n - 1);

  // Fill each row with the sequence 0, 1, ..., n-2 using std::iota
  for (int i = 0; i < n; ++i) {
    std::iota(v_mat.row(i).begin(), v_mat.row(i).end(), 0);
  }

  // Add 1 where column index >= row index
  for (int i = 0; i < n; i++) {
    for (int j = i; j < n - 1; j++) {
      v_mat(i, j) += 1;
    }
  }

  return v_mat;
}

// [[Rcpp::export]]
double calc_distant_cpp(
    NumericVector v1,
    NumericVector v2,
    LogicalVector m1,
    LogicalVector m2
  ) {

  int n = v1.size();
  double dist = 0.0;
  int total_valid = 0;

  for (int i = 0; i < n; ++i) {
    if (!m1[i] && !m2[i]) {
      double diff = v1[i] - v2[i];
      dist += diff * diff;
      total_valid++;
    }
  }

  if (total_valid == 0) {
    return R_PosInf;
  }

  return dist / total_valid;
}

// [[Rcpp::export]]
IntegerVector partial_sort_cpp(NumericVector distances, int k) {
  // initialization
  int n = distances.size();
  std::vector<int> indices(n);
  std::iota(indices.begin(), indices.end(), 0);

  // sort
  std::partial_sort(
    indices.begin(),
    indices.begin() + k,
    indices.end(),
    [&](int a, int b) {
      return distances[a] < distances[b];
    }
  );

  // gather results. Get just the top k neighbor and store as IntegerVector
  IntegerVector result(k);
  for (int i = 0; i < k; ++i) {
    result[i] = indices[i];
  }

  return result;
}

// impute_knn_r <- function(obj, miss, k) {
//   n <- nrow(obj)
//   p <- ncol(obj)
//
//   results <- matrix(0, nrow = n, ncol = p)
//   miss_rows <- rowSums(miss)
//   rows_to_impute <- which(miss_rows > 0)
//   ni <- row_index_r(nrow(miss)) # Rows to calculate distance from
//   dist_matrix <- matrix(NA, nrow = nrow(ni), ncol = ncol(ni)) # Initialize distant matrix

// [[Rcpp::export]]
void impute_knn_cpp(NumericMatrix obj, LogicalMatrix miss, int k) {
  int n = obj.rows();
  int p = obj.cols();

  NumericMatrix results(n, p);
  std::vector<int> rows_to_impute;
  for (int i = 0; i < n; i++) {
    int miss_rows = 0;
    for (int j = 0; j < p; j++) {
      miss_rows += miss(i,j);
    }
    if(miss_rows > 0) {
      rows_to_impute.push_back(i);
    }
  }

  Rcout << "Rows to impute: ";
  for(size_t i = 0; i < rows_to_impute.size(); i++) {
    Rcout << rows_to_impute[i] << " ";
  }
  Rcout << std::endl;
}
