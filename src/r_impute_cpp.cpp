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

// [[Rcpp::export]]
NumericMatrix impute_knn_cpp(NumericMatrix obj, LogicalMatrix miss, int k) {
  int n = obj.nrow();
  int p = obj.ncol();
  int nr = n - 1; // n-1 rows to check for nearest neighbors
  NumericMatrix results(n, p);
  IntegerMatrix dist_idx = row_index(n);
  NumericMatrix dist_matrix(n, nr);
  IntegerMatrix nearest_neighbors(n, k);

  for (int i = 0; i < n; i++) {
    int miss_rows = 0;
    for (int j = 0; j < p; j++) {
      miss_rows += miss(i, j);
    }
    if (miss_rows == 0) {
      continue;
    }
    // KNN step 1: calc distance between a row and all other rows (excluding self)
    for (int j = 0; j < nr; j++) {
      dist_matrix(i, j) = calc_distant_cpp(
        obj(i, _),
        obj(dist_idx(i, j), _),
        miss(i, _),
        miss(dist_idx(i, j), _)
      );
    }
    // KNN step 2: get indices of top k neighbors (indices into dist_matrix row)
    nearest_neighbors(i, _) = partial_sort_cpp(dist_matrix(i, _), k);
    // KNN step 3: mean impute for each missing column
    for (int j = 0; j < p; j++) {
      if (miss(i, j) == true) {
        double sum_val = 0.0;
        int count = 0;
        for (int kk = 0; kk < k; kk++) {
          int pos = nearest_neighbors(i, kk);
          int nn_row = dist_idx(i, pos);
          if (miss(nn_row, j) == false) {
            sum_val += obj(nn_row, j);
            count++;
          }
        }
        if (count == 0) {
          results(i, j) = NA_REAL;
        } else {
          results(i, j) = sum_val / count;
        }
      }
    }
  }

  // Compute imputed matrix: obj + results (handling NA properly)
  NumericMatrix imputed = clone(obj);
  for (int i = 0; i < n; i++) {
    for (int j = 0; j < p; j++) {
      imputed(i, j) += results(i, j);
    }
  }

  return imputed;
}
