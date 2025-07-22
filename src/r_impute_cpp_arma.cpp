#include <RcppArmadillo.h>
#include <Rcpp.h>
// [[Rcpp::depends(RcppArmadillo)]]

arma::imat row_index_arma(int n) {
  // Create matrix of size n x (n-1)
  arma::imat v_mat(n, n - 1);

  // Fill each row with the sequence 0, 1, ..., n-2
  v_mat.each_row() = arma::regspace<arma::irowvec>(0, n - 2);
  // Add 1 where column index is greater than or equal to the row index

  for (int i = 0; i < n; ++i) {
    for (int j = i; j < n - 1; ++j) {
      v_mat(i, j) += 1;
    }
  }

  return v_mat;
}

double calc_distant_arma(
    const arma::vec& v1,
    const arma::vec& v2,
    const arma::uvec& m1,
    const arma::uvec& m2
  ) {
  arma::uword n = v1.n_elem;
  double dist = 0.0;
  arma::uword total_valid = 0;
  for (arma::uword i = 0; i < n; ++i) {
    if (m1[i] == 0 && m2[i] == 0) {
      double diff = v1[i] - v2[i];
      dist += diff * diff;
      total_valid++;
    }
  }
  if (total_valid == 0) {
    return arma::datum::inf;
  }

  return dist / total_valid;

}

arma::uvec partial_sort_arma(const arma::vec& distances, int k) {
  int n = distances.n_elem;
  // Create a vector of indices: 0, 1, ..., n-1
  arma::uvec indices = arma::regspace<arma::uvec>(0, n - 1);

  // Partially sort the indices based on the distance values
  // using a C++ lambda function for the custom comparison.
  std::partial_sort(
    indices.begin(),
    indices.begin() + k,
    indices.end(),
    [&](arma::uword a, arma::uword b) {
      return distances(a) < distances(b);
    }
  );

  // Return the first k sorted indices
  return indices.head(k);
}

//' @export
// [[Rcpp::export]]
arma::mat impute_knn_arma(const arma::mat& obj, const arma::umat& miss, int k) {
  // Initialization
  int n = obj.n_rows;
  int p = obj.n_cols;
  int nr = n - 1;

  arma::mat results(n, p, arma::fill::zeros);
  arma::imat dist_idx = row_index_arma(n);
  arma::mat dist_matrix(n, nr);
  arma::imat nearest_neighbors(n, k);

  // Main loop over each row
  for (int i = 0; i < n; i++) {
    // Check if the row has any missing values by summing the boolean flags.
    if (arma::accu(miss.row(i)) == 0) {
      continue;
    }
    // step 1, Calculate distance between row `i` and all other rows
    for (int j = 0; j < nr; j++) {
      int other_row_idx = dist_idx(i, j);
      dist_matrix(i, j) = calc_distant_arma(
        obj.row(i).t(),
        obj.row(other_row_idx).t(),
        miss.row(i).t(),
        miss.row(other_row_idx).t()
      );
    }

    // step 2, get indices of top k neighbors
    arma::uvec nn_indices_u = partial_sort_arma(dist_matrix.row(i).t(), k);
    nearest_neighbors.row(i) = arma::conv_to<arma::irowvec>::from(nn_indices_u.t());

    // step 3, Mean impute for each missing column, following the original loop structure
    for (int j = 0; j < p; j++) {
      if (miss(i, j) == 1) { // 1 means TRUE for missing
        double sum_val = 0.0;
        int count = 0;
        // Manual loop to calculate mean from valid neighbors
        for (int kk = 0; kk < k; kk++) {
          int pos = nearest_neighbors(i, kk);
          int nn_row = dist_idx(i, pos);
          if (miss(nn_row, j) == 0) { // 0 means FALSE for missing
            sum_val += obj(nn_row, j);
            count++;
          }
        }
        if (count == 0) {
          results(i, j) = arma::datum::nan; // Use NaN for R's NA
        } else {
          results(i, j) = sum_val / count;
        }
      }
    }
  }

  arma::mat imputed = obj;
  for (int i = 0; i < n; i++) {
    for (int j = 0; j < p; j++) {
      if (miss(i, j) == 1) {
        imputed(i, j) = results(i, j);
      }
    }
  }

  return imputed;
}
