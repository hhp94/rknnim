// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>
#include <Rcpp.h>
#include <vector>
#include <algorithm>

using namespace Rcpp;

// [[Rcpp::export]]
arma::uvec col_index(const arma::uword i, const arma::uword p) {
  // The last index of the initial sequence is p-2 because self is excluded.
  arma::uword end = p - 2;
  arma::uvec v = arma::regspace<arma::uvec>(arma::uword(0), end);

  if (i <= end) {
    v.subvec(i, end) += 1;
  }

  return v;
}

// [[Rcpp::export]]
double calc_distant_arma(
    const arma::vec& v1,
    const arma::vec& v2,
    const arma::uvec& m1,
    const arma::uvec& m2
) {
  double dist = 0.0;
  arma::uword total_valid = 0;
  for (arma::uword i = 0; i < v1.n_elem; ++i) {
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

// [[Rcpp::export]]
arma::uvec find_knn_indices_arma(const arma::vec& distances, const arma::uword k) {
  // Create a vector of indices: 0, 1, ..., n-1
  arma::uvec indices = arma::regspace<arma::uvec>(arma::uword(0), distances.n_elem - 1);

  std::nth_element(
    indices.begin(),
    indices.begin() + k,
    indices.end(),
    [&](arma::uword a, arma::uword b) {
      return distances(a) < distances(b);
    }
  );

  // Return the first k indices.
  return indices.head(k);
}

//' @export
// [[Rcpp::export]]
arma::mat impute_knn_arma(const arma::mat& obj, const arma::umat& miss, const arma::uword k) {
  const arma::uword n = obj.n_rows;
  const arma::uword p = obj.n_cols;
  const arma::uword n_neighbors = p - 1;

  arma::mat imputed = obj;

  // Loop over all columns to impute missing values.
  for(arma::uword i = 0; i < p; i++) {
    const arma::uword n_miss = arma::sum(miss.col(i));
    if (n_miss == 0) {
      continue;
    }

    arma::uvec neighbor_indices = col_index(i, p);
    arma::vec working_dist(n_neighbors);

    // Calculate the distance from column 'i' to all other columns.
    for(arma::uword j = 0; j < n_neighbors; j++) {
      working_dist(j) = calc_distant_arma(
        obj.col(i),
        obj.col(neighbor_indices(j)),
        miss.col(i),
        miss.col(neighbor_indices(j))
      );
    }

    // Find the k-nearest neighbors.
    arma::uvec sorted_indices = find_knn_indices_arma(working_dist, k);
    arma::uvec nn_indices = neighbor_indices.elem(sorted_indices);

    // Impute each missing value in column 'i'.
    for(arma::uword r = 0; r < n; r++) {
      if (miss(r, i) == 1) {
        double sum_val = 0.0;
        arma::uword count = 0;
        for (arma::uword kk = 0; kk < k; kk++) {
          if (miss(r, nn_indices(kk)) == 0) {
            sum_val += obj(r, nn_indices(kk));
            count++;
          }
        }

        if (count > 0) {
          imputed(r, i) = sum_val / count;
        } else {
          imputed(r, i) = arma::datum::nan;
        }
      }
    }
  }

  return imputed;
}
