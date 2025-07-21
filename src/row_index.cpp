#include <Rcpp.h>
#include <numeric>

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
