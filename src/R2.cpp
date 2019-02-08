// Implements row and column sums of squared elements of matrices. Only one
//   of these functions is ever called, and it is only called once (during
//   initialization), so simplicity is better than performance.

#include <Rcpp.h>
using namespace Rcpp;

// [[Rcpp::export]]
NumericVector calc_R2_colsums_rcpp (NumericMatrix X) {
  int n = X.nrow();
  int p = X.ncol();

  NumericVector out(p);

  for (int j = 0; j < p; j++) {
    for (int i = 0; i < n; i++) {
      out[j] += X(i, j) * X(i, j);
    }
  }

  return out;
}

// [[Rcpp::export]]
NumericVector calc_R2_rowsums_rcpp (NumericMatrix X) {
  int n = X.nrow();
  int p = X.ncol();

  NumericVector out(n);

  for (int j = 0; j < p; j++) {
    for (int i = 0; i < n; i++) {
      out[i] += X(i, j) * X(i, j);
    }
  }

  return out;
}
