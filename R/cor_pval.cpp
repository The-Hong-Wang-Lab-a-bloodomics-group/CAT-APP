#include <Rcpp.h>
#include <cmath>
using namespace Rcpp;

// [[Rcpp::export]]
NumericMatrix fast_t_stats(NumericMatrix cor, int n) {
  int nrow = cor.nrow();
  int ncol = cor.ncol();
  NumericMatrix t_stats(nrow, ncol);
  double sqrt_df = sqrt(n - 2);
  
  for (int i = 0; i < nrow; ++i) {
    for (int j = 0; j < ncol; ++j) {
      double r = cor(i, j);
      double denom = sqrt(1 - r * r);
      t_stats(i, j) = (denom > 1e-10) ? (r * sqrt_df / denom) : NA_REAL;
    }
  }
  return t_stats;
}
