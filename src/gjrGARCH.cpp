#include <Rcpp.h>
using namespace Rcpp;

// [[Rcpp::export]]
NumericVector gjrGARCH_C(NumericVector y, NumericVector betas, double u, double c) {
  int N = y.size();
  NumericVector var(N);
  var[0] = u;

  for (int i = 1; i < N; ++i) {
    var[i] = betas[0] + betas[1] * var[i - 1] + betas[2] * y[i - 1] + betas[3] * y[i - 1] * (y[i - 1] < 0);
  }

  return var;
}