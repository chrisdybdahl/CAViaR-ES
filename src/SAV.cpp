#include <Rcpp.h>
using namespace Rcpp;

// [[Rcpp::export]]
NumericVector SAV_C(NumericVector y, NumericVector betas, double u, double c) {
  int N = y.size();
  NumericVector var(N);
  var[0] = u;

  for (int i = 1; i < N; ++i) {
    var[i] = betas[0] + betas[1] * var[i - 1] + betas[2] * std::abs(y[i - 1]);
  }

  return var;
}