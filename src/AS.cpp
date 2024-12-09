#include <Rcpp.h>
using namespace Rcpp;

// [[Rcpp::export]]
NumericVector AS_C(NumericVector y, NumericVector betas, double u, double c) {
  int N = y.size();
  NumericVector var(N);
  var[0] = u;

  for (int i = 1; i < N; ++i) {
    double max_val = std::max(y[i - 1], 0.0);
    double min_val = std::min(y[i - 1], 0.0);
    var[i] = betas[0] + betas[1] * var[i - 1] + betas[2] * max_val + betas[3] * min_val;
  }

  return var;
}