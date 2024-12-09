#include <Rcpp.h>
using namespace Rcpp;

// [[Rcpp::export]]
NumericVector indirectGARCH_C(NumericVector y, NumericVector betas, double u, double c) {
  int N = y.size();
  NumericVector var(N);
  var[0] = u;

  for (int i = 1; i < N; ++i) {
    var[i] = (1 - 2 * (c < 0.5)) *
             std::pow(betas[0] + betas[1] * std::pow(var[i - 1], 2) + betas[2] * std::pow(y[i - 1], 2), 0.5);
  }

  return var;
}
