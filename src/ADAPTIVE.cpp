#include <Rcpp.h>
using namespace Rcpp;

// [[Rcpp::export]]
NumericVector adaptive_C(NumericVector y, NumericVector betas, double u, double c, double G = 10) {
  int N = y.size();
  NumericVector var(N);
  var[0] = u;

  for (int i = 1; i < N; ++i) {
    var[i] = var[i - 1] + betas[0] * (c - 1 / (1 + std::exp(G * (y[i - 1] - var[i - 1]))));
  }

  return var;
}
