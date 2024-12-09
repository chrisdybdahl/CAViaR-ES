#include <Rcpp.h>
using namespace Rcpp;

// [[Rcpp::export]]
NumericVector al_log_loss_function_C(NumericVector y, NumericVector Q, NumericVector ES, double c) {
  int N = y.size();
  NumericVector loss_vector(N);

  for (int i = 0; i < N; ++i) {
    loss_vector[i] = -std::log((c - 1) / ES[i]) -
      ((y[i] - Q[i]) * (c - (y[i] <= Q[i]))) / (c * ES[i]);
  }

  return loss_vector;
}
