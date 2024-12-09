#include <Rcpp.h>
using namespace Rcpp;

// [[Rcpp::export]]
NumericVector ar_ES_C(NumericVector y, NumericVector Q, NumericVector gammas, NumericVector x = NumericVector::create()) {
  int n = Q.size();

  // Initialize x if it is not provided
  if (x.size() == 0) {
    x = NumericVector(n);
    x[0] = gammas[0];  // Default initialization
    for (int t = 1; t < n; ++t) {
      if (y[t - 1] <= Q[t - 1]) {
        x[t] = gammas[0] + gammas[1] * (Q[t - 1] - y[t - 1]) + gammas[2] * x[t - 1];
      } else {
        x[t] = x[t - 1];
      }
    }
  } else {
    // Update the last element of x
    int t = n - 1;
    if (y[t - 1] <= Q[t - 1]) {
      x[t] = gammas[0] + gammas[1] * (Q[t - 1] - y[t - 1]) + gammas[2] * x[t - 1];
    } else {
      x[t] = x[t - 1];
    }
  }

  // Compute ES
  NumericVector es(n);
  for (int i = 0; i < n; ++i) {
    es[i] = Q[i] - x[i];
  }

  return es;
}
