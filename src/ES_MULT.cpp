#include <Rcpp.h>
using namespace Rcpp;

// [[Rcpp::export]]
NumericVector mult_ES_C(NumericVector y, NumericVector Q, NumericVector gammas) {
  // Compute es as (1 + exp(gamma0)) * Q
  NumericVector es = (1 + std::exp(gammas[0])) * Q;

  return es;
}
