#include <Rcpp.h>
using namespace Rcpp;
using namespace std;

// [[Rcpp::export]]
Rcpp::NumericVector quantile4LargeVectors(Rcpp::NumericVector x, Rcpp::NumericVector probs) {
  
  NumericVector y = clone(x);
  std::sort(y.begin(), y.end());
  return y[x.size()*(probs - 0.000000001)];
}