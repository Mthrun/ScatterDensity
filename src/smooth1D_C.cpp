#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]

using namespace Rcpp;

// [[Rcpp::export]]
arma::mat smooth1D_C(arma::mat Y, double lambda, bool na_rm = false, bool Silent = false) {
  if (Y.n_cols == 1) {
    if (!Silent) {
      Rcpp::warning("smooth1D expected matrix. Converting to matrix.");
    }
    Y.reshape(Y.n_rows, 1);
  }
  
  if (na_rm) {
    Y.elem(find_nonfinite(Y)).zeros();
  }
  
  int m = Y.n_rows;
  arma::mat E = arma::eye(m, m);
  
  arma::mat D1 = arma::diff(E, 1, 0); // First-order differences
  arma::mat D2 = arma::diff(D1, 1, 0); // Second-order differences
  
  arma::mat P = std::pow(lambda, 2) * (D2.t() * D2) + 2 * lambda * (D1.t() * D1);
  
  arma::mat Z = arma::solve(E + P, Y);
  
  return Z;
}
