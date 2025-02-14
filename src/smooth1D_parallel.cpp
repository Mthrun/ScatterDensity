#include <RcppArmadillo.h>
#include <RcppParallel.h>
// [[Rcpp::depends(RcppArmadillo, RcppParallel)]]

using namespace Rcpp;
using namespace RcppParallel;

struct Smooth1DWorker : public Worker {
  const arma::mat& Y;
  arma::mat& Z;
  const arma::mat& E;
  const arma::mat& P;
  
  Smooth1DWorker(const arma::mat& Y, arma::mat& Z, const arma::mat& E, const arma::mat& P)
    : Y(Y), Z(Z), E(E), P(P) {}
  
  void operator()(std::size_t start, std::size_t end) {
    for (std::size_t i = start; i < end; i++) {
      Z.col(i) = arma::solve(E + P, Y.col(i));
    }
  }
};

// [[Rcpp::export]]
arma::mat smooth1D_parallel(arma::mat Y, double lambda, bool na_rm = false, bool Silent = false) {
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
  
  arma::mat Z(Y.n_rows, Y.n_cols);
  Smooth1DWorker worker(Y, Z, E, P);
  parallelFor(0, Y.n_cols, worker);
  
  return Z;
}
