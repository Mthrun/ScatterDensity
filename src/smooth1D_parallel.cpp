#include <RcppArmadillo.h>
#include <RcppParallel.h>
// [[Rcpp::depends(RcppArmadillo, RcppParallel)]]

using namespace Rcpp;
using namespace RcppParallel;

// Parallel worker for row-wise or column-wise smoothing
struct Smooth2DWorker : public Worker {
  const arma::mat& G;
  arma::mat& Z;
  const arma::mat& P_inv;
  bool smooth_rows;
  
  Smooth2DWorker(const arma::mat& G, arma::mat& Z, const arma::mat& P_inv, bool smooth_rows)
    : G(G), Z(Z), P_inv(P_inv), smooth_rows(smooth_rows) {}
  
  void operator()(std::size_t start, std::size_t end) {
    if (smooth_rows) {
      for (std::size_t i = start; i < end; i++) {
        Z.row(i) = G.row(i) * P_inv;  // Row-wise smoothing
      }
    } else {
      for (std::size_t i = start; i < end; i++) {
        Z.col(i) = P_inv * G.col(i);  // Column-wise smoothing
      }
    }
  }
};

// [[Rcpp::export]]
arma::mat smooth1D_parallel(arma::mat G, double lambda, int nbins, bool smooth_rows = true, bool na_rm = false, bool Silent = false) {
  if (na_rm) {
    G.elem(find_nonfinite(G)).zeros();
  }
  
  int m = smooth_rows ? G.n_cols : G.n_rows;
  arma::mat E = arma::eye(m, m);
  
  arma::mat D1 = arma::diff(E, 1, 0); // First-order differences
  arma::mat D2 = arma::diff(D1, 1, 0); // Second-order differences
  
  arma::mat P = std::pow(lambda, 2) * (D2.t() * D2) + 2 * lambda * (D1.t() * D1);
  arma::mat P_inv = arma::inv(E + P);  // Precompute inverse
  
  arma::mat Z(G.n_rows, G.n_cols);
  Smooth2DWorker worker(G, Z, P_inv, smooth_rows);
  parallelFor(0, smooth_rows ? G.n_rows : G.n_cols, worker);
  
  return Z;
}
