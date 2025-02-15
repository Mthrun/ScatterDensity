#include <RcppArmadillo.h>
#include <RcppParallel.h>
#include <vector>
#include <cmath>

using namespace Rcpp;
using namespace RcppParallel;

//We define a PSphereWorker struct inheriting from RcppParallel::Worker. Its job is to loop over a range of bins (using a flattened index) and, for each bin, identify:
  
//The data points that lie exactly in the bin (center tile).
//The data points in the surrounding bins (neighbors).

// A custom worker to process bins in parallel
struct PSphereWorker : public Worker {
  // Inputs
  const Rcpp::NumericMatrix data;
  const Rcpp::IntegerVector xBinNr;
  const Rcpp::IntegerVector yBinNr;
  const unsigned int nrXBins;
  const unsigned int nrYBins;
  const unsigned int nrData;
  const double paretoRadius;
  
  // Output vector (shared among threads; each index is written only once)
  std::vector<int>& output;
  
  // Constructor: initialize all members by reference or value
  PSphereWorker(const Rcpp::NumericMatrix data,
                const Rcpp::IntegerVector xBinNr,
                const Rcpp::IntegerVector yBinNr,
                unsigned int nrXBins,
                unsigned int nrYBins,
                unsigned int nrData,
                double paretoRadius,
                std::vector<int>& output)
    : data(data), xBinNr(xBinNr), yBinNr(yBinNr),
      nrXBins(nrXBins), nrYBins(nrYBins), nrData(nrData),
      paretoRadius(paretoRadius), output(output) { }
  
  // This operator() will be called in parallel by RcppParallel::parallelFor
  void operator()(std::size_t begin, std::size_t end) {
    // Loop over bins using a flattened index.
    for (std::size_t index = begin; index < end; index++) {
      // Compute 2D bin coordinates from the flattened index.
      int i = index / nrYBins;
      int j = index % nrYBins;
      
      // For the current bin (i,j), collect indices of data points that
      // fall exactly in the center tile and those in its surrounding.
      std::vector<int> pointsInCenterTileInd;
      std::vector<int> pointsInSurroundingInd;
      
      for (unsigned int k = 0; k < nrData; k++) {
        if (xBinNr[k] == i && yBinNr[k] == j) {
          pointsInCenterTileInd.push_back(k);
        } else if (std::abs(xBinNr[k] - i) < 2 && std::abs(yBinNr[k] - j) < 2) {
          pointsInSurroundingInd.push_back(k);
        }
      }
      
      int nrInCenterTile = pointsInCenterTileInd.size();
      if (nrInCenterTile > 0) {
        // Combine the indices of center and surrounding points.
        std::vector<int> points;
        points.reserve(pointsInCenterTileInd.size() + pointsInSurroundingInd.size());
        points.insert(points.end(), pointsInCenterTileInd.begin(), pointsInCenterTileInd.end());
        points.insert(points.end(), pointsInSurroundingInd.begin(), pointsInSurroundingInd.end());
        
        unsigned int numrows = points.size();
        unsigned int numcols = data.ncol();
        
        // Build an Armadillo matrix from the selected rows.
        arma::mat matrixfordist(numrows, numcols);
        for (unsigned int row = 0; row < numrows; row++) {
          int datarow = points[row];
          for (unsigned int col = 0; col < numcols; col++) {
            matrixfordist(row, col) = data(datarow, col);
          }
        }
        
       // For each bin with at least one center point, we merge the indices (center and surrounding), extract the corresponding rows from the data matrix, and compute the pairwise Euclidean distances using Armadillo. For each center point, we then count how many of these distances are within paretoRadius.
        
        // Compute squared Euclidean distances via vectorized operations.
        arma::colvec rowSums = arma::sum(arma::square(matrixfordist), 1);
        arma::mat C = -2 * (matrixfordist * matrixfordist.t());
        C.each_col() += rowSums;
        C.each_row() += rowSums.t();
        arma::mat dists = arma::sqrt(C);
        
        // For each center point (first nrInCenterTile rows), count how many
        // points (including itself and the surrounding points) are within
        // the specified paretoRadius.
        for (int r = 0; r < nrInCenterTile; r++) {
          int count = 0;
          for (unsigned int c = 0; c < dists.n_cols; c++) {
            if (dists(r, c) <= paretoRadius)
              count++;
          }
          // Write the result to the output vector at the index corresponding
          // to the center tile point.
          output[pointsInCenterTileInd[r]] = count;
        }
      }
    }
  }
};

// [[Rcpp::depends(RcppArmadillo, RcppParallel)]]
// [[Rcpp::export]]
IntegerVector c_inPSphere2D_parallel(NumericMatrix data,
                                     IntegerVector xBinNr,
                                     IntegerVector yBinNr,
                                     unsigned int nrXBins,
                                     unsigned int nrYBins,
                                     unsigned int nrData,
                                     double paretoRadius) {
  // Initialize output vector with zeros.
  std::vector<int> output(nrData, 0);
  
  // Total number of bins (flattened loop over 2D bins).
  std::size_t totalBins = nrXBins * nrYBins;
  
  // Instantiate the worker.
  PSphereWorker worker(data, xBinNr, yBinNr, nrXBins, nrYBins, nrData, paretoRadius, output);
  
  //The parallelFor function from RcppParallel splits the bins among available threads. Since each data point appears only in one center tile, concurrent writes to the output vector are safe.
  // Parallel loop over bins.
  parallelFor(0, totalBins, worker);
  
  // Return the results as an R integer vector.
  return wrap(output);
}
