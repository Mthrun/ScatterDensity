#include <Rcpp.h>
using namespace Rcpp;

//this function replaces pracma::accumarray for sum, it is faster for larger dataset

// [[Rcpp::export]]
SEXP accumarray_rcpp(SEXP subs, NumericVector val, Nullable<NumericVector> sz = R_NilValue, double fillval = 0) {
  // Convert subs to a NumericMatrix.
  // If subs is a vector, make it an n x 1 matrix.
  NumericMatrix subsMat;
  bool wasMatrix = true;
  if (Rf_isMatrix(subs)) {
    subsMat = as<NumericMatrix>(subs);
  } else {
    NumericVector subsVec = as<NumericVector>(subs);
    int n = subsVec.size();
    subsMat = NumericMatrix(n, 1);
    for (int i = 0; i < n; i++){
      subsMat(i, 0) = subsVec[i];
    }
    wasMatrix = false;
  }
  
  // Get number of rows (n = number of values) and columns (m = number of dimensions)
  int n = subsMat.nrow();
  int m = subsMat.ncol();
  
  // Floor the subscripts and check that all indices are >= 1.
  for (int i = 0; i < n; i++){
    for (int j = 0; j < m; j++){
      subsMat(i,j) = std::floor(subsMat(i,j));
      if (subsMat(i,j) < 1)
        stop("All indices in 'subs' must be >= 1.");
    }
  }
  
  if (val.size() < n)
    stop("Length of 'val' must not be smaller than the number of rows of 'subs'.");
  
  // Determine the dimensions (dm) for the output.
  // For each column, take the maximum subscript.
  IntegerVector dm(m);
  for (int j = 0; j < m; j++){
    double maxVal = 0;
    for (int i = 0; i < n; i++){
      double cur = subsMat(i, j);
      if (cur > maxVal)
        maxVal = cur;
    }
    dm[j] = (int) maxVal;
  }
  
  // If a size vector (sz) was supplied, check its consistency and override dm.
  if (sz.isNotNull()) {
    NumericVector szVec(sz);
    if (szVec.size() != m)
      stop("Argument 'sz' must have length equal to the number of columns of 'subs'.");
    for (int j = 0; j < m; j++){
      if (szVec[j] < dm[j])
        stop("Argument 'sz' does not fit with 'subs'.");
      dm[j] = (int) szVec[j];
    }
  }
  
  // We'll implement the common case of summing the 'val' entries that map to each bin.
  // For bins that receive no value, we leave the entry as fillval.
  
  if (m == 1) {
    // 1-dimensional case.
    int outLength = dm[0];  // dm[0] is the maximum index along the only dimension
    NumericVector A(outLength, 0.0);
    std::vector<bool> updated(outLength, false);
    
    // Loop over all values and accumulate
    for (int i = 0; i < n; i++){
      int idx = (int)subsMat(i, 0) - 1; // convert to 0-indexed
      A[idx] += val[i];
      updated[idx] = true;
    }
    // For bins that never received a value, set to fillval.
    for (int i = 0; i < outLength; i++){
      if (!updated[i])
        A[i] = fillval;
    }
    // If the original subs was a matrix, return as a column matrix.
    if (wasMatrix)
      A.attr("dim") = Dimension(outLength, 1);
    return A;
    
  } else {
    // Multi-dimensional case.
    // Compute the total number of cells and the multipliers for linear indexing.
    int tot = 1;
    for (int j = 0; j < m; j++){
      tot *= dm[j];
    }
    NumericVector A(tot, 0.0);
    std::vector<bool> updated(tot, false);
    
    // In R arrays (and matrices) indexing is column-major.
    // Compute multipliers such that for a subscript vector (i1, i2, ..., im),
    // the linear (0-indexed) index is:
    //    (i1 - 1) + (i2 - 1)*dm[0] + (i3 - 1)*dm[0]*dm[1] + ... 
    std::vector<int> mult(m);
    mult[0] = 1;
    for (int j = 1; j < m; j++){
      mult[j] = mult[j - 1] * dm[j - 1];
    }
    
    // Loop over each row (each set of subscripts) and accumulate.
    for (int i = 0; i < n; i++){
      int lin = 0;
      for (int j = 0; j < m; j++){
        int sub = (int)subsMat(i, j);
        lin += (sub - 1) * mult[j];
      }
      A[lin] += val[i];
      updated[lin] = true;
    }
    
    // Fill bins that received no value with fillval.
    for (int i = 0; i < tot; i++){
      if (!updated[i])
        A[i] = fillval;
    }
    
    // Set the dimension attribute so that A is an array of dimension dm.
    A.attr("dim") = dm;
    return A;
  }
}
