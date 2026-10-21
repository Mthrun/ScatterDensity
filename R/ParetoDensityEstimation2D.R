# Minimal 2-D Pareto density extension --------------------------------------
#
# Essential steps only:
#   1. published 1-D PDE marginals;
#   2. paired rectangular product-PDE joint seed;
#   3. IPF to make the joint agree with the two 1-D PDE marginals.

.pde2d_weights <- function(grid) {
  
  ###############################################################################
  # .pde2d_weights
  #
  # INTERNAL HELPER
  #
  # Calculates trapezoidal quadrature weights for a strictly increasing,
  # finite one-dimensional grid.
  #
  # INPUT
  # grid              Numeric vector with at least two strictly increasing,
  #                   finite grid locations.
  #
  # OUTPUT
  # Numeric vector of the same length as grid containing trapezoidal integration
  # weights.
  #
  # ERROR HANDLING
  # Stops if grid has fewer than two values, contains non-finite values, or is not
  # strictly increasing.
  #
  # INTERNAL
  # Used by ParetoDensityEstimation2D() and Mutualinfo_pde_2d().
  ###############################################################################
  
  n <- length(grid)
  if (n < 2L || any(!is.finite(grid)) || any(diff(grid) <= 0)) {
    stop("PDE2D: invalid kernel grid.")
  }
  if (n == 2L) {
    return(rep(diff(grid) / 2, 2L))
  }
  c(
    grid[2L] - grid[1L],
    grid[3L:n] - grid[1L:(n - 2L)],
    grid[n] - grid[n - 1L]
  ) / 2
}

# Returns only the fields needed by ScatterDensity and mutual information:
ParetoDensityEstimation2D <- function(X, Y) {
  ###############################################################################
  # ParetoDensityEstimation2D
  #
  # Computes a minimal two-dimensional extension of Pareto Density Estimation
  # for two paired continuous variables.
  #
  # The estimator combines:
  #   1) published one-dimensional Pareto Density Estimates for X and Y;
  #   2) a paired rectangular product-PDE joint seed;
  #   3) iterative proportional fitting (IPF) to match the one-dimensional
  #      marginal probability masses.
  #
  # INPUT
  # X                 Numeric vector containing observations of the first
  #                   continuous variable.
  #
  # Y                 Numeric vector containing observations of the second
  #                   continuous variable. Must have the same length as X.
  #
  # OUTPUT
  # List with:
  #
  # Densities         Numeric vector with one bilinearly interpolated joint
  #                   density value for every observed pair (X[i], Y[i]).
  #
  # GridDensity       Numeric matrix containing the estimated joint density on
  #                   the rectangular grid. Rows correspond to Xkernels and
  #                   columns correspond to Ykernels.
  #
  # Xkernels          Numeric vector containing the grid locations for X.
  #
  # Ykernels          Numeric vector containing the grid locations for Y.
  #
  # DETAILS
  # The function first estimates one-dimensional Pareto densities with:
  #
  #   DataVisualizations::ParetoDensityEstimation()
  #
  # The Pareto radii and kernel grids returned by that function are reused.
  # Published one-dimensional density grids with more than 256 locations are
  # linearly interpolated to 256 grid points to limit quadratic memory growth.
  # This interpolation is not additional density smoothing.
  #
  # The paired rectangular joint seed is:
  #
  #   q(s,t) = sum_i I(|X_i-s| <= r_x) I(|Y_i-t| <= r_y)
  #
  # and therefore preserves the observed pairing between X and Y.
  #
  # The joint seed is converted to probability mass using trapezoidal quadrature.
  # A small positive support floor is added before applying IPF. IPF rescales rows
  # and columns until the joint mass matches the two one-dimensional Pareto
  # Density Estimate marginals, within tolerance 1e-8 or at most 500 iterations.
  #
  # The final grid density is obtained by dividing the fitted probability mass by
  # the two-dimensional quadrature weights. Density values at observed pairs are
  # calculated by bilinear interpolation.
  #
  # REQUIREMENTS
  # - X and Y must be numeric vectors of equal length.
  # - All observations must be finite.
  # - At least 50 paired observations are required.
  # - Each variable must contain at least 12 unique values.
  # - Suggested package DataVisualizations must be installed.
  #
  # ERROR HANDLING
  # Stops if:
  # - X or Y is not a numeric vector;
  # - the vectors have unequal lengths;
  # - non-finite observations are present;
  # - the minimum sample or uniqueness requirements are not met;
  # - DataVisualizations is unavailable;
  # - a Pareto radius, grid, density, or probability mass is invalid.
  #
  # A warning is issued if IPF does not reach a marginal error below 1e-8 within
  # 500 iterations.
  #
  # NOTE
  # This is an explicit paired rectangular two-dimensional extension constrained
  # by published one-dimensional Pareto Density Estimates. It is not claimed to
  # be a published two-dimensional hypersphere Pareto Density Estimator.
  #
  # DEPENDENCIES
  # Suggests: DataVisualizations
  #
  
  if (!is.numeric(X) || !is.numeric(Y) ||
      !is.null(dim(X)) || !is.null(dim(Y)) || length(X) != length(Y)) {
    stop("PDE2D: X and Y must be numeric vectors of equal length.")
  }
  if (any(!is.finite(X)) || any(!is.finite(Y))) {
    stop("PDE2D: remove non-finite paired observations first.")
  }
  if (length(X) < 50L ||
      length(unique(X)) < 12L || length(unique(Y)) < 12L) {
    stop("PDE2D: at least 50 pairs and 12 unique values per variable are required.")
  }
  if (!requireNamespace("DataVisualizations", quietly = TRUE)) {
    stop("PDE2D: package 'DataVisualizations' is required.")
  }
  
  # ParetoDensityEstimation() selects the kernel grid and Pareto radius and
  # returns both. Calling ParetoRadius_fast() separately would duplicate work.
  px_fit <- DataVisualizations::ParetoDensityEstimation(
    Data = X, PlotIt = FALSE, Compute = "Cpp", Silent = TRUE
  )
  py_fit <- DataVisualizations::ParetoDensityEstimation(
    Data = Y, PlotIt = FALSE, Compute = "Cpp", Silent = TRUE
  )
  
  rx <- as.numeric(px_fit$paretoRadius)[1L]
  ry <- as.numeric(py_fit$paretoRadius)[1L]
  if (!is.finite(rx) || rx <= 0 || !is.finite(ry) || ry <= 0) {
    stop("PDE2D: ParetoDensityEstimation() returned an invalid radius.")
  }
  
  # A 2-D grid grows quadratically. Keep the published 1-D estimate but
  # represent it on at most 256 points per axis. This is interpolation only,
  # not additional density smoothing.
  compact <- function(fit) {
    grid <- as.numeric(fit$kernels)
    density <- pmax(0, as.numeric(fit$paretoDensity))
    if (length(grid) != length(density) || any(!is.finite(density))) {
      stop("PDE2D: invalid output from ParetoDensityEstimation().")
    }
    if (length(grid) > 256L) {
      new_grid <- seq(grid[1L], grid[length(grid)], length.out = 256L)
      density <- stats::approx(grid, density, xout = new_grid, rule = 2)$y
      grid <- new_grid
    }
    list(grid = grid, density = density)
  }
  
  px_fit <- compact(px_fit)
  py_fit <- compact(py_fit)
  gx <- px_fit$grid
  gy <- py_fit$grid
  fx <- px_fit$density
  fy <- py_fit$density
  
  wx <- .pde2d_weights(gx)
  wy <- .pde2d_weights(gy)
  px <- fx * wx
  py <- fy * wy
  if (sum(px) <= 0 || sum(py) <= 0) {
    stop("PDE2D: a marginal PDE has zero mass.")
  }
  px <- px / sum(px)
  py <- py / sum(py)
  
  # q(s,t) = sum_i I(|X_i-s| <= rx) I(|Y_i-t| <= ry)
  counts <- matrix(0, length(gx), length(gy))
  chunk <- max(1L, floor(2e6 / max(length(gx), length(gy))))
  for (first in seq.int(1L, length(X), by = chunk)) {
    i <- first:min(length(X), first + chunk - 1L)
    in_x <- abs(outer(gx, X[i], `-`)) <= rx
    in_y <- abs(outer(gy, Y[i], `-`)) <= ry
    counts <- counts + tcrossprod(in_x, in_y)
  }
  
  mass <- counts * outer(wx, wy)
  positive <- mass[mass > 0]
  if (!length(positive)) {
    stop("PDE2D: the paired PDE has zero mass.")
  }
  mass <- mass + mean(positive) * 1e-12
  mass <- mass / sum(mass)
  
  # IPF: match the published 1-D PDE marginals.
  for (iteration in seq_len(500L)) {
    row_scale <- numeric(length(px))
    row_scale[px > 0] <- px[px > 0] / rowSums(mass)[px > 0]
    mass <- sweep(mass, 1L, row_scale, `*`)
    
    col_scale <- numeric(length(py))
    col_scale[py > 0] <- py[py > 0] / colSums(mass)[py > 0]
    mass <- sweep(mass, 2L, col_scale, `*`)
    
    error <- max(
      abs(rowSums(mass) - px),
      abs(colSums(mass) - py)
    )
    if (error < 1e-8) break
  }
  mass <- mass / sum(mass)
  if (error >= 1e-8) {
    warning("PDE2D: IPF marginal error = ", signif(error, 3L))
  }
  
  density <- mass / outer(wx, wy)
  
  # Bilinear density at each observed pair, for ScatterDensity coloring.
  ix <- pmax(1L, pmin(length(gx) - 1L, findInterval(X, gx)))
  iy <- pmax(1L, pmin(length(gy) - 1L, findInterval(Y, gy)))
  tx <- pmax(0, pmin(1, (X - gx[ix]) / (gx[ix + 1L] - gx[ix])))
  ty <- pmax(0, pmin(1, (Y - gy[iy]) / (gy[iy + 1L] - gy[iy])))
  
  point_density <-
    (1 - tx) * (1 - ty) * density[cbind(ix, iy)] +
    tx * (1 - ty) * density[cbind(ix + 1L, iy)] +
    (1 - tx) * ty * density[cbind(ix, iy + 1L)] +
    tx * ty * density[cbind(ix + 1L, iy + 1L)]
  
  list(
    Densities = point_density,
    GridDensity = density,
    Xkernels = gx,
    Ykernels = gy
  )
}