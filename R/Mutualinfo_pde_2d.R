Mutualinfo_pde_2d = function(x, y, base = 2, Jackknife = FALSE) {
  ###############################################################################
  # Mutualinfo_pde_2d
  #
  # Estimates mutual information between two continuous variables from the
  # two-dimensional density returned by ParetoDensityEstimation2D().
  #
  # INPUT
  # x                 Numeric vector containing observations of the first
  #                   continuous variable.
  #
  # y                 Numeric vector containing observations of the second
  #                   continuous variable.
  #
  # OPTIONAL
  # base              Numeric scalar defining the logarithm base.
  #                   Must be positive and different from 1.
  #                   Default: 2, yielding mutual information in bits.
  #
  # OUTPUT
  # Numeric scalar containing the estimated mutual information:
  #
  #   sum_{i,j} p(i,j) log_base(p(i,j)/(p_X(i) p_Y(j)))
  #
  # DETAILS
  # The joint density is converted to a normalized probability-mass matrix using
  # trapezoidal quadrature. Marginal masses are calculated from row and column
  # sums, and mutual information is evaluated on cells with positive joint mass.
  #
  # Small negative values caused only by floating-point error are truncated to 0.
  #
  # DEPENDENCIES
  # - ParetoDensityEstimation2D()
  # - .pde2d_weights()
  # - Suggested package DataVisualizations
  #
  # author: not specified
  ###############################################################################
  
  if (!is.logical(Jackknife) || length(Jackknife) != 1L || is.na(Jackknife)) {
    stop("Mutualinfo_pde_2d: Jackknife must be one non-missing logical value.")
  }
  if (isTRUE(Jackknife)) {
    return(Mutualinfo_pde_2d_jackknife(x, y, base = base))
  }
  
  if (!is.numeric(base) || length(base) != 1L ||
      !is.finite(base) || base <= 0 || base == 1) {
    stop("Mutualinfo_pde_2d: invalid logarithm base.")
  }
  
  fit <- ParetoDensityEstimation2D(x, y)
  wx <- .pde2d_weights(fit$Xkernels)
  wy <- .pde2d_weights(fit$Ykernels)
  joint <- fit$GridDensity * outer(wx, wy)
  joint <- joint / sum(joint)
  independent <- tcrossprod(rowSums(joint), colSums(joint))
  keep <- joint > 0
  
  value <- sum(joint[keep] * log(joint[keep] / independent[keep])) / log(base)
  if (value < 0 && value > -1e-12) value <- 0
  value
}
# Delete-group Jackknife ------------------------------------------------------
#
# The two-dimensional estimator contains a nonlinear IPF step.  To keep this
# function small and maintainable, the Jackknife treats Mutualinfo_pde_2d() as a
# black box and re-estimates the complete 2-D density after deleting each group.

Mutualinfo_pde_2d_jackknife = function(x,
                                       y,
                                       base = 2,
                                       groups = 20L,
                                       seed = NULL) {
  if (!is.numeric(x) || !is.numeric(y) ||
      !is.null(dim(x)) || !is.null(dim(y)) || length(x) != length(y)) {
    stop("Mutualinfo_pde_2d_jackknife: x and y must be numeric vectors of equal length.")
  }
  if (any(!is.finite(x)) || any(!is.finite(y))) {
    stop("Mutualinfo_pde_2d_jackknife: remove non-finite paired observations first.")
  }
  if (!is.numeric(base) || length(base) != 1L ||
      !is.finite(base) || base <= 0 || base == 1) {
    stop("Mutualinfo_pde_2d_jackknife: invalid logarithm base.")
  }
  
  n = length(x)
  if (!is.numeric(groups) || is.complex(groups) || length(groups) != 1L ||
      is.na(groups) || !is.finite(groups) || groups != floor(groups) ||
      groups < 2L || groups >= n) {
    stop("Mutualinfo_pde_2d_jackknife: groups must be an integer between 2 and n - 1.")
  }
  groups = as.integer(groups)
  
  groupIndex = rep(seq_len(groups), length.out = n)
  groupSizes = tabulate(groupIndex, nbins = groups)
  if (n - max(groupSizes) < 50L) {
    stop(paste0(
      "Mutualinfo_pde_2d_jackknife: each reduced sample must retain at least ",
      "50 observations; reduce the number of deleted observations per group."
    ))
  }
  
  if (!is.null(seed)) {
    if (!is.numeric(seed) || is.complex(seed) || length(seed) != 1L ||
        is.na(seed) || !is.finite(seed) || seed != floor(seed)) {
      stop("Mutualinfo_pde_2d_jackknife: seed must be NULL or one integer.")
    }
    seed = as.integer(seed)
    set.seed(seed)
    groupIndex = sample(groupIndex, length(groupIndex), replace = FALSE)
  }
  
  mi = Mutualinfo_pde_2d(x, y, base = base)
  deleted = numeric(groups)
  
  for (g in seq_len(groups)) {
    omit = which(groupIndex == g)
    deleted[g] = tryCatch(
      Mutualinfo_pde_2d(x[-omit], y[-omit], base = base),
      error = function(e) {
        stop(
          paste0(
            "Mutualinfo_pde_2d_jackknife: group ", g,
            " could not be evaluated: ", conditionMessage(e)
          ),
          call. = FALSE
        )
      }
    )
  }
  
  # Generalized group pseudovalues.  For equal group sizes this reduces to
  # groups * mi - (groups - 1) * deleted[g].
  pseudoValues = (
    n * mi - (n - groupSizes) * deleted
  ) / groupSizes
  weights = groupSizes / n
  bcmi = sum(weights * pseudoValues)
  
  weightSquares = sum(weights^2)
  denominator = 1 - weightSquares
  se = if (denominator > 0) {
    sqrt(
      sum(weights * (pseudoValues - bcmi)^2) *
        weightSquares / denominator
    )
  } else {
    NA_real_
  }
  zvalue = if (is.finite(se) && se > 0) bcmi / se else NA_real_
  
  list(
    mi = mi,
    bcmi = bcmi,
    bias = mi - bcmi,
    se = se,
    zvalue = zvalue,
    n = n,
    base = base,
    groups = groups,
    groupSizes = groupSizes,
    method = "delete-group Jackknife with complete 2-D PDE and IPF refitting"
  )
}