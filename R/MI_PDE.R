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

# Unified mutual information estimator ---------------------------------------
#
# Continuous-continuous:
#   delegates to Mutualinfo_pde_2d(), which uses ParetoDensityEstimation2D().
#
# Continuous-discrete:
#   uses class-conditional one-dimensional Pareto density estimates on one
#   common grid and with one class-independent Pareto radius.  The densities
#   are converted to a joint probability-mass matrix before MI is calculated.
#
# Discrete-discrete:
#   uses the empirical contingency table.
#
# The same probability-mass MI calculation and logarithm base are used in all
# branches.  The old compiled c_mutualinfo() routine is therefore unnecessary.

.mi_from_mass = function(joint, base = 2,
                          tol = .Machine$double.eps * 1000) {
  joint = as.matrix(joint)
  storage.mode(joint) = "double"

  if (!length(joint) || any(!is.finite(joint))) {
    stop("mutualinfo: joint mass must be a finite, non-empty matrix.")
  }

  scale = max(1, max(abs(joint)))
  if (any(joint < -tol * scale)) {
    stop("mutualinfo: substantial negative probability mass was produced.")
  }
  joint[joint < 0] = 0

  total = sum(joint)
  if (!is.finite(total) || total <= 0) {
    stop("mutualinfo: joint probability mass is empty.")
  }

  joint = joint / total
  independent = tcrossprod(rowSums(joint), colSums(joint))
  keep = joint > 0 & independent > 0

  value = sum(
    joint[keep] * (log(joint[keep]) - log(independent[keep]))
  ) / log(base)

  if (value < 0 && value > -sqrt(.Machine$double.eps)) value = 0
  value
}

.mi_grid_weights = function(grid) {
  grid = as.numeric(grid)
  n = length(grid)

  if (n < 2L || any(!is.finite(grid)) || any(diff(grid) <= 0)) {
    stop("mutualinfo: invalid PDE kernel grid.")
  }

  if (n == 2L) return(rep(diff(grid) / 2, 2L))

  c(
    grid[2L] - grid[1L],
    grid[3L:n] - grid[1L:(n - 2L)],
    grid[n] - grid[n - 1L]
  ) / 2
}

.mi_histogram_joint = function(continuous, discrete) {
  discrete = droplevels(factor(discrete))

  if (length(unique(continuous)) < 2L) {
    return(matrix(table(discrete), nrow = 1L))
  }

  h = graphics::hist(
    continuous,
    breaks = "Scott",
    plot = FALSE,
    right = TRUE,
    include.lowest = TRUE
  )

  bin = findInterval(
    continuous,
    h$breaks,
    rightmost.closed = TRUE,
    all.inside = TRUE
  )

  as.matrix(table(
    factor(bin, levels = seq_len(length(h$breaks) - 1L)),
    discrete
  ))
}

.mi_histogram_mass_on_grid = function(z, grid, weights) {
  m = length(grid)

  if (length(unique(z)) < 2L) {
    mass = numeric(m)
    mass[which.min(abs(grid - z[1L]))] = 1
    return(mass)
  }

  h = graphics::hist(
    z,
    breaks = "Scott",
    plot = FALSE,
    probability = TRUE,
    right = TRUE,
    include.lowest = TRUE
  )

  bin = findInterval(
    grid,
    h$breaks,
    rightmost.closed = TRUE,
    all.inside = FALSE
  )

  density = numeric(m)
  valid = bin >= 1L & bin <= length(h$density)
  density[valid] = h$density[bin[valid]]
  mass = density * weights

  # A very narrow class can fall between all requested kernel locations.
  # In that exceptional case, use the nearest common-grid cells directly.
  if (!is.finite(sum(mass)) || sum(mass) <= 0) {
    middle = (grid[-1L] + grid[-m]) / 2
    index = findInterval(
      z,
      c(-Inf, middle, Inf),
      rightmost.closed = TRUE,
      all.inside = TRUE
    )
    mass = tabulate(index, nbins = m)
  }

  mass / sum(mass)
}

.mi_class_mass = function(z, radius, grid, weights,
                           tol = .Machine$double.eps * 1000) {
  # The published practical rule uses PDE only with more than 50 observations
  # and at least 12 distinct values.  Otherwise use the histogram fallback.
  if (length(z) <= 50L || length(unique(z)) < 12L) {
    return(.mi_histogram_mass_on_grid(z, grid, weights))
  }

  fit = tryCatch(
    DataVisualizations::ParetoDensityEstimation(
      Data = z,
      paretoRadius = radius,
      kernels = grid,
      PlotIt = FALSE,
      Compute = "Cpp",
      Silent = TRUE
    ),
    error = function(e) NULL
  )

  if (is.null(fit)) {
    return(.mi_histogram_mass_on_grid(z, grid, weights))
  }

  density = as.numeric(fit$paretoDensity)
  if (length(density) != length(grid) || any(!is.finite(density))) {
    return(.mi_histogram_mass_on_grid(z, grid, weights))
  }

  scale = max(1, max(abs(density)))
  if (any(density < -tol * scale)) {
    stop("mutualinfo: ParetoDensityEstimation() returned negative density.")
  }
  density[density < 0] = 0

  mass = density * weights
  if (!is.finite(sum(mass)) || sum(mass) <= 0) {
    return(.mi_histogram_mass_on_grid(z, grid, weights))
  }

  mass / sum(mass)
}

.mi_continuous_discrete_joint = function(continuous, discrete,
                                          tol = .Machine$double.eps * 1000) {
  discrete = droplevels(factor(discrete))

  if (nlevels(discrete) < 2L) {
    return(matrix(1, nrow = 1L, ncol = 1L))
  }

  # Below the range in which PDE was empirically investigated, use a direct
  # Scott-histogram joint estimate for the whole mixed pair.
  if (length(continuous) <= 50L || length(unique(continuous)) < 12L) {
    return(.mi_histogram_joint(continuous, discrete))
  }

  if (!requireNamespace("DataVisualizations", quietly = TRUE)) {
    stop("mutualinfo: package 'DataVisualizations' is required.")
  }

  # Estimate the grid and the single class-independent Pareto radius from all
  # observations of the continuous feature.
  global = tryCatch(
    DataVisualizations::ParetoDensityEstimation(
      Data = continuous,
      PlotIt = FALSE,
      Compute = "Cpp",
      Silent = TRUE
    ),
    error = function(e) NULL
  )

  if (is.null(global)) {
    return(.mi_histogram_joint(continuous, discrete))
  }

  grid = as.numeric(global$kernels)
  radius = as.numeric(global$paretoRadius)[1L]

  if (length(grid) < 2L || any(!is.finite(grid)) || any(diff(grid) <= 0) ||
      !is.finite(radius) || radius <= 0) {
    return(.mi_histogram_joint(continuous, discrete))
  }

  weights = .mi_grid_weights(grid)
  class_index = as.integer(discrete)
  class_names = levels(discrete)
  class_counts = tabulate(class_index, nbins = length(class_names))
  priors = class_counts / length(discrete)

  joint = matrix(
    0,
    nrow = length(grid),
    ncol = length(class_names),
    dimnames = list(NULL, class_names)
  )

  for (j in seq_along(class_names)) {
    conditional_mass = .mi_class_mass(
      continuous[class_index == j],
      radius = radius,
      grid = grid,
      weights = weights,
      tol = tol
    )
    joint[, j] = priors[j] * conditional_mass
  }

  joint
}

MI_PDE = function(x,
                       y,
                       isXDiscrete = FALSE,
                       isYDiscrete = FALSE,
                       eps = .Machine$double.eps * 1000,
                       useMPMI = FALSE,
                       na.rm = FALSE,
                       base = 2,
                       Jackknife = FALSE) {
  if (!is.logical(Jackknife) || length(Jackknife) != 1L || is.na(Jackknife)) {
    stop("mutualinfo: Jackknife must be one non-missing logical value.")
  }
  if (isTRUE(Jackknife)) {
    if (isTRUE(useMPMI)) {
      stop("mutualinfo: use mpmi's bcmi output directly when useMPMI = TRUE.")
    }
    return(MI_PDE_jackknife(
      x, y,
      isXDiscrete = isXDiscrete,
      isYDiscrete = isYDiscrete,
      eps = eps,
      na.rm = na.rm,
      base = base
    ))
  }

  if (!is.atomic(x) || !is.atomic(y) ||
      !is.null(dim(x)) || !is.null(dim(y))) {
    stop("mutualinfo: x and y must be vectors.")
  }
  if (length(x) != length(y)) {
    stop("mutualinfo: x and y must have the same length.")
  }

  isXDiscrete = isTRUE(isXDiscrete)
  isYDiscrete = isTRUE(isYDiscrete)

  if (!is.numeric(base) || length(base) != 1L ||
      !is.finite(base) || base <= 0 || base == 1) {
    stop("mutualinfo: base must be positive and different from 1.")
  }
  if (!is.numeric(eps) || length(eps) != 1L ||
      !is.finite(eps) || eps < 0) {
    stop("mutualinfo: eps must be a finite non-negative scalar.")
  }

  if (!isXDiscrete && !is.numeric(x)) {
    stop("mutualinfo: x must be numeric when treated as continuous.")
  }
  if (!isYDiscrete && !is.numeric(y)) {
    stop("mutualinfo: y must be numeric when treated as continuous.")
  }

  observed = function(z, discrete) {
    if (!discrete) return(is.finite(z))
    if (is.numeric(z) && !is.factor(z)) return(is.finite(z))
    !is.na(z)
  }

  complete = observed(x, isXDiscrete) & observed(y, isYDiscrete)
  if (!all(complete)) {
    if (isTRUE(na.rm)) {
      x = x[complete]
      y = y[complete]
    } else {
      stop(paste0(
        "mutualinfo: missing or non-finite paired observations are present; ",
        "use na.rm = TRUE to remove them."
      ))
    }
  }

  if (length(x) < 3L ||
      length(unique(x)) < 2L || length(unique(y)) < 2L) {
    return(0)
  }

  if (isTRUE(useMPMI)) {
    if (!requireNamespace("mpmi", quietly = TRUE)) {
      stop("mutualinfo: package 'mpmi' is required for useMPMI = TRUE.")
    }

    raw = if (!isXDiscrete && !isYDiscrete) {
      mpmi::cmi.pw(x, y)$mi
    } else if (isXDiscrete && isYDiscrete) {
      mpmi::dmi.pw(x, y)$mi
    } else if (!isXDiscrete) {
      mpmi::mmi.pw(x, y)$mi
    } else {
      mpmi::mmi.pw(y, x)$mi
    }

    # mpmi evaluates the logarithm naturally; convert nats to the chosen base.
    return(as.numeric(raw)[1L] / log(base))
  }

  if (!isXDiscrete && !isYDiscrete) {
    if (!exists("Mutualinfo_pde_2d", mode = "function", inherits = TRUE)) {
      stop(paste0(
        "mutualinfo: source the file containing Mutualinfo_pde_2d() and ",
        "ParetoDensityEstimation2D() before using two continuous variables."
      ))
    }
    return(Mutualinfo_pde_2d(x, y, base = base))
  }

  if (isXDiscrete && isYDiscrete) {
    joint = table(x, y, useNA = "no")
    return(.mi_from_mass(joint, base = base, tol = eps))
  }

  if (!isXDiscrete) {
    joint = .mi_continuous_discrete_joint(x, y, tol = eps)
  } else {
    joint = .mi_continuous_discrete_joint(y, x, tol = eps)
  }

  .mi_from_mass(joint, base = base, tol = eps)
}

# Jackknife correction --------------------------------------------------------

.mi_jk_mass_from_pde = function(V, grid, weights,
                                tol = .Machine$double.eps * 1000) {
  density = as.numeric(V$paretoDensity)
  returnedGrid = as.numeric(V$kernels)

  if (length(density) != length(grid) ||
      length(returnedGrid) != length(grid) ||
      any(!is.finite(density)) || any(!is.finite(returnedGrid)) ||
      !isTRUE(all.equal(returnedGrid, grid, tolerance = 0))) {
    stop("MI_PDE_jackknife: invalid Pareto-density result on the common grid.")
  }

  scale = max(1, max(abs(density)))
  if (any(density < -tol * scale)) {
    stop("MI_PDE_jackknife: ParetoDensityEstimation returned negative density.")
  }
  density[density < 0] = 0

  mass = density * weights
  total = sum(mass)
  if (!is.finite(total) || total <= 0) {
    stop("MI_PDE_jackknife: the Pareto density has zero probability mass.")
  }
  mass / total
}

.mi_jk_summary = function(mi, deleted, base, method) {
  n = length(deleted)
  bcmi = n * mi - (n - 1L) * mean(deleted)
  pseudoValues = n * mi - (n - 1L) * deleted
  se = stats::sd(pseudoValues) / sqrt(n)

  list(
    mi = mi,
    bcmi = bcmi,
    bias = mi - bcmi,
    se = se,
    zvalue = if (is.finite(se) && se > 0) bcmi / se else NA_real_,
    n = n,
    base = base,
    method = method
  )
}

MI_PDE_jackknife = function(x,
                             y,
                             isXDiscrete = FALSE,
                             isYDiscrete = FALSE,
                             eps = .Machine$double.eps * 1000,
                             na.rm = FALSE,
                             base = 2,
                             groups2d = 20L,
                             seed = NULL) {
  if (!is.atomic(x) || !is.atomic(y) ||
      !is.null(dim(x)) || !is.null(dim(y)) || length(x) != length(y)) {
    stop("MI_PDE_jackknife: x and y must be vectors of equal length.")
  }
  
  isXDiscrete = isTRUE(isXDiscrete)
  isYDiscrete = isTRUE(isYDiscrete)
  
  if (!is.numeric(base) || length(base) != 1L ||
      !is.finite(base) || base <= 0 || base == 1) {
    stop("MI_PDE_jackknife: base must be positive and different from 1.")
  }
  if (!is.numeric(eps) || length(eps) != 1L ||
      !is.finite(eps) || eps < 0) {
    stop("MI_PDE_jackknife: eps must be a finite non-negative scalar.")
  }
  if (!isXDiscrete && !is.numeric(x)) {
    stop("MI_PDE_jackknife: x must be numeric when treated as continuous.")
  }
  if (!isYDiscrete && !is.numeric(y)) {
    stop("MI_PDE_jackknife: y must be numeric when treated as continuous.")
  }
  
  observed = function(z, discrete) {
    if (!discrete) return(is.finite(z))
    if (is.numeric(z) && !is.factor(z)) return(is.finite(z))
    !is.na(z)
  }
  complete = observed(x, isXDiscrete) & observed(y, isYDiscrete)
  if (!all(complete)) {
    if (isTRUE(na.rm)) {
      x = x[complete]
      y = y[complete]
    } else {
      stop(paste0(
        "MI_PDE_jackknife: missing or non-finite paired observations are ",
        "present; use na.rm = TRUE to remove them."
      ))
    }
  }
  
  n = length(x)
  if (n < 3L || length(unique(x)) < 2L || length(unique(y)) < 2L) {
    stop("MI_PDE_jackknife: at least three non-constant complete pairs are required.")
  }
  
  if (!isXDiscrete && !isYDiscrete) {
    if (!exists("Mutualinfo_pde_2d_jackknife", mode = "function", inherits = TRUE)) {
      stop("MI_PDE_jackknife: source Mutualinfo_pde_2d_jackknife() first.")
    }
    return(Mutualinfo_pde_2d_jackknife(
      x, y, base = base, groups = groups2d, seed = seed
    ))
  }
  
  if (isXDiscrete && isYDiscrete) {
    xf = factor(x)
    yf = factor(y)
    xi = as.integer(xf)
    yi = as.integer(yf)
    joint = table(xf, yf, useNA = "no")
    mi = .mi_from_mass(joint, base = base, tol = eps)
    deleted = numeric(n)
    
    for (i in seq_len(n)) {
      jointMinusI = joint
      jointMinusI[xi[i], yi[i]] = jointMinusI[xi[i], yi[i]] - 1L
      deleted[i] = .mi_from_mass(jointMinusI, base = base, tol = eps)
    }
    
    return(.mi_jk_summary(
      mi, deleted, base,
      "exact leave-one-out on the empirical contingency table"
    ))
  }
  
  if (!exists("ParetoDensityEstimationFromFit", mode = "function", inherits = TRUE)) {
    stop(paste0(
      "MI_PDE_jackknife: ParetoDensityEstimationFromFit() must be exported ",
      "or sourced before the fast mixed Jackknife is used."
    ))
  }
  if (!requireNamespace("DataVisualizations", quietly = TRUE)) {
    stop("MI_PDE_jackknife: package 'DataVisualizations' is required.")
  }
  
  if (!isXDiscrete) {
    continuous = x
    discrete = y
  } else {
    continuous = y
    discrete = x
  }
  discrete = droplevels(factor(discrete))
  classIndex = as.integer(discrete)
  classNames = levels(discrete)
  classCounts = tabulate(classIndex, nbins = length(classNames))
  
  if (length(classNames) < 2L || length(continuous) <= 50L ||
      length(unique(continuous)) < 12L) {
    stop(paste0(
      "MI_PDE_jackknife: the fast PDE Jackknife requires more than 50 ",
      "observations and at least 12 distinct continuous values."
    ))
  }
  
  uniquePerClass = vapply(
    seq_along(classNames),
    function(j) length(unique(continuous[classIndex == j])),
    integer(1L)
  )
  if (any(classCounts <= 50L) || any(uniquePerClass < 12L)) {
    stop(paste0(
      "MI_PDE_jackknife: every class must contain more than 50 observations ",
      "and at least 12 distinct continuous values."
    ))
  }
  
  global = DataVisualizations::ParetoDensityEstimation(
    Data = continuous,
    PlotIt = FALSE,
    Compute = "Cpp",
    Silent = TRUE
  )
  grid = as.numeric(global$kernels)
  radius = as.numeric(global$paretoRadius)[1L]
  weights = .mi_grid_weights(grid)
  
  if (length(grid) < 2L || any(!is.finite(grid)) || any(diff(grid) <= 0) ||
      !is.finite(radius) || radius <= 0) {
    stop("MI_PDE_jackknife: invalid global Pareto-density fit.")
  }
  
  classFits = vector("list", length(classNames))
  conditionalMass = matrix(0, length(grid), length(classNames))
  localIndex = integer(n)
  
  for (j in seq_along(classNames)) {
    indices = which(classIndex == j)
    localIndex[indices] = seq_along(indices)
    classFits[[j]] = DataVisualizations::ParetoDensityEstimation(
      Data = continuous[indices],
      paretoRadius = radius,
      kernels = grid,
      PlotIt = FALSE,
      Compute = "Cpp",
      Silent = TRUE
    )
    if (is.null(classFits[[j]]$Fit)) {
      stop("MI_PDE_jackknife: a class returned no Pareto-density fit context.")
    }
    conditionalMass[, j] = .mi_jk_mass_from_pde(
      classFits[[j]], grid, weights, tol = eps
    )
  }
  
  priors = classCounts / n
  joint = sweep(conditionalMass, 2L, priors, `*`)
  mi = .mi_from_mass(joint, base = base, tol = eps)
  deleted = numeric(n)
  
  for (i in seq_len(n)) {
    j = classIndex[i]
    VminusI = ParetoDensityEstimationFromFit(
      classFits[[j]],
      Omit = localIndex[i]
    )
    massMinusI = .mi_jk_mass_from_pde(
      VminusI, grid, weights, tol = eps
    )
    
    countsMinusI = classCounts
    countsMinusI[j] = countsMinusI[j] - 1L
    priorsMinusI = countsMinusI / (n - 1L)
    jointMinusI = sweep(conditionalMass, 2L, priorsMinusI, `*`)
    jointMinusI[, j] = priorsMinusI[j] * massMinusI
    deleted[i] = .mi_from_mass(jointMinusI, base = base, tol = eps)
  }
  
  result = .mi_jk_summary(
    mi,
    deleted,
    base,
    "fast fixed-context leave-one-out using ParetoDensityEstimationFromFit()"
  )
  result$paretoRadius = radius
  result
}
