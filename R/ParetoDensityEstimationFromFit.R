ParetoDensityEstimationFromFit = function(V,Omit = NULL) {
  #  Vnew = ParetoDensityEstimationFromFit(V,Omit)
  #  Re-evaluates a Pareto density after omitting one or more observations.
  #  Radius, full-sample boundaries, counting grid, interpolation decision,
  #  and interpolation target remain fixed.
  #
  #  INPUT
  #  V				Result returned by ParetoDensityEstimation(). V$Fit must
  #				be non-NULL.
  #
  #  OPTIONAL
  #  Omit				Finite integer indices of observations in V$Fit$Data to omit.
  #				Repeated indices are ignored. NULL performs no omission.
  #
  #  OUTPUT
  #  List with
  #  kernels			Data values at which paretoDensity is returned.
  #  paretoDensity		Pareto-density values after omission.
  #  paretoRadius			Unchanged Pareto radius from V.
  #
  #  Revision: 2026
  
  requiredResultNames = c("kernels", "paretoRadius", "Fit")
  if (!is.list(V) || !all(requiredResultNames %in% names(V))) {
    stop("V is not a valid ParetoDensityEstimation result.", call. = FALSE)
  }
  
  Fit = V$Fit
  if (is.null(Fit)) {
    stop("A numerical Dirac representation has no Pareto-density fit context.", call. = FALSE)
  }
  
  requiredFitNames = c("Data", "paretoCounts", "countingKernels")
  if (!is.list(Fit) || !all(requiredFitNames %in% names(Fit))) {
    stop("V$Fit is not a valid Pareto-density fit context.", call. = FALSE)
  }
  
  Data = Fit$Data
  paretoRadius = V$paretoRadius
  paretoCounts = Fit$paretoCounts
  kernels = V$kernels
  countingKernels = Fit$countingKernels
  interpolationUsed = if ("interpolationUsed" %in% names(Fit)) {
    Fit$interpolationUsed
  } else {
    # Compatibility with fit objects created before interpolationUsed was
    # stored explicitly. This infers all cases except the rare situation in
    # which compatible_2020 interpolated onto an identical grid.
    !is.null(countingKernels)
  }
  if (!is.logical(interpolationUsed) || length(interpolationUsed) != 1L ||
      is.na(interpolationUsed)) {
    stop("V$Fit contains an invalid interpolation flag.", call. = FALSE)
  }
  minData = min(Data)
  maxData = max(Data)
  # Older fits used the empirical boundaries for reflection in every case.
  boundaryReflection = if ("boundaryReflection" %in% names(Fit)) Fit$boundaryReflection else TRUE
  if (!is.logical(boundaryReflection) || length(boundaryReflection) != 1L || is.na(boundaryReflection)) {
    stop("V$Fit contains an invalid boundary-reflection flag.", call. = FALSE)
  }
  
  if (is.null(countingKernels)) {
    evaluationKernels = kernels
  } else {
    evaluationKernels = countingKernels
  }
  
  if (length(paretoCounts) != length(evaluationKernels) ||
      any(!is.finite(paretoCounts)) || any(paretoCounts < 0)) {
    stop("V$Fit contains invalid Pareto-density counts.", call. = FALSE)
  }
  
  if (is.null(Omit) || length(Omit) == 0L) {
    Omit = integer(0L)
  } else {
    if (!is.numeric(Omit) || is.complex(Omit) || anyNA(Omit) ||
        any(!is.finite(Omit)) || any(Omit != floor(Omit))) {
      stop("Omit must contain finite integer observation indices.", call. = FALSE)
    }
    Omit = unique(as.integer(Omit))
    if (any(Omit < 1L) || any(Omit > length(Data))) {
      stop("Omit contains observation indices outside V$Fit$Data.", call. = FALSE)
    }
  }
  
  lowerBounds = evaluationKernels - paretoRadius
  upperBounds = evaluationKernels + paretoRadius
  
  for (i in Omit) {
    observation = Data[i]
    isInParetoSphere = observation >= lowerBounds & observation <= upperBounds
    paretoCounts[isInParetoSphere] = paretoCounts[isInParetoSphere] - 1
    
    if (boundaryReflection && observation < minData + paretoRadius) {
      lowerReflection = 2 * minData - observation
      isInParetoSphere = lowerReflection >= lowerBounds & lowerReflection <= upperBounds
      paretoCounts[isInParetoSphere] = paretoCounts[isInParetoSphere] - 1
    }
    if (boundaryReflection && observation > maxData - paretoRadius) {
      upperReflection = 2 * maxData - observation
      isInParetoSphere = upperReflection >= lowerBounds & upperReflection <= upperBounds
      paretoCounts[isInParetoSphere] = paretoCounts[isInParetoSphere] - 1
    }
  }
  
  if (any(!is.finite(paretoCounts)) || any(paretoCounts < 0)) {
    stop("Omitting observations produced invalid Pareto-density counts.", call. = FALSE)
  }
  
  dx = diff(evaluationKernels)
  area = sum(dx * (head(paretoCounts, -1L) + tail(paretoCounts, -1L)) / 2)
  if (!is.finite(area) || area <= 0) {
    stop("The Pareto density could not be normalized on the selected grid.", call. = FALSE)
  }
  densityOnEvaluationGrid = paretoCounts / area
  
  if (isTRUE(interpolationUsed)) {
    paretoDensity = stats::approx(
      x = evaluationKernels,
      y = densityOnEvaluationGrid,
      xout = kernels,
      rule = 1,
      ties = "ordered"
    )$y
    paretoDensity[!is.finite(paretoDensity)] = 0
  } else {
    paretoDensity = densityOnEvaluationGrid
  }
  
  return(list(
    kernels = kernels,
    paretoDensity = paretoDensity,
    paretoRadius = paretoRadius
  ))
}