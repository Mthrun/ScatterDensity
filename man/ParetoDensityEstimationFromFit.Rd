\name{ParetoDensityEstimationFromFit}
\alias{ParetoDensityEstimationFromFit}
\title{Re-evaluate a Pareto Density After Omitting Observations}
\description{
Re-evaluates a previously fitted one-dimensional Pareto density after omitting
one or more observations. The original Pareto radius, full-sample boundaries,
counting grid, interpolation decision, and interpolation target remain fixed.
}
\usage{
ParetoDensityEstimationFromFit(V, Omit = NULL)
}
\arguments{
  \item{V}{
    Result returned by \code{ParetoDensityEstimation()}. It must contain
    \code{kernels}, \code{paretoRadius}, and a non-\code{NULL} \code{Fit}
    component. The fit context must contain \code{Data},
    \code{paretoCounts}, and \code{countingKernels}.
  }
  \item{Omit}{
    Optional numeric vector of finite integer indices referring to observations
    in \code{V$Fit$Data}. Repeated indices are ignored. \code{NULL} or an empty
    vector performs no omission.
  }
}
\value{
A list with components:

\item{kernels}{
  Numeric vector of evaluation locations, unchanged from \code{V$kernels}.
}
\item{paretoDensity}{
  Numeric vector of Pareto-density values after removing the requested
  observations and renormalizing the remaining counts.
}
\item{paretoRadius}{
  Unchanged Pareto radius from \code{V}.
}
}
\details{
The function supports leave-one-out or leave-several-out re-evaluation without
re-estimating the structural settings of the original Pareto density.

If \code{V$Fit$countingKernels} is \code{NULL}, \code{V$kernels} are used as
the counting grid. Otherwise, the stored counting kernels are used.

For each omitted observation, one count is removed from every Pareto sphere
whose interval contains that observation. If the observation lies within one
Pareto radius of the lower or upper full-sample boundary, the corresponding
reflected contribution is removed as well.

The updated counts are normalized by trapezoidal integration over the fixed
evaluation grid. If the original fit used interpolation, the normalized
density is interpolated back to \code{V$kernels} with
\code{stats::approx()}. Non-finite interpolated values are replaced by zero.

For compatibility with older fit objects that do not store
\code{interpolationUsed}, interpolation is inferred from whether
\code{countingKernels} is non-\code{NULL}.
}
\note{
A numerical Dirac representation has no usable fit context and cannot be
re-evaluated by this function.
}
\examples{
\dontrun{
fit <- DataVisualizations::ParetoDensityEstimation(
  Data = rnorm(200),
  PlotIt = FALSE
)

fit_minus_one <- ParetoDensityEstimationFromFit(fit, Omit = 1)
fit_minus_many <- ParetoDensityEstimationFromFit(fit, Omit = c(1, 2, 2, 10))
}
}
\seealso{
\code{\link[DataVisualizations]{ParetoDensityEstimation}}
}
\keyword{distribution}
\keyword{nonparametric}
