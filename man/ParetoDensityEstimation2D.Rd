\name{ParetoDensityEstimation2D}
\alias{ParetoDensityEstimation2D}
\title{Paired Two-Dimensional Pareto Density Extension}
\description{
Constructs a two-dimensional density estimate for two paired continuous
variables. The estimator combines published one-dimensional Pareto Density
Estimates with a paired rectangular joint-density seed and iterative
proportional fitting.

The resulting joint distribution preserves the observed pairing between
\code{X} and \code{Y} while matching the two one-dimensional Pareto Density
Estimate marginals.
}
\usage{
ParetoDensityEstimation2D(X, Y)
}
\arguments{
  \item{X}{
    Numeric vector containing observations of the first continuous variable.
  }

  \item{Y}{
    Numeric vector containing observations of the second continuous variable.
    Must have the same length as \code{X}.
  }
}
\details{
The function performs the following steps:

\enumerate{
  \item Estimate the one-dimensional Pareto densities of \code{X} and \code{Y}
  using \code{DataVisualizations::ParetoDensityEstimation()}.

  \item Extract the Pareto radii, kernel locations, and density estimates
  returned by the one-dimensional estimator.

  \item If a one-dimensional kernel grid contains more than 256 locations,
  interpolate the published density to a grid of 256 locations. This limits
  quadratic memory growth and is not an additional smoothing operation.

  \item Convert the one-dimensional density estimates to normalized probability
  masses using trapezoidal quadrature.

  \item Construct a paired rectangular joint-density seed:
  \deqn{
  q(s,t) =
  \sum_i I\{|X_i-s| \leq r_X\}
         I\{|Y_i-t| \leq r_Y\}.
  }

  \item Convert the joint seed to probability mass and add a very small positive
  support floor.

  \item Apply iterative proportional fitting so that the row and column
  marginals agree with the two one-dimensional Pareto Density Estimate
  probability masses.

  \item Convert the fitted probability mass back to a density and calculate
  density values at the observed pairs using bilinear interpolation.
}

The IPF algorithm stops when the maximum marginal error is less than
\code{1e-8}, or after 500 iterations. A warning is issued if the tolerance is
not reached.

This function is an explicit paired rectangular two-dimensional extension
constrained by published one-dimensional Pareto Density Estimates. It is not
claimed to be a published two-dimensional hypersphere Pareto Density
Estimator.
}
\value{
A list with components:

\item{Densities}{
  Numeric vector containing the bilinearly interpolated joint-density value at
  each observed pair \code{(X[i], Y[i])}.
}

\item{GridDensity}{
  Numeric matrix containing the estimated joint density on the rectangular
  grid. Rows correspond to \code{Xkernels}; columns correspond to
  \code{Ykernels}.
}

\item{Xkernels}{
  Numeric vector containing the kernel-grid locations for \code{X}.
}

\item{Ykernels}{
  Numeric vector containing the kernel-grid locations for \code{Y}.
}
}
\note{
At least 50 paired observations and at least 12 unique values in each variable
are required.

All observations must be finite. Remove missing, infinite, or non-finite paired
observations before calling this function.

The suggested package \pkg{DataVisualizations} is required.
}
\examples{
\dontrun{
set.seed(1)

X <- rnorm(200)
Y <- 0.7 * X + rnorm(200, sd = 0.6)

fit <- ParetoDensityEstimation2D(X, Y)

length(fit$Densities)
dim(fit$GridDensity)
length(fit$Xkernels)
length(fit$Ykernels)
}
}
\seealso{
\code{\link{Mutualinfo_pde_2d}},
\code{\link[DataVisualizations]{ParetoDensityEstimation}}
}
\keyword{distribution}
\keyword{multivariate}
\keyword{nonparametric}
