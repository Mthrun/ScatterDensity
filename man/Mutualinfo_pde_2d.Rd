\name{Mutualinfo_pde_2d}
\alias{Mutualinfo_pde_2d}
\alias{Mutualinfo_pde_2d_jackknife}
\title{Mutual Information from a Two-Dimensional Pareto Density Estimate}
\description{
Estimates mutual information between two continuous variables from the joint
density returned by \code{\link{ParetoDensityEstimation2D}}.
}
\usage{
Mutualinfo_pde_2d(x, y, base = 2)
}
\arguments{
  \item{x}{
    Numeric vector containing observations of the first continuous variable.
  }

  \item{y}{
    Numeric vector containing observations of the second continuous variable.
    It must satisfy the input requirements of
    \code{\link{ParetoDensityEstimation2D}}.
  }

  \item{base}{
    Numeric scalar defining the logarithm base. It must be finite, positive,
    and different from 1. Default is \code{2}, yielding mutual information in
    bits.
  }
}
\details{
The function first estimates the joint density using
\code{\link{ParetoDensityEstimation2D}}.

The density grid is converted to a probability-mass matrix using trapezoidal
quadrature:

\deqn{
p_{ij} =
f(x_i,y_j) w_i^{(X)} w_j^{(Y)}.
}

The mass matrix is normalized to sum to one. Marginal probability masses are
then calculated from its row and column sums.

Mutual information is evaluated as:

\deqn{
I(X;Y) =
\sum_{i,j} p_{ij}
\log_b\left(
\frac{p_{ij}}{p_{i\cdot}p_{\cdot j}}
\right),
}

where \eqn{b} is the logarithm base supplied through \code{base}. Only cells
with positive joint probability mass are included.

A small negative result caused by floating-point error is replaced by zero.
}
\value{
Numeric scalar containing the estimated mutual information.

With the default \code{base = 2}, the result is expressed in bits. With
\code{base = exp(1)}, it is expressed in nats.
}
\note{
The estimator inherits all sample-size, uniqueness, and dependency requirements
of \code{\link{ParetoDensityEstimation2D}}.

In particular, at least 50 paired observations and at least 12 unique values
per variable are required.
}
\examples{
\dontrun{
set.seed(1)

x <- rnorm(200)
y <- x + rnorm(200, sd = 0.5)

Mutualinfo_pde_2d(x, y)

# Result in natural-logarithm units
Mutualinfo_pde_2d(x, y, base = exp(1))
}
}
\seealso{
\code{\link{ParetoDensityEstimation2D}}
}
\keyword{distribution}
\keyword{multivariate}
\keyword{nonparametric}