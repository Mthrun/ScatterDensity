\name{Mask2Polygon}
\alias{Mask2Polygon}
\title{Construct a Convex Geometry from a Logical Population Mask}
\description{
Constructs a convex two-dimensional geometry for observations selected by
a logical mask.

The mask defines membership in a population, cluster, or hierarchy node,
while \code{GateVars} defines the two variables used for its geometric
representation.
}
\usage{
Mask2Polygon(
  Data,
  Mask,
  ClusterNo,
  GateVars,
  PlotIt = FALSE
)
}
\arguments{
  \item{Data}{
    Matrix or data frame containing observations in rows and variables in
    columns.
  }

  \item{Mask}{
    Logical vector of length \code{nrow(Data)} defining which observations
    belong to the population.

    \code{TRUE} selects an observation and \code{FALSE} excludes it.
    Missing values in \code{Mask} are treated as \code{FALSE}.
  }

  \item{ClusterNo}{
    Scalar identifier of the population, cluster, or hierarchy node.

    The value is passed to \code{\link{Cluster2Polygon}} and included in the
    resulting geometry representation.
  }

  \item{GateVars}{
    Vector defining the two columns of \code{Data} used as X and Y
    coordinates.

    Elements may be column names or column indices. At least two variables
    must be supplied. If more than two variables are supplied, only the first
    two are used and a warning is issued.

    The two selected variables must be different.
  }

  \item{PlotIt}{
    Logical. If \code{TRUE}, the observations selected by \code{Mask} and
    the resulting convex geometry are plotted by
    \code{\link{Cluster2Polygon}}. Defaults to \code{FALSE}.
  }
}
\value{
List V with following elements

Polygon a numerical matrix/data.frame with X and Y,
ClusterNo, Ratio-
An open polygon is closed internally.
}
\details{
Population membership is determined exclusively by \code{Mask}. For example,
a hierarchical child population may be defined as:

\preformatted{
ChildMask <- ParentMask & LocalGateRule
}

In contrast, \code{GateVars} only determines the local two-dimensional
coordinate system in which the selected population is represented.

Missing values in \code{Mask} are replaced by \code{FALSE}. If the resulting
mask does not select any observations, the function stops with an error.

The selected rows and the two gate variables are extracted from
\code{Data}. A local cluster vector is then constructed by assigning
\code{ClusterNo} to every selected observation.

The resulting data are passed to:

\preformatted{
Cluster2Polygon(
  Data = GateData,
  Cls = LocalCls,
  ClusterNo = ClusterNo,
  Ratio = 1,
  PlotIt = PlotIt
)
}

Consequently, \code{Mask2Polygon} always constructs a convex hull or, for a
degenerate set of selected observations, the corresponding point or line
geometry.

The geometric hull and the original logical mask are not necessarily
equivalent. In particular, observations that were not selected by
\code{Mask} may geometrically lie inside the convex hull of the selected
observations.

The returned geometry can be used as input for downstream geometric gating
functions such as \code{PointsInPolygon} or \code{PolygonGate}.
}
\examples{
Data <- matrix(
  c(
    0, 0,
    1, 0,
    1, 1,
    0, 1,
    5, 5
  ),
  ncol = 2,
  byrow = TRUE
)

colnames(Data) <- c("X", "Y")

Mask <- c(TRUE, TRUE, TRUE, TRUE, FALSE)

Polygon <- Mask2Polygon(
  Data = Data,
  Mask = Mask,
  ClusterNo = 1,
  GateVars = c("X", "Y")
)

Polygon

\dontrun{
# Optional visual inspection
Mask2Polygon(
  Data = Data,
  Mask = Mask,
  ClusterNo = 1,
  GateVars = c("X", "Y"),
  PlotIt = TRUE
)
}
}
\author{
Michael Thrun, IAP-GmbH Intelligent Analytics Projects
}
\seealso{
\code{\link{Cluster2Polygon}},
\code{\link{PointsInPolygon}},
\code{\link{PolygonGate}}
}
\keyword{cluster}
\keyword{geometry}
\keyword{utilities}