\name{Cluster2Polygon}
\alias{Cluster2Polygon}
\title{Construct a Convex or Concave Boundary for One Cluster}
\description{
Constructs a two-dimensional geometric boundary for one selected cluster using
the first two columns of \code{Data}.

With \code{Ratio = 1}, the function returns the convex hull. With
\code{Ratio < 1}, it requests a concave hull using
\code{sf::st_concave_hull()} when the selected cluster contains at least three
non-collinear hull points.
}
\usage{
Cluster2Polygon(
  Data,
  Cls,
  ClusterNo,
  Ratio = 1,
  PlotIt = FALSE
)
}
\arguments{
  \item{Data}{
    Matrix or data frame containing observations in rows and variables in
    columns. At least two columns are required. Only the first two columns are
    used for geometry construction, and both must be numeric.
  }
  \item{Cls}{
    Vector of cluster assignments. Its length must equal
    \code{nrow(Data)}. Missing cluster labels are ignored when selecting the
    requested cluster.
  }
  \item{ClusterNo}{
    Scalar cluster identifier specifying which observations are used.
    Comparison with \code{Cls} is performed after conversion to character.
  }
  \item{Ratio}{
    Numeric scalar in \code{[0,1]} controlling the requested geometry.
    \code{Ratio = 1} returns the convex hull; \code{Ratio < 1} requests a
    concave hull using \code{sf::st_concave_hull()}. Smaller values generally
    allow a more concave boundary. Default is \code{1}.
  }
  \item{PlotIt}{
    Logical scalar. If \code{TRUE}, the finite unique observations belonging
    to the selected cluster and the resulting geometry are plotted with base R
    graphics. Default is \code{FALSE}.
  }
}
\value{
A list containing exactly four elements:

\item{Vertices}{
  Numeric matrix with columns \code{VertexNo}, \code{X}, and \code{Y}.
  Each row is one vertex of the returned geometry. Convex polygon output is
  explicitly closed by repeating the first polygon vertex at the end.
}

\item{ClusterNo}{
  The cluster identifier supplied to the function.
}

\item{Ratio}{
  The ratio used to construct the geometry.
}

\item{GeometryType}{
  Character scalar describing the geometry returned for the selected cluster.
  Possible values include \code{"POINT"}, \code{"LINESTRING"}, and
  \code{"POLYGON"}.
}
}
\details{
Observations are selected using

\code{!is.na(Cls) & as.character(Cls) == as.character(ClusterNo)}.

Only the first two columns of \code{Data} are used. Rows containing
\code{NA}, \code{NaN}, \code{Inf}, or \code{-Inf} in either coordinate are
removed before constructing the geometry. Duplicate coordinate pairs are also
removed.

A convex hull is always determined first and is used to detect degenerate
clusters:
\itemize{
  \item one hull point gives \code{"POINT"};
  \item two hull points give \code{"LINESTRING"};
  \item three or more hull points give \code{"POLYGON"}.
}

For one or two unique finite observations, the observations themselves define
the hull. Otherwise hull indices are calculated with \code{grDevices::chull()}.

If \code{Ratio = 1}, the convex geometry is returned directly. Convex polygon
paths are explicitly closed by appending the first vertex.

If \code{Ratio < 1} and at least three non-collinear hull points are available,
package \pkg{sf} is required. The selected coordinates are converted to a
MULTIPOINT geometry and passed to

\code{sf::st_concave_hull(..., ratio = Ratio, allow_holes = FALSE)}.

If \code{Ratio < 1} but the cluster is degenerate, a concave polygon is not
defined. The corresponding POINT or LINESTRING is returned instead and a
warning is issued. In this case \pkg{sf} is not required.

When \code{PlotIt = TRUE}, selected observations are plotted and the resulting
geometry is overlaid in red.
}
\note{
The function validates that \code{Data} has at least two columns,
\code{length(Cls) == nrow(Data)}, \code{ClusterNo} contains exactly one
non-missing value, \code{Ratio} is finite and lies in \code{[0,1]}, and the
first two columns used for geometry construction are numeric.
}
\examples{
Data <- data.frame(
  X = c(0, 1, 1, 0, 3, 4, 4, 3),
  Y = c(0, 0, 1, 1, 0, 0, 1, 1)
)

Cls <- c(1, 1, 1, 1, 2, 2, 2, 2)

poly <- Cluster2Polygon(
  Data = Data,
  Cls = Cls,
  ClusterNo = 1,
  Ratio = 1
)

poly$Vertices
poly$GeometryType

\dontrun{
concave <- Cluster2Polygon(
  Data = Data,
  Cls = Cls,
  ClusterNo = 1,
  Ratio = 0.5,
  PlotIt = TRUE
)
}
}
\seealso{
\code{\link[grDevices]{chull}}
}
\author{
Michael Thrun
}
\keyword{cluster}
\keyword{geometry}
