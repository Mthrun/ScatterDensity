\name{PointsInPolygon}
\alias{PointsInPolygon}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
PointsInPolygon
}
\description{
Defines a Cls based on points in a given polygon.
}
\usage{
PointsInPolygon(Points, Polygon, PlotIt = FALSE, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{Points}{
    [1:n,1:2] xy cartesian coordinates of a projection
}
  \item{Polygon}{
    Numerical matrix of 2 columns defining a closed polygon
}
  \item{PlotIt}{
    TRUE: Plots marked points
}
  \item{\dots}{
    \code{BMUorProjected}: Default == FALSE, If TRUE assuming BestMatches of ESOM instead of Projected Points
    
    \code{main}: title of plot

    Further Plotting Arguments,xlab etc used in \code{\link[DataVisualizations]{Classplot}}

    }
}
\details{
We assume that polygon is closed, i.e., that the last point connects to the fist point
}
\value{
Numerical classification vector Cls with 1 = outside polygon and 2 = inside polygon
}

\author{
Michael Thrun
}

\seealso{
\code{\link[DataVisualizations]{Classplot}}
}

\examples{
XY=cbind(runif(80,min = -1,max = 1),rnorm(80))
#closed polygon
polymat <- cbind(x = c(0,1,1,0), y = c(0,0,1,1))
#takes sometimes more than 5 sec
\donttest{
Cls=PointsInPolygon(XY,polymat,PlotIt = TRUE)
}
}

\keyword{Gate}
\concept{Manual Gating}

