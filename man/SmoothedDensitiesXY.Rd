\name{SmoothedDensitiesXY}
\alias{SmoothedDensitiesXY}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Smoothed Densities X with Y
}
\description{
Density is the smothed histogram density at [X,Y] of [Eilers/Goeman, 2004]
}
\usage{
SmoothedDensitiesXY(X, Y, nbins, lambda, Xkernels, Ykernels,

Compute="Cpp", PlotIt = FALSE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{X}{
Numeric vector [1:n], first feature (for x axis values)
}
  \item{Y}{
Numeric vector [1:n], second feature (for y axis values),
}
  \item{nbins}{
Optional, number of bins, default is 250, if scalar nbins=nx=ny otherwise c(nx,ny) vector of length 2
}
  \item{lambda}{
Optional, smoothing factor used by the density estimator, default is lambda = 20
}
  \item{Xkernels}{
Optional, bin kernels in x direction are given, if not set nx bins in range of data are generated with \code{\link{seq}}
}
  \item{Ykernels}{
Optional, bin kernels y direction are given, if not set ny bins in range of data are generated with \code{\link{seq}}
}
  \item{Compute}{
Either \code{R}, \code{Cpp} or \code{Parallel}, defines how the heavy-load computation is performed
}
  \item{PlotIt}{
FALSE: no plotting, TRUE: simple plot
}
}
\details{
\code{lambda} has to chosen by the user and is a sensitive parameter and a \code{lambda = 20} roughly means that the smoothing is over 20 bins around a given point.
}
\value{
 List of:
  \item{Densities}{numeric vector [1:n] is the smothed density in 3D  }
  \item{Xkernels}{numeric vector [1:nx], nx defined by \code{nbins}, such that mesh(Xkernels,Ykernels,F) form the ( not NaN) smothed densisties}
  \item{Ykernels}{numeric vector [1:ny], nx defined by \code{nbins}, such that mesh(Xkernels,Ykernels,F) form the ( not NaN) smothed densisties}
  \item{GridDensity}{matrix [1:nx,1:ny] beeing the smoothed 2D histogram }
  \item{Points2GridInd}{an index such that \code{Densities = GridDensity[Points2GridInd]}}

}
\references{
[Eilers/Goeman, 2004]  Eilers, P. H., & Goeman, J. J.: Enhancing scatterplots with smoothed densities, Bioinformatics, Vol. 20(5), pp. 623-628.DOI: \doi{10.1093/bioinformatics/btg454}, 2004.
}
\author{
Michael Thrun
}

\examples{
\dontshow{
c1 = rnorm(n=40)
c2 = rnorm(n=40,1,2)
densities = SmoothedDensitiesXY(c1,c2)$GridDensity
}


\donttest{
if(requireNamespace("DataVisualizations")){
data("ITS",package = "DataVisualizations")
data("MTY",package = "DataVisualizations")
Inds=which(ITS<900&MTY<8000)
V=SmoothedDensitiesXY(ITS[Inds],MTY[Inds])
}else{
#sample random data
ITS=rnorm(1000)
MTY=rnorm(1000)
V=SmoothedDensitiesXY(ITS,MTY)
}
}

}
\keyword{scatter}
\keyword{SDH}
\concept{Density Estimation}
\concept{scatter plot}
\concept{scatter density plot}
