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
nbins= nxy      => the nr of bins in x and y is nxy
nbins = c(nx,ny) => the nr of bins in x is nx and for y is ny
}
  \item{nbins}{
 number of bins, nbins =200 (default)
}
  \item{lambda}{
smoothing factor used by the density estimator or c()
default: lambda = 20 which roughly  means that the smoothing is over 20 bins around a given point.
}
  \item{Xkernels}{
bin kernels in x direction are given
}
  \item{Ykernels}{
bin kernels y direction are given
}
  \item{Compute}{
Either \code{R}, \code{Cpp} or \code{Parallel}, defines where heavy-load computation is performed
}
  \item{PlotIt}{
FALSE: no plotting, TRUE: simple plot
}
}
\details{
lambda has to chosen by the user and is a sensitive parameter.
}
\value{
 List of:
   \item{Densities}{numeric vector [1:n] is the smothed density in 3D  }
  \item{Xkernels}{numeric vector [1:nx], nx defined by \code{nbins}, such that mesh(Xkernels,Ykernels,F) form the ( not NaN) smothed densisties}
  \item{Ykernels}{numeric vector [1:ny], nx defined by \code{nbins}, such that mesh(Xkernels,Ykernels,F) form the ( not NaN) smothed densisties}
  \item{GridDensity}{matrix [1:nx,1:ny] beeing the smoothed 2D histogram }
  \item{Points2GridInd}{an index such that \code{Densities = GridDensity[Points2GridInd]}}
   \item{MaxInGrid}{scalar, maximum number of points in grid that was used to normalize \code{GridDensity} }
  

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
