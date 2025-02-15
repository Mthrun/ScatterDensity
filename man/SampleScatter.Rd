\name{SampleScatter}
\alias{SampleScatter}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
takes a sample for a scatter plot
}
\description{
Given 2D points having X and Y coordinates takes a sample, such that these points are
is optimally visualized if a plot function is called.
}
\usage{
SampleScatter(X, Y, ThresholdPoints = 20,

DensityThreshold, nbins = 100,na.rm=TRUE, PlotIt = FALSE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{X}{
[1:n] x coordinate
}
  \item{Y}{
[1:n] y coordinate
}
  \item{ThresholdPoints}{
how many points per grid are allowed to overlap, default 20
}
  \item{DensityThreshold}{
SDH threshold for density computed with \code{SmoothedDensitiesXY}
}
  \item{nbins}{
number of bins in grid, default 100x100
}

  \item{na.rm}{Function may not work with non finite values. If these cases should be automatically removed, set parameter TRUE}
  
  \item{PlotIt}{
Plots the remaining points
}
}
\details{
"Optimally"" visualized in the sense that not too much point overap visually
}
\value{
SubsampleInd[1:m] indices of m points, m<n, that will be relevant for a optimal scatter plot
}

\author{
Michael Thrun
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
\code{\link{SmoothedDensitiesXY}}
}
\examples{
\dontshow{
c1 = rnorm(n=40)
c2 = rnorm(n=40,1,2)
sample_ind = SampleScatter(c1,c2)

}

\donttest{
if(requireNamespace("DataVisualizations")){
data("ITS",package = "DataVisualizations")
data("MTY",package = "DataVisualizations")
sample_ind=SampleScatter(ITS,MTY,PlotIt=TRUE)
}else{
#sample random data
ITS=rnorm(10000)
MTY=rnorm(10000)
sample_ind=SampleScatter(ITS,MTY,ThresholdPoints = 5)
del_ind=setdiff(1:length(ITS),sample_ind)
plot(ITS,MTY,type="p",pch=20,col="grey",main="Grey=full data, red=overlapping data points")
points(ITS[del_ind],MTY[del_ind],type="p",pch=20,col="red")
}
}


}

% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory (show via RShowDoc("KEYWORDS")):
\keyword{SampleScatter}
\keyword{Scatter}
% Use only one keyword per line.
% For non-standard keywords, use \concept instead of \keyword:
% \concept{ ~cpt1 }
\concept{Sample Scatter}
% Use only one concept per line.
