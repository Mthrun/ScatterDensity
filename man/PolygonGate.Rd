\name{PolygonGate}
\alias{PolygonGate}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
PolygonGate
}
\description{
A specific Gate defined by xy coordinates that result in a closed polygon is applied to the flowcytometry data.
}
\usage{
PolygonGate(Data, Polygon, GateVars,  PlotIt = FALSE, PlotSampleSize = 1000)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{Data}{
numerical matrix n x d
}
  \item{Polygon}{
numerical marix of two columns defining the coordiantes of the polygon. polygon assumed to be closed, i.e.,last coordinate connects to first coordinate.
}
  \item{GateVars}{
vector, either column index in Data of X and Y coordinates of gate or its variable names as string
}


  \item{PlotIt}{
if TRUE: plots a sample of data in the two selected variables and marks point inside the gate as yellow and outside as magenta
}
  \item{PlotSampleSize}{
size pof the plottet sample
}
}
\details{
Gates are alwaxs two dimensional, i.e., require two filters, although all dimensions of data are filted by the gates. Only high-dimensional points inside the polygon (gate) are given back
}
\value{
list of
\item{DataInGate }{m x d numerical matrix with m<=n of data points inside the gate}
\item{InGateInd }{index of length m for the datapoints in original matrix}

}
\author{
Michael Thrun
}
\note{
if GateVars is not found a text is given back which will state this issue
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
\code{\link{PointsInPolygon}}
}
\examples{
Data <- matrix(runif(1000), ncol = 10)
colnames(Data)=paste0("GateVar",1:ncol(Data))
poly <- cbind(x = c(0.2,0.5,0.8), y = c(0.2,0.8,0.2))
#set PlotIt TRUE for understanding the example

\donttest{
#Select index
V=PolygonGate(Data,poly,c(5,8),PlotIt=FALSE,100)

#select var name
V=PolygonGate(Data,poly,c("GateVar5","GateVar8"),PlotIt=FALSE,100)
}
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory (show via RShowDoc("KEYWORDS")):
 \keyword{Gate}
% \keyword{ ~kwd2 }
% Use only one keyword per line.
% For non-standard keywords, use \concept instead of \keyword:
\concept{Manual Gating }
% \concept{ ~cpt2 }
% Use only one concept per line.
