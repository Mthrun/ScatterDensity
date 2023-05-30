\name{DDCAL}
\alias{DDCAL}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Density Distribution Cluster Algorithm of [Lux and Rinderle-Ma, 2023].
}
\description{
DDCAL is a clustering-algorithm for one-dimensional data, which heuristically finds clusters to evenly distribute the data points in low variance clusters.
}
\usage{
DDCAL(data, nClusters, minBoundary = 0.1, maxBoundary = 0.45,

numSimulations = 20, csTolerance = 0.45, csToleranceIncrease = 0.5)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{data}{
  [1:n] Numeric vector, with the data values
}
  \item{nClusters}{
  Scalar, number of clusters to be found
}
  \item{minBoundary}{
  Scalar, in the range (0,1), gives the lower boundary (in percent), for the simulation. Default is 0.1
}
  \item{maxBoundary}{
  Scalar, in the range (0,1), gives the upper boundary (in percent), for the simulation. Default is 0.45
}
  \item{numSimulations}{
  Scalar, number of simulations/iterations of the algorithm
}
  \item{csTolerance}{
  Scalar, in the range (0,1). Gives cluster size tolerance factor. The necessary cluster size is defined by (dataSize/nClusters - dataSize/nClusters * csTolerance). Default is 0.45
}
  \item{csToleranceIncrease}{
  Scalar, in the range (0,1), gives the procentual increase of the csTolerance-factor, if some clusters did not reach the necessary size. Default is 0.5
}
}
\details{
DDCAL creates a evenly spaced division of the min-max-normalized data from minBoundary to maxBoundary. Those divisions will be used as boundaries. The first initial clusters will be the data from min(data) to minBoundary and maxBoundary to max(data). The clusters will be extended to neighboring points, as long as the standard deviations of the clusters will be reduced. A potential clusters will be used, if they have the necessary size, given as (dataSize/nClusters - dataSize/nClusters * csTolerance). If both clusters can be used, the left cluster (which is the cluster from min(data) to minBoundary or above) is preferred. If no clusters can be found with the necessary size, then the csTolerance-factor and with it the necessary cluster size will be lowered.
If a clusters is used, the next boundaries are found, which are not in the already existing clusters and the procedure is repeated with the not already clustered data, until all points are assigned to clusters.


If a matrix is given as input data, the first column of the matrix will be used as data for the clustering

Non-finite values will not be clustered, but instead will get the cluster label \code{NaN}.

The algorithm is not garantueed to produce the given number of clusters, given in nClusters. The found number of clusters can be lower, depending on the data and input parameters.

}
\value{
\item{labels}{[1:n] Numeric vector, containing the labels for the input data points}

}
\references{
[Lux and Rinderle-Ma, 2023] Lux, M., Rinderle-Ma, S.: DDCAL: Evenly Distributing Data into Low Variance Clusters Based on Iterative Feature Scaling; Springer Journal of Classification, Vol. 40, pp. 106-144, DOI: \doi{10.1007/s00357-022-09428-6}, 2023.
}
\author{
Luca Brinkmann
}

\examples{
\dontshow{
c1 = rnorm(n=40)
c2 = rnorm(n=40,1,2)
densities = SmoothedDensitiesXY(c1,c2)$Densities
labels = ScatterDensity::DDCAL(densities, 2)
}
\donttest{
# Load data
if(requireNamespace("FCPS")){
data(EngyTime, package = "FCPS")
engyTimeData = EngyTime$Data
c1 = engyTimeData[,1]
c2 = engyTimeData[,2]
}else{
c1 = rnorm(n=4000)
c2 = rnorm(n=4000,1,2)
}
# Calculate Densities

densities = SmoothedDensitiesXY(c1,c2)$Densities
# Use DDCAL to cluster the densities
labels = DDCAL(densities, 9)

# Plot Densities according to labels
my_colors = c("#000066", "#3333CC", "#9999FF", "#00FFFF", "#66FF33",
                       "#FFFF00", "#FF9900", "#FF0000", "#990000")
labels = as.factor(labels)
df = data.frame(c1, c2, labels)

if(requireNamespace("ggplot2")){
ggplot2::ggplot(df, ggplot2::aes(c1, c2, color = labels)) +
  ggplot2::geom_point() +
  ggplot2::scale_color_manual(values = my_colors)
}
}
}

\keyword{Clustering}
\keyword{Density}
\keyword{Cluster}
\concept{One-Dimensional Clustering}
\concept{Iterative Clustering}
\concept{Low-Variance Clustering}

% Use only one concept per line.
