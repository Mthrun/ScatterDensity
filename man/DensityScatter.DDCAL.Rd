\name{DensityScatter.DDCAL}
\alias{DensityScatter.DDCAL}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Scatter density plot using DDCAL
}
\description{
Density estimation (PDE) [Ultsch, 2005] or "SDH" [Eilers/Goeman, 2004] used for a scatter density plot, with clustering of densities with DDCAL [Lux/Rinderle-Ma, 2023]
}
\usage{
DensityScatter.DDCAL(X, Y, xlab, ylab, SDHorPDE = TRUE, Plotter = "native",

Silent = FALSE, Marginals = FALSE, pch = 10, Size = 1, BW = TRUE,

PDEsample = 5000, lwd = 2, na.rm=TRUE,Polygon, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{X}{
Numeric vector [1:n], first feature (for x axis values)
}
  \item{Y}{
Numeric vector [1:n], second feature (for y axis values)
}
  \item{xlab}{
String, title of the x axis. Default: "X",   see \code{plot()} function
}
  \item{ylab}{
String, title of the y axis. Default: "Y",   see \code{plot()} function
}
  \item{SDHorPDE}{
(Optional) Boolean, if TRUE SDH is used to calculate density, if FALSE PDE is used
}
  \item{Plotter}{
(Optional) String, name of the plotting backend to use. Possible values are: "\code{native}" and "\code{ggplot2}", default is "\code{native}""
}
  \item{Silent}{
(Optional) Boolean, if TRUE no messages will be printed, default is FALSE
}
  \item{Marginals}{
(Optional) Boolean, if TRUE the marginal distributions of X and Y will be plotted together with the 2D density of X and Y. Default is FALSE
}
  \item{pch}{
(Optional) Scalar or character. Indicates the shape of data points, see \code{plot()} function or the \code{shape} argument in ggplot2. Default is 10
}
  \item{Size}{
(Optional) Scalar, size of data points in plot, default is 1
}
  \item{BW}{
(Optional) Boolean, if TRUE ggplot2 will use a white background, if FALSE the typical ggplot2 backgournd is used. Not needed if "\code{native}" as Plotter is used. Default is TRUE
}
  \item{PDEsample}{
(Optional) Scalar, Sample size for PDE and/or for ggplot2 plotting. Default is 5000
}
  \item{lwd}{
(Optional) Scalar, thickness of the lines used for the marginal distributions (only needed if \code{Marginals=TRUE}), see \code{plot()}. Default = 2
}
  \item{na.rm}{
(Optional) Boolean, if TRUE non finite values will be removed
}
  \item{Polygon}{
(Optional) [1:p,1:2] numeric matrix that defines for x and y coordinates a polygon in magenta, works only for native method
}
  \item{\dots}{
Further plot arguments
}
}
\details{
The \code{DensityScatter.DDCAL} function generates the density of the xy data as a z coordinate. Afterwards xyz will be plotted as a contour plot. It assumens that the cases of x and y are mapped to each other meaning that a \code{cbind(x,y)} operation is allowed.
The colors for the densities in the contour plot are calculated with DDCAL, which produces clusters to evenly distribute the densities in low variance clusters.

In the case of "\code{native}" as Plotter, the handle returns \code{NULL} because the basic R functon \code{plot}() is used
}
\value{
If "\code{ggplot2}" as Plotter is used, the ggobj is returned

}
\references{
[Ultsch, 2005]  Ultsch, A.: Pareto density estimation: A density estimation for knowledge discovery, In Baier, D. & Werrnecke, K. D. (Eds.), Innovations in classification, data science, and information systems, (Vol. 27, pp. 91-100), Berlin, Germany, Springer, 2005.

[Eilers/Goeman, 2004]  Eilers, P. H., & Goeman, J. J.: Enhancing scatterplots with smoothed densities, Bioinformatics, Vol. 20(5), pp. 623-628. 2004.

[Lux/Rinderle-Ma, 2023]  Lux, M. & Rinderle-Ma, S.: DDCAL: Evenly Distributing Data into Low Variance Clusters Based on Iterative Feature Scaling, Journal of Classification vol. 40, pp. 106-144, 2023.

}
\author{
Luca Brinkmann, Michael Thrun
}
\note{
Support for plotly will be implemented later
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\examples{
\dontshow{
# Create two bimodial distributions
x1=rnorm(n = 75,mean = 0,sd = 1)
y1=rnorm(n = 75,mean = 0,sd = 1)
x2=rnorm(n = 75,mean = 2.5,sd = 1)
y2=rnorm(n = 75,mean = 2.5,sd = 1)
x=c(x1,x2)
y=c(y1,y2)

DensityScatter.DDCAL(x, y, Marginals = TRUE)
}

\donttest{
# Create two bimodial distributions
x1=rnorm(n = 7500,mean = 0,sd = 1)
y1=rnorm(n = 7500,mean = 0,sd = 1)
x2=rnorm(n = 7500,mean = 2.5,sd = 1)
y2=rnorm(n = 7500,mean = 2.5,sd = 1)
x=c(x1,x2)
y=c(y1,y2)

DensityScatter.DDCAL(x, y, Marginals = TRUE)
}
}
\keyword{scatter}
\keyword{PDE}
\keyword{SDH}
\concept{DDCAL}
\concept{Density Estimation}
\concept{scatter plot}
\concept{scatter density plot}