\name{PDEscatter}
\alias{PDEscatter}
\title{Scatter Density Plot}
\description{
Concept of Pareto density estimation (PDE) proposed for univsariate data by [Ultsch, 2005] and comparet to varius density estimation techniques by [Thrun et al., 2020] for univariate data is here applied for a scatter density plot. It was also applied in [Thrun and Ultsch, 2018] to bivariate data, but is not yet compared to other techniques.
}
\usage{
PDEscatter(x,y,SampleSize,

na.rm=FALSE,PlotIt=TRUE,ParetoRadius,Compute="Cpp", 

sampleParetoRadius, NrOfContourLines=20,Plotter='native',

DrawTopView = TRUE, xlab="X", ylab="Y", main="PDEscatter",
                              
xlim, ylim, Legendlab_ggplot="value")
}
\details{
The \code{PDEscatter} function generates the density of the xy data as a z coordinate. Afterwards xyz will be plotted either as a contour plot or a 3d plot. It assumens that the cases of x and y are mapped to each other meaning that a \code{cbind(x,y)} operation is allowed.
This function plots the PDE on top of a scatterplot. Variances of x and y should not differ by extreme numbers, otherwise calculate the percentiles on both first. If \code{DrawTopView=FALSE} only the plotly option is currently available. If another option is chosen, the method switches automatically there.

The method was succesfully used in [Thrun, 2018; Thrun/Ultsch 2018].

\code{PlotIt=FALSE} is usefull if one likes to perform adjustements like axis scaling prior to plotting with \pkg{ggplot2} or  \pkg{plotly}. In the case of "\code{native}"" the handle returns \code{NULL} because the basic R functon \code{plot}() is used
}
\arguments{
  \item{x}{Numeric vector [1:n], first feature (for x axis values)}
  \item{y}{Numeric vector [1:n], second feature (for y axis values)}
  \item{SampleSize}{Numeric m, positiv scalar, maximum size of the sample used for calculation. High values increase runtime significantly. The default is that no sample is drawn}
  \item{na.rm}{Function may not work with non finite values. If these cases should be automatically removed, set parameter TRUE}
  \item{ParetoRadius}{Numeric, positiv scalar, the Pareto Radius. If omitted (or 0), calculate by paretoRad.}
  \item{Compute}{
      Either \code{Cpp} or \code{Parallel}, defines if multicore computation should be performed
    }
   \item{sampleParetoRadius}{Numeric, positiv scalar, maximum size of the sample used for estimation of "kernel", should be significantly lower than SampleSize because requires distance computations which is memory expensive}
    \item{PlotIt}{\code{TRUE}: plots with function call
    
    \code{FALSE}: Does not plot, plotting can be done using the list element \code{Handle}
    
    \code{-1}: Computes density only, does not perfom any preperation for plotting meaning that \code{Handle=NULL}
    }

  \item{NrOfContourLines}{Numeric, number of contour lines to be drawn. 20 by default.}
  \item{Plotter}{String, name of the plotting backend to use. Possible values are: "\code{native}", "\code{ggplot}", "\code{plotly}"}
  \item{DrawTopView}{Boolean, True means contur is drawn, otherwise a 3D plot is drawn. Default: TRUE}
    \item{xlab}{String, title of the x axis. Default: "X",   see \code{plot()} function}
  \item{ylab}{String, title of the y axis. Default: "Y",   see \code{plot()} function}
    \item{main}{
    string, the same as "main" in \code{plot()} function
  }
  \item{xlim}{
  see \code{plot()} function
  }
    \item{ylim}{
    see \code{plot()} function
  }
    \item{Legendlab_ggplot}{String, in case of \code{Plotter="ggplot"} label for the legend. Default: "value"}
}
\note{
MT contributed with several adjustments
}

\value{
 List of:
    \item{X}{Numeric vector [1:m],m<=n, first feature used in the plot or the kernels used}
  \item{Y}{Numeric vector [1:m],m<=n, second feature used in the plot or the kernels used}
  \item{Densities}{Numeric vector [1:m],m<=n, Number of points within the ParetoRadius of each point, i.e. density information}
   \item{Matrix3D}{1:n,1:3] marix of x,y and density information}
  
  \item{ParetoRadius}{ParetoRadius used for PDEscatter}
  \item{Handle}{Handle of the plot object. Information-string if native R plot is used.}

}
\references{
[Thrun/Ultsch, 2018] Thrun, M. C., & Ultsch, A. : Effects of the payout system of income taxes to municipalities in Germany, in Papiez, M. & Smiech,, S. (eds.), Proc. 12th Professor Aleksander Zelias International Conference on Modelling and Forecasting of Socio-Economic Phenomena, pp. 533-542, Cracow: Foundation of the Cracow University of Economics, Cracow, Poland, 2018. 

[Ultsch, 2005]  Ultsch, A.: Pareto density estimation: A density estimation for knowledge discovery, In Baier, D. & Werrnecke, K. D. (Eds.), Innovations in classification, data science, and information systems, (Vol. 27, pp. 91-100), Berlin, Germany, Springer, 2005.

[Thrun et al., 2020] Thrun, M. C., Gehlert, T. & Ultsch, A.: Analyzing the Fine Structure of Distributions, PLoS ONE, Vol. 15(10), pp. 1-66, DOI \doi{10.1371/journal.pone.0238835}, 2020. 

}

\examples{
#taken from [Thrun/Ultsch, 2018]
if(requireNamespace("DataVisualizations")){
data("ITS",package = "DataVisualizations")
data("MTY",package = "DataVisualizations")
Inds=which(ITS<900&MTY<8000)
plot(ITS[Inds],MTY[Inds],main='Bimodality is not visible in normal scatter plot')

\donttest{
PDEscatter(ITS[Inds],MTY[Inds],xlab = 'ITS in EUR',

ylab ='MTY in EUR' ,main='Pareto Density Estimation indicates Bimodality' )
}
\dontshow{

x=rnorm(10)
y=rnorm(10)

PDEscatter(x,y,xlab = 'x in EUR',

ylab ='y in EUR' ,main='Pareto Density Estimation indicates Bimodality' )
}
}
}

\author{Felix Pape}

\keyword{scatter}
\keyword{PDE}

\concept{Density Estimation}
\concept{scatter plot}
\concept{scatter density plot}