\name{InteractiveGate}
\alias{InteractiveGate}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Interactive manual gating of scatter pplots
}
\description{
Allows to draw a polygon around a specifc area in a two-dimensional scatter plot
}
\usage{
InteractiveGate(Data, GateVars, Gatename = "Leukocytes", PlotIt = FALSE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{Data}{
[1:n,1:d] numerical matrix of n cases and d variables. each variable has a name accessible via colnames(Data)
}
  \item{GateVars}{
[1:2] character vector of columns that are selected
}
  \item{Gatename}{
name of the gate
}
  \item{PlotIt}{
TRUE: plots the result
}
}
\details{
Only polygon gates can be used in the interactive tool.
}
\value{
a LIST,
\item{Polygon }{[1:p,1:2] numerical matrix of coordinates for variables named by\code{GateVars} defining the polygon  }
\item{Index }{indices of data points within the polygon}
\item{DataInGate }{[1:m,1:d] numerical matrix of n cases and d variables of data points within the polygon}
}

\author{
Michael Thrun
}


\seealso{
\code{\link[flowGate]{gs_gate_interactive}}
}
\examples{
Disk=""
#path=ReDi("AutoGating/09Originale",Disk)
#V=ReadLRN("110001_T1_NB_d11_N33k",path)
#Data=V$Data
#V=InteractiveGate(Data,GateVars = c( "SS", "FS"),PlotIt = T)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory (show via RShowDoc("KEYWORDS")):
\keyword{ Gate }
\keyword{ Gating }
% Use only one keyword per line.
% For non-standard keywords, use \concept instead of \keyword:
\concept{Manual Gating}

