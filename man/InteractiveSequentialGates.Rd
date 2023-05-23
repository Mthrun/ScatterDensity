\name{InteractiveSequentialGates}
\alias{InteractiveSequentialGates}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Interactive Sequential Gates
}
\description{
Performs several consecutive gates interactively on the same data. Only polygon gates can be used in the interactive tool.
}
\usage{
InteractiveSequentialGates(Data, GateVarsMat, Gatenames, PlotIt = FALSE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{Data}{
[1:n,1:d] numerical matrix of n cases and d variables. each variable has a name accessible via colnames(Data)
}
  \item{GateVarsMat}{
[1:2,1:c] character matrix of columns that are selected, each rwo represents one gate
}
  \item{Gatenames}{
[1:c] character vector of names of the sequential gates
}
  \item{PlotIt}{
TRUE: plots the last gate on the prior of last data
}
}
\details{
In flow cytometry, the process of setting sequential or consecutive gates is known as interactive gating. Interactive gating involves analyzing and selecting specific populations of cells or particles based on their characteristics in a step-by-step manner.

Once the initial gate is set, you can further refine the analysis by examining the population within that gate. By assessing the patterns and characteristics of the cells within the initial gate, you can determine which parameters are relevant for subsequent gating.

Based on the analysis of the first gate, you can then create a second gate or subsequent gates to isolate more specific subsets of cells. This process involves creating new scatter plots based on the relevant parameters identified earlier. The successive gates are usually set based on the characteristics observed within the previously gated populations.
}
\value{
a LIST,
\item{Polygons }{[1:c] list named by Gatenames, each element is a [1:p,1:2] numerical matrix of coordinates for variables named by\code{GateVars} defining the polygon  }
\item{Indices }{[1:c] list named by Gatenames, each element defines the indices of data points within the polygon with regards to the data of the prior polygon (gate)}
\item{DataInGates }{[1:c] list named by Gatenames, each element is a [1:m,1:d] numerical matrix of n cases and d variables of data points within the current polygon}
}

\author{
Michael Thrun
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
\code{\link[flowGate]{gs_gate_interactive}}
}
\examples{
Disk=""
#path=ReDi("AutoGating/09Originale",Disk)
#read data
#filename="110001_T1_NB_d11_N33k"
#V=dbt.DataIO::ReadLRN(filename,path)
#Data=V$Data
#find out names of given parameters
#colnames(Data)
#select parameters, each two define one gate
#GateVarsMat=rbind(c("SS", "FS"),c("CD45", "CD19"))
#name the gates, the first one is overall data, every next one is a consecutive gate
#base on the gated information of the previous gate
#Gatenames=c("Lymphocytes","Leukocytes")
#apply interactive function
#GateInfo=ScatterDensity::InteractiveSequentialGates(Data,GateVarsMat =GateVarsMat,Gatenames =Gatenames#,PlotIt = T)

}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory (show via RShowDoc("KEYWORDS")):
\keyword{ Gate }
\keyword{ Gating }
% Use only one keyword per line.
% For non-standard keywords, use \concept instead of \keyword:
\concept{Manual Gating}
