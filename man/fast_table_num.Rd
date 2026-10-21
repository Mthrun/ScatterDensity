\name{fast_table_num}
\alias{fast_table_num}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
fast_table_num
}
\description{
Faster table computation in R compared to table() by omitting as.character and as.factor
}
\usage{
fast_table_num(x, y, edges_x, edges_y, redefine = TRUE,

byrow = FALSE, all.inside = FALSE, rightmost.closed = FALSE,

sort = FALSE, na.rm = FALSE, names = FALSE,extendOutput=FALSE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{
[1:n] numerical vector 
}
  \item{y}{
[1:n] numerical vector
}
  \item{edges_x}{
 Optional, [1:(\eqn{k_1}+1)] numerical  vector defining the specific borders in \code{x} default \code{unique(x)} for categorical scale
}
  \item{edges_y}{
 Optional, [1:(\eqn{k_2}+1)] numerical  vector defining the specific borders in \code{y}, default \code{unique(y)} for categorical scale
}
  \item{redefine}{
 Optional, boolean TRUE: resets counts in \code{y} direction in order from 1:\eqn{k_2} to \eqn{k_2}:1
}
  \item{byrow}{
 Optional, boolean, If FALSE (the default) the count matrix is filled by columns, otherwise the matrix is filled by rows.
}
  \item{all.inside}{
 Optional, boolean, if TRUE, the returned indices are coerced into 1,...,N-1, i.e., 0 is mapped to 1 and N to N-1
}
  \item{rightmost.closed}{
 Optional,  boolean, if TRUE, the rightmost interval, vec[N-1] .. vec[N] is treated as closed
}
  \item{sort}{
  Optional,  boolean, if TRUE, \code{edges_x}, \code{edges_y} are sorted non-decreasingly, NA/NaN are gthen ignored
}
  \item{na.rm}{
 Optional, boolean, if TRUE, only complete observations are taken into account
}
  \item{names}{
 Optional, boolean, if TRUE, output matrix is named by \code{edges_x}[1:\eqn{k_1}] and \code{edges_y}[1:\eqn{k_2}] (left-sided)
}
  \item{extendOutput}{
 Optional, boolean, default FALSE, if TRUE, list is the output, otherwise numerical matrix, see below.
}
}
\details{
edges_x and edges_y must be sorted non-decreasingly. Beware that kernels are centers of bins, \code{edges_x}, \code{edges_y} are borders of bins.
If edges are given, \code{edges_x}, \code{edges_y} can contain Inf,-Inf borders. In that case,
edges always define n-1 bins lying within the edges. data outside first edge or last edge are ignored.
Edges have either to be sorted non-decreasingly or set sort=T.


If edges are not given, set sort=T. In this case, they define the unique number of points. 
Then the number of edges internally sets the number of bins.

Beware that in matrix notation, \code{count matrix} would be expected be ordered [1:\eqn{k_1},1:\eqn{k_2}] instead of [1:\eqn{k_2},1:\eqn{k_1}]. Here we use the ordering that intuitively is given in plot(x,y), i.e. x are columns and y are rows. 
}
\value{
\code{extendOutput==FALSE}: [1:\eqn{k_2},1:\eqn{k_1}] numerical matrix of counts

If \code{extendOutput==TRUE}, then list of named elements

\code{count matrix}: [1:\eqn{k_2},1:\eqn{k_1}]    numerical matrix of counts

\code{x_idx}:        [1:\eqn{k_1}]         numerical vector of counts for \code{x} based on \code{edges_x}   

\code{y_idx}:        [1:\eqn{k_2}]         numerical vector of counts for \code{y} based on \code{edges_y}       

}

\author{
Michael Thrun
}

\seealso{
\code{\link{table}}, \code{\link{findInterval}}
}
\examples{
if(requireNamespace("FCPS")){
data(Hepta,package ="FCPS")
Cls=Hepta$Cls
Cls1=Cls+1
#k unqiue points define k bins
fast_table_num(Cls,Cls1,

redefine = FALSE,names=TRUE)==as.matrix(table(Cls,Cls1))
}
#k unqiue points define k bins
tab=fast_table_num(rnorm(100),rnorm(100),redefine=FALSE,sort=TRUE)

#set k+1 edges to get k bins
x=rnorm(100)
y=rnorm(100)
binsxy=5
edgex=seq(from=min(x),to=max(x),length.out=binsxy+1)
edgesy=seq(from=min(y),to=max(y),length.out=binsxy+1)
fast_table_num(x,y,edgex,edgesy,

redefine=FALSE,names=TRUE,rightmost.closed =TRUE)

#definition of counts analog to plotting
x = c(rnorm(1000, mean=-5), rnorm(1000, mean=5))
y = rnorm(2000)
edgesx = seq(min(x), max(x), length.out=512+1)
edgesy = seq(min(y), max(y), length.out=256+1)
joint_table = fast_table_num(x, y, edgesx, edgesy)
\donttest{
plot(x,y)
plot(colSums(joint_table),xlab="x marginal",

ylab="sum of counts",main="x-values are stored in columns")

plot(rowSums(joint_table),xlab="y marginal",

ylab="sum of counts",main="y-values are stored in rows")
}
}
\keyword{utilities}
