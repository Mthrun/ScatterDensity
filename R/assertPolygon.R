assertPolygon=function(Polygon,Eps=1e-03){
# Polygon       Polygon representation to normalize and close. Accepted inputs are:
#               - a numerical matrix/data.frame with X and Y columns;
#               - a matrix/data.frame with optional VertexNo plus named X and Y columns;
#                     or list containing Vertices or X and Y.
#                     VertexNo, ClusterNo, Ratio, and GeometryType are optional; if VertexNo is supplied,
#                     it is used to restore the path order. An open polygon is closed internally.
#               ClusterNo and Ratio are ignored. GeometryType is optional, but must be "POLYGON"
#               when supplied. If VertexNo is supplied, it is used to restore the path order.
# Eps           Non-negative tolerance used to decide whether first and last vertex already close
#               the polygon. Default: 1e-03. Use Eps=0 for exact closure.
#
# OUTPUT
# Polygon       Numerical matrix with the normalized, ordered and closed X/Y polygon path.
#
# DETAILS
# X and Y are mandatory. At least three unique finite coordinate pairs are required. For backwards
# compatibility, a path whose endpoints differ by more than Eps is closed by appending the first
# point shifted by c(Eps,-Eps). Therefore, Eps=0 appends the first point exactly.
# Michael Thrun

  if(length(Eps)!=1L||!is.numeric(Eps)||!is.finite(Eps)||Eps<0) stop("Eps must be one finite non-negative numeric value.")
  GeometryType=NULL

  # Extract vertices and optional geometry metadata from Cluster2Polygon() output or a compact X/Y list.
  if(is.list(Polygon)&&!is.data.frame(Polygon)){
    if(!is.null(Polygon$GeometryType)) GeometryType=as.character(Polygon$GeometryType)
    if(!is.null(Polygon$Vertices)){
      Polygon=Polygon$Vertices
    }else if(!is.null(Polygon$X)&&!is.null(Polygon$Y)){
      if(length(Polygon$X)!=length(Polygon$Y)) stop("Polygon$X and Polygon$Y must have the same length.")
      VertexNo=Polygon$VertexNo
      Polygon=data.frame(X=Polygon$X,Y=Polygon$Y)
      if(!is.null(VertexNo)){
        if(length(VertexNo)!=nrow(Polygon)) stop("Polygon$VertexNo must have the same length as Polygon$X and Polygon$Y.")
        Polygon=cbind(VertexNo=VertexNo,Polygon)
      }
    }else{
      stop("Polygon must contain Vertices or X and Y.")
    }
  }

  # Normalize matrix/data.frame input before using optional metadata and coordinate column names.
  if(!is.matrix(Polygon)&&!is.data.frame(Polygon)) stop("Polygon has to be a matrix, data.frame, or Cluster2Polygon output.")
  Polygon=as.data.frame(Polygon)

  # Validate optional GeometryType metadata; POINT and LINESTRING do not define an interior.
  if(is.null(GeometryType)&&"GeometryType"%in%colnames(Polygon)) GeometryType=unique(as.character(Polygon$GeometryType))
  if(!is.null(GeometryType)&&(length(GeometryType)!=1L||is.na(GeometryType)||toupper(GeometryType)!="POLYGON")) stop("assertPolygon requires GeometryType='POLYGON'.")

  # Restore the polygon path order when optional VertexNo information is available.
  if("VertexNo"%in%colnames(Polygon)){
    if(!is.numeric(Polygon$VertexNo)||any(!is.finite(Polygon$VertexNo))) stop("Polygon$VertexNo must contain finite numeric values.")
    if(anyDuplicated(Polygon$VertexNo)) stop("Polygon$VertexNo must not contain duplicate values.")
    Polygon=Polygon[order(Polygon$VertexNo),,drop=FALSE]
  }

  # Prefer named X/Y columns; retain backwards compatibility with an unnamed two-column polygon.
  if(all(c("X","Y")%in%colnames(Polygon))){
    Polygon=Polygon[,c("X","Y"),drop=FALSE]
  }else if(ncol(Polygon)==2L){
    Polygon=Polygon[,1:2,drop=FALSE]
  }else{
    stop("Polygon must contain columns X and Y or exactly two coordinate columns.")
  }

  # Convert the validated coordinates to the numerical matrix required by point-in-polygon routines.
  if(!is.numeric(Polygon[[1L]])||!is.numeric(Polygon[[2L]])) stop("Polygon X and Y coordinates must be numeric.")
  Polygon=as.matrix(Polygon)
  colnames(Polygon)=c("X","Y")
  if(any(!is.finite(Polygon))) stop("Polygon X and Y coordinates must be finite.")
  if(nrow(unique(Polygon))<3L) stop("Polygon must contain at least three unique X/Y points.")

  # Close the path only when first and last vertex differ by more than Eps; Eps=0 closes it exactly.
  Diff=abs(Polygon[1L,]-Polygon[nrow(Polygon),])
  if(any(Diff>Eps)) Polygon=rbind(Polygon,Polygon[1L,]+c(Eps,-Eps))
  return(Polygon)
}
