Mask2Polygon=function(Data,Mask,ClusterNo,GateVars,PlotIt=FALSE){
  # Polygon=Mask2Polygon(Data,Mask,ClusterNo,GateVars,PlotIt=FALSE)
  # Constructs a convex polygon for observations selected by a logical mask in two gate variables.
  #
  # INPUT
  # Data        Numerical matrix or data.frame with observations in rows and variables in columns.
  # Mask        Logical vector of length nrow(Data). TRUE selects an observation; NA is treated as FALSE.
  # ClusterNo   Scalar identifier of the population or hierarchy node, stored in Polygon$ClusterNo.
  # GateVars    Two column indices or two column names defining the X and Y coordinates. If more than
  #             two variables are supplied, only the first two are used with a warning.
  # PlotIt      If TRUE, Cluster2Polygon() plots the selected observations and convex polygon.
  #
  # OUTPUT
  # Polygon     Standardized four-element list Vertices, ClusterNo,
  #             Ratio and GeometryType. Ratio is always 1.
  #
  # DETAILS
  # The hierarchy is defined by Mask, e.g. ChildMask=ParentMask&LocalGateRule. GateVars only define
  # the local two-dimensional representation. The result can be passed directly to PointsInPolygon()
  # or PolygonGate(). Applying the convex hull again can include additional observations that were not
  # selected by Mask but lie inside the convex hull.
  # Michael Thrun
  
  if(is.null(dim(Data))) stop("Data must be a matrix or data.frame.")
  if(!is.logical(Mask)||length(Mask)!=nrow(Data)) stop("Mask must be a logical vector with length nrow(Data).")
  Mask[is.na(Mask)]=FALSE
  if(!any(Mask)) stop("Mask does not select any observations.")
  if(length(GateVars)<2L) stop("Mask2Polygon: Within GateVars two variables were not provided.")
  if(length(GateVars)>2L){
    warning("Mask2Polygon: Within GateVars more than two variables were provided. Using the first two.")
    GateVars=GateVars[1:2]
  }
  if(length(unique(GateVars))!=2L) stop("GateVars must identify two different variables.")
  GateData=Data[Mask,GateVars,drop=FALSE]
  if(ncol(GateData)!=2L) stop("One or both GateVars were not found in Data.")
  LocalCls=rep(ClusterNo,nrow(GateData))
  Polygon=Cluster2Polygon(Data=GateData,Cls=LocalCls,ClusterNo=ClusterNo,Ratio=1,PlotIt=PlotIt)
  return(Polygon)
}
