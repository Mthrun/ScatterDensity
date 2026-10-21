PolygonGate=function(Data,Polygon,GateVars,PlotIt=FALSE,PlotSampleSize=1000){
  #V=PolygonGate(Data,Polygon,GateVars,PlotIt=FALSE,PlotSampleSize=1000)
  #A specific two-dimensional polygon gate is applied to flow-cytometry data.
  #
  # INPUT
  # Data                numerical matrix n x d
  # Polygon             Polygon a numerical matrix/data.frame with X and Y,
  #                     optionally VertexNo, or list containing Vertices or X and Y.
  #                     VertexNo, ClusterNo, Ratio, and GeometryType are optional; if VertexNo is supplied,
  #                     it is used to restore the path order. An open polygon is closed internally.
  # GateVars            two column indices or two variable names defining the X and Y coordinates of the gate
  # PlotIt              if TRUE: plots a sample of data in the two selected variables and marks point inside the gate as yellow and outside as magenta
  # PlotSampleSize      size of the plottet sample
  #
  #OUTPUT
  #   list of
  # DataInGate        m x d numerical matrix with m<=n of data points inside the gate
  # InGateInd         index of length m for the datapoints in original matrix
  # InGateMask        logical vector of length nrow(Data), suitable for hierarchical mask combinations
  
  #DETAILS
  # Gates are always two dimensional although all dimensions of points inside the gate are returned.
  # Polygon input is normalized by PointsInPolygon() through assertPolygon(); X and Y are mandatory,
  # whereas VertexNo and the Cluster2Polygon metadata are optional. A hierarchical child mask can be
  # constructed by ChildMask=ParentMask&PolygonGate(...)$InGateMask.
  
  #   Michael Thrun
  #EXAMPLE
  if(length(GateVars)<2){
    stop("PolygonGate: Within GateVars two variables were not provided")
  }
  if(length(GateVars)>2){
    warning("PolygonGate: Within GateVars more than two variables were provided. Using the first two.")
    # Keep only the first two variables because a polygon gate is two-dimensional.
    GateVars=GateVars[1:2]
  }
  GateVar1=GateVars[1]
  GateVar2=GateVars[2]
  if(!is.numeric(GateVar1) & !is.numeric(GateVar2)){
    if(is.null(colnames(Data))){
      text="One or two of the Gate variables not found in data, because strings are used for 'GateVar1'/'GateVar1' but 'Data' has no column names. Please provide column indices instead or name columns of data."
      message(text)
      return(text)
    }
  }
  
  X=Data[,GateVar1]
  Y=Data[,GateVar2]
  
  GateData=cbind(X,Y)
  
  if(ncol(GateData)==2){
    # Reuse PointsInPolygon()/assertPolygon() so all standardized polygon representations are accepted.
    InGateCls=ScatterDensity::PointsInPolygon(GateData,Polygon=Polygon,PlotIt=FALSE)
    if(is.null(InGateCls)) return(NULL)
    # Preserve the logical result so gates can be combined hierarchically with parent masks.
    InGateMask=InGateCls==2L
    ind=which(InGateMask)
    DataInGate=Data[ind,,drop=FALSE]
    if(isTRUE(PlotIt)){
      GateDataSample=GateData[sample(1:nrow(GateData),replace = T,size = PlotSampleSize),]
      if(is.numeric(GateVar1))
        xlab=paste("Index",GateVar1,"specifying Variable", colnames(Data)[GateVar1])
      else
        xlab=GateVar1
      
      if(is.numeric(GateVar2))
        ylab=paste("Index",GateVar2,"specifying Variable", colnames(Data)[GateVar2])
      else
        ylab=GateVar2
      
      XYlab=colnames(GateDataSample)
      ScatterDensity::PointsInPolygon(GateDataSample,Polygon=Polygon,PlotIt=TRUE,main="Yellow = In Gate",xlab=xlab,ylab=ylab)
    }
    return(list(DataInGate=DataInGate,InGateInd=ind,InGateMask=InGateMask))
  }else{
    text="One or two of the Gate variables not found in data. Please check manually"
    message(text)
    return(text)
  }
}
