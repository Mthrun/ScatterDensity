PolygonGate=function(Data,Polygon,GateVars,PlotIt=FALSE,PlotSampleSize=1000){
  #V= PolygonGate(Data, GateVar1, GateVar2, Polygon, PlotIt = FALSE, PlotSampleSize = 1000)
  #A specific Gate defined by xy coordinates that result in a closed polygon is applied to the flowcytometry data.
  #
  # INPUT
  # Data                numerical matrix n x d
  # GateVar1            either column index in Data of X coordinate of gate or variable name as string if matrix has named cols.
  # GateVar2            either column index in Data of Y coordinate of gate or variable name as string if matrix has named cols.
  # Polygon             numerical marix of two columns defining the coordiantes of the polygon, polygon has to be closed, i.e., first and last coordinate has to be the same
  # PlotIt              if TRUE: plots a sample of data in the two selected variables and marks point inside the gate as yellow and outside as magenta
  # PlotSampleSize      size of the plottet sample
  #
  #OUTPUT
  #   list of
  # DataInGate        m x d numerical matrix with m<=n of data points inside the gate
  # InGateInd         index of length m for the datapoints in original matrix
  
  #DETAILS
  # Gates are alwaxs two dimensional, i.e., require two filters, although all dimensions of data are filted by the gates. Only high-dimensional points inside the polygon (gate) are given back
  #   if GateVar1 or GateVar2 is not found a text is given back which will state this issue
  
  #   Michael Thrun
  #EXAMPLE
  if(length(GateVars)<2){
    stop("PolygonGate: Within GateVars two variables were not provided")
  }
  if(length(GateVars)>2){
    warning("PolygonGate: Within GateVars more than two variables were  provided. Using the first two.")
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
  
  if(dim(GateData)[2]){
    InGateCls=ScatterDensity::PointsInPolygon(GateData,Polygon = Polygon,PlotIt = F)
    ind=which(InGateCls==2)
    DataInGate=Data[ind,]
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
      ScatterDensity::PointsInPolygon(GateDataSample,Polygon = Polygon,
                                      PlotIt = T,main = "Yellow = In Gate",xlab=xlab,ylab=ylab)
    }
    return(list(DataInGate=DataInGate,InGateInd=ind))
  }else{
    text="One or two of the Gate variables not found in data. Please check manually"
    message(text)
    return(text)
  }
  
  
}
