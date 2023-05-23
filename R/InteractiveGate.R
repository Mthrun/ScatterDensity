InteractiveGate=function(Data,GateVars,Gatename="Leukocytes",PlotIt=FALSE){
  
  if (!requireNamespace('flowCore',quietly = TRUE)) {
    warning(
      'Subordinate "flowCore" package is missing. No computations are performed.
            Please install the package which is defined in "Suggests".'
    )
  }
  if (!requireNamespace('flowWorkspace',quietly = TRUE)) {
    warning(
      'Subordinate "flowWorkspace" package is missing. No computations are performed.
            Please install the package which is defined in "Suggests".'
    )
  }
  if (!requireNamespace('flowGate',quietly = TRUE)) {
    warning(
      'Subordinate "flowGate" package is missing. No computations are performed.
            Please install the package which is defined in "Suggests".'
    )
  }
  if (!is.matrix(Data)){
    warning('InteractiveGate:: Data has to be a matrix. Calling as.matrix()')
    Data=as.matrix(Data)
  }
  if(mode(Data)!="numeric"){
    warning('InteractiveGate:: Data has to be numeric.Transforming data to numeric')
    mode(Data)="numeric"
  }
  if(length(GateVars)<2){
    stop("InteractiveGate: Within GateVars two variables were not provided")
  }
  if(length(GateVars)>2){
    warning("InteractiveGate: Within GateVars more than two variables were  provided. Using the first two.")
  }
  if(is.null(colnames(Data))){
    stop("InteractiveGate: colnames of Data has to be defined otherwise names of variabes are unknown.")
  }
  if(length(intersect(colnames(Data),GateVars))<2){
    stop("InteractiveGate: GateVars were not found in Data. Names of selected variabes are unknown.")
  }
  frame=flowCore::flowFrame(Data)
  fs=flowCore::flowSet(frame)
  gs <- flowWorkspace::GatingSet(fs)
  
  mod=flowGate::gs_gate_interactive(gs,
                          filterId = Gatename,#gatename for reference
                          dims = list(GateVars[1], GateVars[2]))
  Gate=mod$Gate
  Polygon=Gate@boundaries
  GatedV=PolygonGate(Data,Polygon = Polygon,GateVars)
  
  DataInGate=GatedV$DataInGate
  if(isTRUE(PlotIt)){
    DensityScatter.DDCAL(X = Data[,GateVars[1]],Y = Data[,GateVars[2]],
                         xlab = GateVars[1],ylab = GateVars[2],SDHorPDE = FALSE,
                         Silent = FALSE,Marginals = TRUE,Plotter = "native",Polygon = Polygon)

  }
  Index=GatedV$DataInGate
  
  V=list(Polygon=Polygon,Index=Index,DataInGate=DataInGate)
  names(V)[3]=Gatename
  
  return(V)
}