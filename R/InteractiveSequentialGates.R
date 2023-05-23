InteractiveSequentialGates=function(Data,GateVarsMat,Gatenames,PlotIt=FALSE){
  
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
    warning('InteractiveSequentialGates:: Data has to be a matrix. Calling as.matrix()')
    Data=as.matrix(Data)
  }
  if(mode(Data)!="numeric"){
    warning('InteractiveSequentialGates:: Data has to be numeric.Transforming data to numeric')
    mode(Data)="numeric"
  }
  if(dim(GateVarsMat)[2]<2){
    stop("InteractiveSequentialGates: Within GateVarsMat two variables were not provided")
  }
  if(dim(GateVarsMat)[2]>2){
    warning("InteractiveSequentialGates: Within GateVarsMat more than two variables were  provided. Using the first two.")
  }
  if(is.null(colnames(Data))){
    stop("InteractiveSequentialGates: colnames of Data has to be defined otherwise names of variabes are unknown.")
  }
  if(length(intersect(colnames(Data),GateVarsMat))<2){
    stop("InteractiveSequentialGates: GateVarsMat were not found in Data. Names of selected variabes are unknown.")
  }
  if(missing(Gatenames)){
    Gatenames=paste0("Gate",1:nrow(GateVarsMat))
  }
  if(dim(GateVarsMat)[1]!=length(Gatenames)){
    warning("InteractiveSequentialGates: Number of rows of GateVarsMat is unequal to number of Gatenames. Generating arbitrary gatenames.")
    Gatenames=paste0("Gate",1:nrow(GateVarsMat))
  }
  
  frame=flowCore::flowFrame(Data)
  fs=flowCore::flowSet(frame)
  gs <- flowWorkspace::GatingSet(fs)
  
  GateVarsMat=rbind(c("SS", "FS"),c("CD45", "CD19"))
  Gatenames=c("Lymphocytes","Leukocytes")
  FilterIDs=c("root",Gatenames[1:(length(Gatenames)-1)])
  
  FullMat=cbind(Gatenames,GateVarsMat,FilterIDs)
  
  strategy <- tibble::tibble(
    filterId = FullMat[1, 1],        # First column as is
    dims = list(FullMat[1, 2:3]),  # Combine columns 2 and 3 into a list
    subset = FullMat[1, 4]
  )
  for(i in 2:nrow(GateVarsMat)){
    strategy=dplyr::bind_rows(strategy,tibble::tibble(
      filterId = FullMat[i, 1],        # First column as is
      dims = list(FullMat[i, 2:3]),  # Combine columns 2 and 3 into a list
      subset = FullMat[i, 4]
    ))
  }
  V=flowGate::gs_apply_gating_strategy(gs, gating_strategy = strategy)
  
  Gates=lapply(V,"[[","Gate")
  PolygonList=c()
  for (i in seq_along(Gates)) {
    Polygon_cur <- Gates[[i]]@boundaries
    PolygonList[[i]]=Polygon_cur
  }
  names(PolygonList)=Gatenames
  
  
  InGateInd=list()
  DataGatedV=list()
  for(i in 1:nrow(GateVarsMat)){
    if(i==1)
      GatedV=PolygonGate(Data,Polygon =  PolygonList[[i]],GateVars = GateVarsMat[i,])
    else
      GatedV=PolygonGate(DataInGate,Polygon =  PolygonList[[i]],GateVars = GateVarsMat[i,])
  
    DataInGate=GatedV$DataInGate
    DataGatedV[[i]]=DataInGate
    InGateInd[[i]] =GatedV$InGateInd
  }
  names(InGateInd)=Gatenames
  names(DataGatedV)=Gatenames
  
  if(isTRUE(PlotIt)){
    #plot prior to last one
    DataInGate=DataGatedV[[(length(DataGatedV)-1)]]
    DensityScatter.DDCAL(X = DataInGate[,GateVarsMat[i,1]],Y = DataInGate[,GateVarsMat[i,2]],
                         xlab = GateVarsMat[i,1],ylab = GateVarsMat[i,2],SDHorPDE = FALSE,
                         Silent = FALSE,Marginals = TRUE,Plotter = "native",Polygon = PolygonList[[i]])
    
  }
  
  V=list(Polygons=PolygonList,Indices=InGateInd,DataInGates=DataGatedV)
  return(V)
}