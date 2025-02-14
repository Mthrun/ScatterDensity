SampleScatter=function(X,Y,ThresholdPoints=20,DensityThreshold,nbins=100,PlotIt=FALSE){

  
V=SmoothedDensitiesXY(X,Y = Y,nbins = nbins,Compute = "Cpp")

range(V$GridDensity)

dd=as.vector(V$GridDensity)
if(missing(DensityThreshold)){
  
  if(!requireNamespace("ABCanalysis")){
    warning("SampleScatter: Please install ABCanalysis in order to use SampleScatter().")
    return("SampleScatter: Please install ABCanalysis in order to use SampleScatter().")
  }
  
  
V2=ABCanalysis::ABCanalysis(dd,PlotIt = F)
nprop=length(V2$Cind)/length(as.vector(dd))
ddnew=dd
ddnew_A=dd
while (nprop>0.8) {
  V2=ABCanalysis::ABCanalysis(ddnew,PlotIt = F)
  nprop=length(V2$Cind)/length(as.vector(dd))
  ddnew_A=ddnew[V2$Aind]
  ddnew=ddnew[V2$Cind]
}
max(ddnew,na.rm = T)
max(ddnew_A,na.rm = T)


ind=which(V$GridDensity>max(ddnew_A,na.rm = T),arr.ind = T)
}else{
  ind=which(V$GridDensity>DensityThreshold,arr.ind = T)
}
V$Xkernels[ind[1,2]] & V$Ykernels[ind[1,1]]

x_min = V$Xkernels[ind[,2]-1]
x_max = V$Xkernels[ind[,2]]
y_min = V$Ykernels[ind[,1]-1]
y_max = V$Ykernels[ind[,1]]

ss=c()
new_ind=c()
ind_all=c()
for(i in 1:length(x_min)){
  ind_cur= X >= x_min[i] & X <= x_max[i] &  Y >= y_min[i] & Y <= y_max[i]
  ss[i]=sum(ind_cur)
  if(ss[i]>ThresholdPoints){
    ind_all_cur=which(ind_cur)
    new_ind=c(new_ind,sample(ind_all_cur,ThresholdPoints))
    ind_all=c(ind_all_cur,ind_all)
  }
}
sum(ss)/dim(Data_pb_sample)[1]
length(new_ind)
length(unique(new_ind))
length(ind_all)
length(unique(ind_all))

ind_sub=sort(c(new_ind,setdiff(1:nrow(Data_pb_sample),ind_all)),decreasing = F)
length(ind_sub)
length(unique(ind_sub))
if(isTRUE(PlotIt)){
  xlab = deparse1(substitute(X))
  ylab = deparse1(substitute(Y))
  plot(X[ind_sub],Y[ind_sub],type="p",pch=20,xlab=xlab,ylab=ylab)
}
SubsampleInd=ind_sub
return(SubsampleInd)
}