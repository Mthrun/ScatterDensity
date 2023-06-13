assertPolygon=function(Polygon,Eps=1e-03){
  #make sure the polygon is closed except of epsilon
  Diff=abs(head(Polygon,1)-tail(Polygon,1))
  if(any(Diff>Eps))
    Polygon=rbind(Polygon,head(Polygon,1)+c(Eps,-Eps))
  
  return(Polygon)
}