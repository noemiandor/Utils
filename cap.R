cap<-function(X,q=NULL,minmax=NULL){
  if(is.null(minmax)){
    minmax=quantile(as.numeric(X),q,na.rm=T); 
  }
  X[!is.na(X) & X>minmax[2]]=minmax[2]
  X[!is.na(X) & X<minmax[1]]=minmax[1]
  return(X)
}
