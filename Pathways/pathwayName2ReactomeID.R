pathwayName2ReactomeID<-function(targets, pnames, verbose=T){
  library("stringdist"); 
  ids=rep(NA, length(targets));   names(ids)=tolower(targets)
  
  x=intersect_MatlabV(tolower(targets), tolower(pnames$`Event Name`)); 
  ids[x$ia]=pnames$ReactomeStableidentifier[x$ib] 
  
  ##Not all targets have perfect matches
  if(length(x$a)<length(targets)){
    x=setdiff(tolower(targets), tolower(pnames$`Event Name`)); 
    print(paste("Perfect match not found for all targets. Looking for best approximate match for",length(x),"targets..."))
    tmp=sapply(x, stringdist, tolower(pnames$`Event Name`))
    d=apply(tmp,2,min);    
    
    ids[tolower(x)]=pnames$`Event Name`[apply(tmp,2,which.min)]
    print(paste("Found approximate matches for",sum(d<=10),"targets:"))
    if(verbose){
      print(ids[tolower(x)[d<=10]])
    }
    print(paste("Found no good matches for",sum(d>10),"targets:"))
    if(verbose){
      print(ids[tolower(x)[d>10]])
    }
    
    ids[tolower(x)]=pnames$ReactomeStableidentifier[apply(tmp,2,which.min)]
    ids[tolower(x)[d>10]]=NA
    
  }
  return(ids)
}