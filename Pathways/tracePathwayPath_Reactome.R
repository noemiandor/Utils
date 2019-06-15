tracePathwayPath_Reactome<-function(target, pr=NULL){
  if(is.null(pr)){
    pr=read.table("~/Downloads/ReactomePathwaysRelation.txt",sep="\t")
  }
  path=c()
  while(!is.na(target) && !isempty(target) && !target %in% path){
    y=target
    path=c(y, path)
    target=as.character(pr$V1[pr$V2==target])[1]; ##@TODO: why first hit only?
  }
  return(path)
}