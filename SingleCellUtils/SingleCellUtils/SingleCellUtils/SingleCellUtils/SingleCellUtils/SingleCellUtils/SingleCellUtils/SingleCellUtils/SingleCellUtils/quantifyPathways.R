quantifyPathways<-function(X, pathways=NULL, min.genes=5){
  library("reactome.db")
  
  if(is.null(pathways)){
    pathways=getAllPathways()
  }
  
  M=matrix(NA,length(pathways),ncol(X)); colnames(M)=colnames(X); rownames(M)=pathways
  for(p in pathways){ 
    genes<-getGenesInvolvedIn(p)
    N=length(genes)
    if(N>=min.genes){
      genes=intersect(genes,rownames(X))
      M[p,]=apply(X[genes,,drop=F],2,sum)/N
    }
  }
  
  M=M[apply(!is.na(M),1,any),]
  
  return(M)
}