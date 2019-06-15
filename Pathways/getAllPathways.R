getAllPathways<-function(include_genes=F, loadPresaved=F){
  if(loadPresaved){
    if(isempty( grep('linux', R.Version()$os))){
      load("~/Projects/GenHetAcrossCancers/data/Pathways/REACTOME/Pathway_GeneSets_Reactome.RObj")
    }else{
      load("/mnt/ix2/nandor/Projects/GenHetAcrossCancers/data/Pathways/REACTOME/Pathway_GeneSets_Reactome.RObj")
    }
    pathways <- gs
  }else{
    library(ReactomePA)
    library(matlab)
    xx <- as.list(reactome.db::reactomePATHID2NAME)
    tmp=sapply(xx,function(x) strsplit(x,": ")[[1]])
    tmp=tmp[sapply(tmp,"[[",1)=="Homo sapiens"]; ##Keep only pathways in human
    pathways=sapply(tmp,"[[",2)  
    pathways=pathways[!duplicated(pathways)]
    
    if(include_genes){
      names=pathways[T]
      pathways=sapply(pathways, getGenesInvolvedIn)
      names(pathways)=names
    }
  }
  return(pathways)
}