getAllPathways<-function(include_genes=F, loadPresaved=F){
  library(matlab)
  if(loadPresaved){
    if(isempty( grep('linux', R.Version()$os))){
      load("~/Projects/code/RCode/github/Utils/Pathways/Pathway_GeneSets_Reactome.RObj")
    }else{
      load("/mnt/ix1/Resources/data/Pathways/REACTOME/Pathway_GeneSets_Reactome.RObj")
    }
    pathways <- gs
  }else{
    library(ReactomePA)
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
