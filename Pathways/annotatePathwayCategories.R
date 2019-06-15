annotatePathwayCategories<-function(targets){
  pr=read.table("~/Projects/GenHetAcrossCancers/data/Pathways/REACTOME/ReactomePathwaysRelation.txt",sep="\t", stringsAsFactors = F)
  pnames=read.table("~/Projects/GenHetAcrossCancers/data/Pathways/REACTOME/NCBI2Reactome_All_Levels_HSapiens.txt",sep="\t", comment.char = "", stringsAsFactors = F, check.names = F)
  colnames(pnames)=c("SourceDBidentifier", "ReactomeStableidentifier","URL","Event Name","Evidence Code","Species")
  pnames=pnames[pnames$Species=="Homo sapiens",]
  pnames=pnames[!duplicated(pnames$ReactomeStableidentifier),]
  rownames(pnames)=pnames$ReactomeStableidentifier
  pr=pr[pr$V1 %in% pnames$ReactomeStableidentifier & pr$V2 %in% pnames$ReactomeStableidentifier,]
  
  l1=pathwayName2ReactomeID(targets=targets, pnames)
  l0=cbind(l1,repmat(NA,length(l1),2)); colnames(l0)[2:3]=c("l0","l0_name")
  ii=which(!is.na(l0[,"l1"]));      print(paste("Tracing",length(ii),"targets..."))
  l0[ii,"l0"]= sapply(l0[ii,"l1"], function(x) tracePathwayPath_Reactome(x, pr=pr)[1])
  l0[ii,"l0_name"]= sapply(l0[ii,"l0"], function(x) pnames[x,"Event Name"])
  return(l0)
}