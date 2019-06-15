annotateWithGeneLoc<-function(sc){
  geneQ=read.table('/mnt/ix2/nandor/Projects/GenHetAcrossCancers/data/Annotation/HG19_GeneCoordinates_UCSC_Autosomes.txt',sep="\t",header=T,check.names = F,stringsAsFactors = F,as.is = T,fill=T,strip.white = T);
  ##Prep gene coor
  geneQ=geneQ[!duplicated(geneQ$GeneSymbol),]
  rownames(geneQ)=geneQ$GeneSymbol; geneQ=geneQ[,c("chr","startpos","endpos")]
  geneQ=geneQ[!is.na(geneQ[,"startpos"]),]
  geneQ=geneQ[!is.na(as.numeric(geneQ$chr)),]
  geneQ$chr=as.numeric(geneQ$chr)
  
  ##keep measured genes only
  goi=intersect(sc$genes,rownames(geneQ))
  geneQ=geneQ[goi,]
  
  ##sort genes
  tmp=sort(geneQ$startpos,index.return=T)
  geneQ=geneQ[tmp$ix,]
  tmp=sort(geneQ$chr,index.return=T)
  geneQ=geneQ[tmp$ix,]  
  
  ##map genes
  ia=match(rownames(geneQ),sc$genes)
  sc$genes=sc$genes[ia]
  sc$mat=sc$mat[,ia]
  
  sc$geneLoc=geneQ;
  return(sc)
}