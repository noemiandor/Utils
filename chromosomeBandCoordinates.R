chromosomeBandCoordinates<-function(host="dec2017.archive.ensembl.org"){
  devtools::source_url("https://github.com/noemiandor/Utils/blob/master/grpstats.R?raw=TRUE")
  ensembl=try(biomaRt::useMart("ENSEMBL_MART_ENSEMBL",dataset="hsapiens_gene_ensembl" ,host=host))
  mart = biomaRt::useDataset("hsapiens_gene_ensembl",mart=ensembl)
  segments<-biomaRt::getBM( c("band" ,"chromosome_name","start_position","end_position"), mart=mart)
  segments = segments[segments[,"chromosome_name"] %in% 1:22,]
  segments$band[grep("^p",segments$band)]="p"
  segments$band[grep("^q",segments$band)]="q"
  segments[,"start_position"]=as.numeric(segments[,"start_position"])
  segments[,"end_position"]=as.numeric(segments[,"end_position"])
  segments[,"chromosome_name"]=as.numeric(segments[,"chromosome_name"])
  segments$band = paste0(segments$chromosome_name,segments$band)
  segments = grpstats(segments[,c("chromosome_name","start_position","end_position")], segments$band,c("min","max"))
  segments = cbind(segments$min, segments$max[rownames(segments$min),])
  segments=segments[sort(segments[,"start_position"], index.return=T)$ix,]
  segments=segments[sort(segments[,"chromosome_name"], index.return=T)$ix,]
  segments= segments[,c(1,2,6)]
  segments=as.data.frame(segments)
  
  segments$segmentLength=1+segments$end_position-segments$start_position    
  segments$segmentLength_Mb=segments$segmentLength/1E6
  colnames(segments)[1:3]=c("chr","startpos","endpos")
  return(segments)
}
