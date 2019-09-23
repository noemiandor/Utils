getGenesInvolvedIn<-function(pathway, dbs="ReactomePA"){
  goi=c()
  if("ReactomePA" %in% dbs){
    library(ReactomePA)
    library(igraph)
    jpeg(filename = paste0("~/Downloads/",filesep,"tmp.jpg"), width=10.5, height=6.5, units="in", res=200)
    tmp=try(viewPathway(pathway),silent=T)
    dev.off()
    if(class(tmp)!="try-error"){
      goi=c(goi,V(tmp)$name)
    }
  }
  if("INPATH" %in% dbs){
    mydb = dbConnect(MySQL(), user='noemi', password='lala', dbname='INPATH')
    genes<-dbGetQuery(mydb,paste('select gene_names from pathway_genes where pathway_name=\'',pathway,'\'',sep=""))$gene_names
    .disconnectAll()
    goi=c(goi,genes)
  }
  goi=gsub("^SYMBOL:","",goi, fixed = F)
  return(goi)
}
