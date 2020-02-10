getGenesInvolvedIn<-function(pathway, dbs="ReactomePA"){
  goi=c()
  if("ReactomePA" %in% dbs){
    library(ReactomePA)
    library(igraph)

    g <- try(ReactomePA::viewPathway(pathway))
    #p <- pathways("hsapiens", 'reactome')[[pathway]]
    #g <- try(pathwayGraph(p),silent=T)
    if(class(g)!="try-error"){
      goi = g$data$name
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
