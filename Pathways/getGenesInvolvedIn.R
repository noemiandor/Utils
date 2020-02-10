getGenesInvolvedIn<-function(pathway, dbs="ReactomePA"){
  goi=c()
  if("ReactomePA" %in% dbs){
    library(ReactomePA)
    library(igraph)

    goi <- ReactomePA::viewPathway(pathway)$data$name
    #p <- pathways("hsapiens", 'reactome')[[pathway]]
    #g <- try(pathwayGraph(p),silent=T)
    #if(class(g)!="try-error"){
    #  gg <- igraph.from.graphNEL(g)
    #  gg <- as.undirected(gg)
    #  V(gg)$name <- sub("[^:]+:", "", V(gg)$name)
    #  goi=c(goi,V(gg)$name)
    #}
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
