getGenesInvolvedIn<-function(pathways, dbs="ReactomePA"){
  if("ReactomePA" %in% dbs){
    allGoi=list()
    library(ReactomePA)
    library(igraph)
    allP =  try(convertIdentifiers(pathways("hsapiens", 'reactome'), "symbol"),silent=T)
    for(pathway in pathways){
      goi=c()
      p = allP[[pathway]]
      g <- try(pathwayGraph(p),silent=T)
      if(class(g)!="try-error"){
        gg <- igraph.from.graphNEL(g)
        gg <- as.undirected(gg)
        V(gg)$name <- sub("[^:]+:", "", V(gg)$name)
        goi=c(goi,V(gg)$name)
      }
      allGoi[[pathway]] = goi;
    }
  }
  if("INPATH" %in% dbs){
    allGoi=c()
    mydb = dbConnect(MySQL(), user='noemi', password='lala', dbname='INPATH')
    genes<-dbGetQuery(mydb,paste('select gene_names from pathway_genes where pathway_name=\'',pathway,'\'',sep=""))$gene_names
    .disconnectAll()
    allGoi=c(allGoi,genes)
    allGoi=gsub("^SYMBOL:","",allGoi, fixed = F)
  }
  return(allGoi)
}
