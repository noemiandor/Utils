colorCellsByPathway<-function(seu, pathways, cells=NULL, pt.size=0.5, outF=NULL){
  seu@data=as.matrix(seu@data)
  
  if(is.null(cells)){
    cells=colnames(seu@data)
  }
  
  ##Matrix of pathway features
  m=matrix(NA,length(pathways),length(cells));
  colnames(m)=cells; rownames(m)=pathways;
  
  gpp=rep(NA,length(pathways)); names(gpp)=pathways; ##Genes per pathway
  for(p in pathways){
    ##Get genes of interest, involved in pathway
    goi=getGenesInvolvedIn(p)
    goi=intersect(goi,rownames(seu@data))
    
    gpp[p]=length(goi)
    if(isempty(goi)){
      print(paste(p,"-genes not in",seu@project.name))
    }
    ##Calculate cummulative expression of goi per cell
    ib=match(goi,rownames(seu@data))
    e=apply(seu@data[ib,cells,drop=F],2,sum)
    m[p,names(e)]=e
  }
  
  ii=which(apply(!is.na(m) & m>0,1,any))
  if(isempty(ii)){
    print(paste("None of the pathways recorded in",seu@project.name))
  }
  m=m[ii,,drop=F]
  rownames(m)=paste0(seu@project.name," > ",rownames(m),": ",gpp[rownames(m)]) ##Add # genes/pathway to row names
  
  
  ##Color every cell by expression
  seu_ <- CreateSeuratObject(raw.data = m, min.cells = 0, min.genes = -1, project = seu@project.name)
  seu_@dr$tsne=seu@dr$tsne;   seu_@dr$tsne@cell.embeddings=seu_@dr$tsne@cell.embeddings[cells,]
  try(graphics.off())
  if(!is.null(outF)){
    file.remove(outF); pdf(file = outF)
  }
  for(p in rownames(seu_@data)){
    FeaturePlot(seu_,features.plot=p,no.legend=F, max.cutoff="q95",pt.size=pt.size)
  }
  if(!is.null(outF)){
    dev.off()
  }
}