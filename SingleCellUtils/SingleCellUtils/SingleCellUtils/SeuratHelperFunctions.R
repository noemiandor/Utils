
write10Xtrio<-function(X,prefix,overwrite=T){
  if(!overwrite && file.exists(paste0(prefix,"_matrix.mtx"))){
    print(paste0(prefix,"_matrix.mtx already exists and overwrite set to FALSE."))
    return()
  }
  if(isempty(grep("::",rownames(X)[1]))){
    genes=rownames(X)
  }else{
    genes=sapply(rownames(X),strsplit,"::")
    genes=t(as.data.frame(genes));
  }
  writeMM(as(X,"sparseMatrix"), file=paste0(prefix,"_matrix.mtx") )
  write.table(genes,file=paste0(prefix,"_genes.tsv"),quote = F,row.names = F,col.names = F, sep="\t" )
  write.table(colnames(X),file=paste0(prefix,"_barcodes.tsv"),quote = F,row.names = F,col.names = F )
}

read10Xtrio<-function(dirname,regex,geneID="ensembl"){
  f=list.files(dirname,pattern=regex,full.names = T,include.dirs = F)
  ##Read
  # print(paste("Reading",f,"..."))
  counts_nm=readMM(grep('.mtx$',f,value=T)); ##(Genes x Cells) 
  counts_nm=as.matrix(counts_nm)
  genes=read.table(grep('.genes.tsv$',f,value=T),sep="\t")
  if(geneID=="hugo"){
    if(ncol(genes)<2){##Hugo symbol not included in _genes.tsv files --> need to read it externally
      genes$V2=mapEnsemble2Hugo(genes$V1)
    }
    rownames(counts_nm)=genes$V2
  }else if(geneID=="ensembl"){
    rownames(counts_nm)=genes$V1
  }else if(geneID=="both"){
    rownames(counts_nm)=paste(genes$V1,genes$V2,sep="::")
  }
  colnames(counts_nm)=read.table(grep('.barcodes.tsv$',f,value=T))$V1
  return(counts_nm)
}

addMitochondialGeneInfo<-function(seu){
  mito.genes <- grep(pattern = "^MT-", x = rownames(x = seu@assays$RNA@data), value = TRUE)
  percent.mito <- Matrix::colSums(seu@assays$RNA@data[mito.genes, ])/Matrix::colSums(seu@assays$RNA@data)
  seu <- AddMetaData(object = seu, metadata = percent.mito, col.name = "percent.mito")
  VlnPlot(object = seu, features = c("nGene", "nFeature_RNA", "percent.mito"), nCol = 3)
  return(seu)
}

addApoptoticGeneInfo<-function(seu){
  apo.genes <- intersect(rownames(seu@assays$RNA@data), getGenesInvolvedIn("Apoptosis"))
  percent.apo <- Matrix::colSums(seu@assays$RNA@data[apo.genes, ])/Matrix::colSums(seu@assays$RNA@data)
  seu <- AddMetaData(object = seu, metadata = percent.apo, col.name = "percent.apo")
  VlnPlot(object = seu, features = c("nGene", "nFeature_RNA", "percent.apo"), nCol = 3)
  return(seu)
}

##Standard analysis with Seurat consisting of Log-normalization & Scaling by UMIs
standard_SeuratCreate<-function(X, projectName, findVarGenes=F, scale=T, min.cells = 3, min.genes = 500){
  library(matlab)
  seu <- CreateSeuratObject(counts = X, min.cells = min.cells, min.features = min.genes, project = projectName)
  seu <- NormalizeData(object = seu, normalization.method = "LogNormalize", scale.factor = 10000)
  # seu <- addMitochondialGeneInfo(seu)
  if(findVarGenes){
    seu <- FindVariableFeatures(object = seu, mean.function = ExpMean, dispersion.function = LogVMR)
  }
  if(scale){
    seu <- ScaleData(object = seu); #, vars.to.regress = scale
  }
  return(seu)
}

##Seurat dimensionality reduction:
pca_tsne_Seurat<-function(seu_,goi=NULL, do.fast = TRUE, clustering=F, savePCvarPlot=F, resolution = 0.4, MINPCT_VAREXPLAINED=90, dim.embed=2, pcs.compute = 40, perplexity=30){
  if(is.null(seu_@assays$RNA@scale.data)){
    print("Scaling data first before performing PCA...")
    seu_ <- ScaleData(object = seu_);
  }
  if(is.null(goi)){
    goi = seu_@assays$RNA@var.features
  }
  seu_ <- RunPCA(object = seu_, features = goi, do.print = F, npcs = pcs.compute)
  ## Set k explaining >=MINPCT_VAREXPLAINED variance
  sdev=seu_@reductions$pca@stdev; sdev=100*sdev/sum(sdev);    k=1
  while(sum(sdev[1:k])<MINPCT_VAREXPLAINED){
    k=k+1
  }
  print(paste("Using",k,"PCs for tSNE..."))
  ## Variance / PC plot
  if(savePCvarPlot){  pdf(paste0(RESDIR,filesep,seu_@project.name,"_VarPerPC.pdf"))  }
  plot(sdev,xlab="PC",ylab="% variance explained",pch=20,cex=1.3,main=paste(MINPCT_VAREXPLAINED,"% k =",k))
  if(savePCvarPlot){    dev.off()  }
  ## TSNE & clustering
  seu_ <- RunTSNE(object = seu_, dims = 1:k, do.fast = do.fast, dim.embed=dim.embed,  perplexity= perplexity)
  if(clustering){
    print("Finding clusters...")
    seu_ <- FindNeighbors(object = seu_, dims = 1:k, k.param = 30, reduction = "pca" )
    seu_ <- FindClusters(object = seu_, resolution = resolution)
  }
  return(seu_)
}
