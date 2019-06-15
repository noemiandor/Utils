mergeRuns<-function(x,FILEEXT,excludeCells=NULL){
  X=list()
  for(i in x){
    indata <- t(readMM(i))
    loci <- read.table(file = gsub(FILEEXT,'_genes.tsv',i),check.names = F,stringsAsFactors = F)
    iK=which(!duplicated(loci$V2)); indata=indata[,iK,drop=F]; loci=loci[iK,,drop=F]
    barcodes <- read.table(file = gsub(FILEEXT,'_barcodes.tsv',i),check.names = F,stringsAsFactors = F)
    if(!is.null(excludeCells)){
      iK=which(!barcodes$V1 %in% excludeCells)
      print(paste("Keeping",length(iK),"out of",nrow(barcodes),"cells that are not in black list."))
      indata=indata[iK,,drop=F];  barcodes=barcodes[iK,,drop=F]
    }
    
    gpc=apply(indata>0,1,sum)
    hist(gpc,main=fileparts(i)$name,400)
    if(isempty(X)){
      X=list(mat=indata, genes=as.character(loci$V2), genes_ensembl=as.character(loci$V1),barcodes=as.character(barcodes$V1))
    }else{
      b=intersect_MatlabV(X$barcodes,barcodes$V1)
      g=intersect_MatlabV(X$genes,as.character(loci$V2))
      X$mat=X$mat[b$ia,g$ia]+indata[b$ib,g$ib]
      X$genes=X$genes[g$ia];    X$genes_ensembl=X$genes_ensembl[g$ia]
      X$barcodes=X$barcodes[b$ia]
    }
    
  }
  return(X)
}
