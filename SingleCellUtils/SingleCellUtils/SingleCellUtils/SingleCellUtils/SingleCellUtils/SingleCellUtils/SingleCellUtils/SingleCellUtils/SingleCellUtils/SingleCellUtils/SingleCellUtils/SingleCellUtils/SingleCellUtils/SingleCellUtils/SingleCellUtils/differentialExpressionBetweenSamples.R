differentialExpressionBetweenSamples<-function(sampleFile1, sampleFile2){
  
  MINEXP=30; ##Min expression in all cells
  MINEPERCELL=2;
  NCELLS=1000
  d_cell1=.readAndFilterSingleCellData(sampleFile1,NCELLS)
  d_cell2=.readAndFilterSingleCellData(sampleFile2,NCELLS)
  
  expGenesI=which(apply(d_cell1$mat, 2, sum)>=MINEXP | apply(d_cell2$mat, 2, sum)>=MINEXP)
  COLS=c("P-value",paste(c(d_cell1$sample,d_cell2$sample),"_MedianExpression",sep=""),paste(c(d_cell1$sample,d_cell2$sample),"_MeanExpression",sep=""),paste(c(d_cell1$sample,d_cell2$sample),"_PercentCells>=",MINEPERCELL,sep=""))
  dE=matrix(NA,length(d_cell1$genes),length(COLS)); rownames(dE)=d_cell1$genes; colnames(dE)=COLS
  
  for (gI in expGenesI){
    test=try(wilcox.test(d_cell1$mat[,gI], d_cell2$mat[,gI]),silent=F)
    print(paste("Testing ",d_cell1$genes[gI],"..."))
    m1=mean(d_cell1$mat[,gI]); med1=median(d_cell1$mat[,gI]); 
    m2=mean(d_cell2$mat[,gI]); med2=median(d_cell2$mat[,gI])
    pE1=100*length(which(d_cell1$mat[,gI]>=MINEPERCELL))/nrow(d_cell1$mat);
    pE2=100*length(which(d_cell2$mat[,gI]>=MINEPERCELL))/nrow(d_cell2$mat);
    dE[gI,]=c(test$p.value, med1,med2, m1,m2, pE1,pE2)
  }  
  
  ii=which(!apply(is.na(dE),1,all));    dE=dE[ii,];
  x=sort(dE[,1],index.return=T);      dE=dE[x$ix,]
  write.table(dE,file=paste(d_cell1$sample,"_",d_cell2$sample,'_DifferentialExpression.txt',sep=""  ),sep="\t",quote=F)
#   plot(dE[,4],dE[,5],xlab=paste("Mean expression in ",d_cell1$sample),ylab=paste("Mean expression in ",d_cell2$sample),log='xy')
#   lines(c(0,max(dE[,4],na.rm=T)),c(0,max(dE[,4],na.rm=T)))
  tiff(filename = paste(d_cell1$sample,"_",d_cell2$sample,'_DifferentialExpression.tiff',sep=""  ), width=6.55, height=8.65, units="in", res=200)
  plot(log(dE[,4]),log(dE[,5]),pch=20, cex=0.5,xlab=paste("log Mean expression in ",d_cell1$sample),ylab=paste("log Mean expression in ",d_cell2$sample))
  lines(log(c(0.001,max(dE[,4],na.rm=T))),log(c(0.001,max(dE[,4],na.rm=T))),col='red',cex=3)
  dev.off()
  return(dE)
}

.readAndFilterSingleCellData<-function(f,nCells){
  library("Matrix")
  library("readr")
  source("~/Desktop/temp/SingleCellRnaSeq_10X/code/filter_dataset.R")
  # rows are barcodes, and columns are genes.
  d=list(mat=t(readMM(f)), genes=read_tsv(gsub(".mtx", "_genes.txt",f), col_names=F, col_types='c')$X1, barcodes=read_tsv(gsub(".mtx", "_barcodes.txt",f), col_names=F, col_types='c')$X1) 
  
  d_cell<-filter_dataset(d,nCells)
  d_cell$sample=gsub("_mex_","",gsub("_hg19_mRNA","",gsub(".mtx", "",f)))
  return(d_cell)
}