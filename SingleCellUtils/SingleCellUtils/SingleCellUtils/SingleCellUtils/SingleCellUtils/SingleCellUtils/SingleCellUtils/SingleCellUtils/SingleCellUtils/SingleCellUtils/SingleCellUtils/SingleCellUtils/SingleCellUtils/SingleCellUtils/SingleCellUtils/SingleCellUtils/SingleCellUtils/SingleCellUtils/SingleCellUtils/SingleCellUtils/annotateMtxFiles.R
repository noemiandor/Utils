annotateMtxFiles<-function(d_cell=NULL,MATSUFF= '_hg19.mtx', NCELLS=1000,MINCELLS=30){
  library("Matrix")
  library("readr")
  #library('fpc')
  library("Rtsne")
  library("plyr")
  #root='~/Desktop/temp/'; root2="~/"; 
  root="/mnt/ix2/nandor/"; root2=root;
  source(paste(root,"SingleCellRnaSeq_10X/code/filter_dataset.R",sep=""))
  # MATSUFF= '_hg19_mRNA.mtx';
  files=list.files(".",pattern=paste('*',MATSUFF,sep=""))
  
  if(is.null(d_cell)){
    for (f in files[1:length(files)]){
      sample=gsub(MATSUFF,"",f)
      print(paste("Processing ",sample))
      #load(paste(sample,".RData",sep=""))
      d=list(mat=t(readMM(f)), genes=read_tsv(gsub(".mtx", "_genes.tsv",f), col_names=F, col_types='cc')$X2, barcodes=read_tsv(gsub(".mtx", "_barcodes.tsv",f), col_names=F, col_types='c')$X1) # Note: I transposed the matrix here, so rows are barcodes, and columns are genes.
      d_cell<-filter_dataset(d,NCELLS)
      
      dm=d_cell$mat; 
      colnames(dm)=d_cell$genes;
      rownames(dm)=d_cell$barcodes
      eGenes=apply(dm>0,2,sum);
      dm=dm[,which(eGenes>MINCELLS)]
      write.table(as.matrix(t(dm)),gsub(".mtx",".anno.mtx",f,fixed = T),sep = "\t",quote = F)
    }
  }else{
    dm=d_cell$mat; 
    colnames(dm)=d_cell$genes;
    rownames(dm)=d_cell$barcodes
    eGenes=apply(dm>0,2,sum);
    dm=dm[,which(eGenes>MINCELLS)]
    write.table(as.matrix(t(dm)),gsub(".mtx",".anno.mtx",d_cell$mtxFile,fixed = T),sep = "\t",quote = F)
  }
}
