group10XcellsIntoTSNEclusters<-function(mtxFile,NCELLS=1000){
  library("Rtsne")
  library("Matrix")
  library("readr")
  d=list(mat=t(readMM(mtxFile)), genes=read_tsv(gsub(".mtx", "_genes.txt",mtxFile), col_names=F, col_types='c')$X1, barcodes=read_tsv(gsub(".mtx", "_barcodes.txt",mtxFile), col_names=F, col_types='c')$X1) # Note: I transposed the matrix here, so rows are barcodes, and columns are genes.
  d_cell<-.filter_dataset(d,NCELLS)
  d_cell$mtxFile=mtxFile
  ##Tsne
  print("Finding principal components...")
  pca<-.do_pca(d_cell)
  tsne<-Rtsne(pca$pca$x[,1:50],pca=F) 
  d_cell$tsne=tsne;
  return(d_cell)
}


.filter_dataset<-function(d, n_cell) {
  ##10X code (from Grace)
  # d is the list of matrix, genes and barcode
  # This function selects the top n_cell barcodes
  m <- d$mat
  g <- d$genes
  b <- d$barcodes
  
  bc_sum <- rowSums(m)
  bc_sum_ord <- order(bc_sum, decreasing=T)
  use_bc <- head(bc_sum_ord, n_cell)
  
  m_cell <- m[use_bc,]
  
  list(mat=m_cell, genes=g, barcodes=b[use_bc])
}

.do_pca<-function(x) {
  ##10X code (from Grace)
  use_genes <- which(colSums(x$mat) > 0)
  m <- x$mat[,use_genes]
  bc_tot <- rowSums(m)
  median_tot <- median(bc_tot)
  m <- sweep(m, 1, median_tot/bc_tot, '*')
  m <- log(1+m)
  m <- sweep(m, 2, colMeans(m), '-')
  m <- sweep(m, 2, apply(m, 2, sd), '/')
  pca <- prcomp(as.matrix(m))
  list(pca=pca, use_genes=use_genes)
}
