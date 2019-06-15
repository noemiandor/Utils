correlateCellExpression2TsneCoord<-function(X, minCells=100){
  ##X - list with fields: mat, genes, barcodes and tSNE
  ##minCells - minimum number of cells per cluster
  
  library(fpc)
  print("Running dbscan clustering algorithm...")
  cl = dbscan(X$tsne$Y,eps = 1.1)
  dm=cl$cluster
  
  ##Assign clusters<minCells as outliers
  dm=dm-min(dm)+1
  fr=count(dm)
  dm[dm %in% fr$x[fr$freq<minCells]]=0; ##Cluster id 0 is reserved for outliers
  ##Now rename clusters
  fr=count(dm)
  for(i in 1:nrow(fr)){
    dm[dm==fr$x[i]]=9999+i
  }
  dm=dm-10000
  print(paste("Found",max(dm),"clusters with >=",minCells,"cells"))
  
  ##Now calculate correlation betweengene expression and cell coordinates within each cluster
  corrC=matrix(NA,length(X$genes),max(dm)); rownames(corrC)=X$genes; colnames(corrC)=1:max(dm)
  corrP=corrC
  for(c in 1:max(dm)){
    clI=which(dm==c)
    print(paste("Calculating intra-cluster correlations for cluster # ", c, "out of",max(dm)))
    for(g in 1:nrow(X$mat)){
      t1=cor.test(X$tsne$Y[clI,1],X$mat[clI,g])
      t2=cor.test(X$tsne$Y[clI,2],X$mat[clI,g])
      corrC[g,c]=max(abs(c(t1$estimate,t2$estimate)),na.rm=T)
      corrP[g,c]=min(t1$p.value,t2$p.value,na.rm=T)
    }
  }
  
  return(list(coef=corrC, p=corrP, clusters=dm))
}