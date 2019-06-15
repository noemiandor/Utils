findSignificantPrincipalComponents<-function (allbcells, genesToIncl=NULL,REP=100, nCores=2,stdOUT="log.findSignificantPrincipalComponents.txt"){
  MAXP=0.01;
  NCELLS=4000;
  if(!any(names(allbcells)=="cpg")){
    allbcells$cpg=apply(allbcells$mat>0,2,sum); 
  }
  iC=sample(1:nrow(allbcells$mat),NCELLS,replace = F);## random sample of cells
  iG=which(allbcells$cpg>0)
  X=allbcells;   X$mat=X$mat[iC,iG];
  X$barcodes=X$barcodes[iC]; X$genes=X$genes[iG]; 
  genesToIncl=intersect(genesToIncl,X$genes);
  ##Approach from Macosko et al.: performing REP independent randomizations of the data such that within 
  # each realization, the values along every gene of the expression matrix are randomly permuted. 
  # This operation randomizes the pairwise correlations between genes while leaving the expression distribution 
  # of every gene unchanged. 
  iGOI=1:ncol(X$mat);
  if (!is.null(genesToIncl)){
    iGOI=match(genesToIncl,X$genes);
  }
  

  #########################
  ##Pararllel processing ##
  library(parallel)
  # Initiate cluster
  if(nCores>round(detectCores()*0.5)){
    nCores=round(detectCores()*0.5);
    print(paste("Attempted to use more than 50% of the available cores. Using only ",nCores,"cores."))
  }
  print(paste("Using ",nCores," cores to find number of significant PCs"))
  print(paste("stdout and stderr connections will be redirected to ",stdOUT))
  cl <- makeCluster(nCores,outfile=stdOUT)
  ## Split jobs
  input=list()
  for(j in 1:REP){
    Y=X;
    for (i in 1:ncol(Y$mat)){
      pI=sample(1:length(iC),length(iC),replace = F);
      Y$mat[,i]=Y$mat[pI,i]
    }
    input[[j]]=list(Y=Y, genesToIncl=genesToIncl);
  }
  ## Distribute jobs
  results=clusterApply(cl,input,interface_callPCA)
  # Gather results
  pcas=list();
  eigenVs=matrix(NA,min(dim(X$mat[,iGOI])),REP)
  for(i in 1:length(results)){
    pca=results[[i]]
    pcas[[i]]=pca;
    eigenVs[,i]=pca$pca$sdev^2;
  }
  stopCluster(cl)
  closeAllConnections() 
  
  
  #########################
  ##Now calculate true PCA on un-permuted data: Significant PCs in
  # the un-permuted data were identified as those with larger eigenvalues compared to the highest eigenvalues across 
  # the 10000 randomized datasets (p < 0.01, Bonferroni corrected). 
  pcaU=.do_pca(X,genesToIncl=genesToIncl)
  pcSignificance=matrix(NA,nrow(eigenVs),1)
  for (numPC in 1:nrow(eigenVs)){
    test=t.test(pcaU$pca$sdev[numPC]^2-eigenVs[numPC,],alternative = "greater")
    pcSignificance[numPC]=test$p.value
  }
  pcSignificance_Adj=p.adjust(pcSignificance, method ="fdr");
  nSignificantComponents=max(which(pcSignificance_Adj<=MAXP));#The number of significant PCs
  print(paste("Found ",nSignificantComponents," significant PCs at FDR adjusted P-value ",MAXP))
  
  ############################################################
  ####Visualize transition between significant and non-significant PCs
  eVdiff=pcaU$pca$sdev^2-eigenVs
  col=rep("blue",nrow(eVdiff)); col[apply(eVdiff,1,median)>0]="red"
  ii=1:200
  outF="findSignificantPrincipalComponents.tiff"
  tiff(filename = outF, width=9.55, height=8.65, units="in", res=200)
  par(mai=c(1.5,1.5,1,0.5))
  boxplot(abs(t(eVdiff[ii,])),log="y",col=col[ii],xlab="PC",ylab=expression(paste("log| ",lambda["Original"]," - ",lambda["Random"]," |",sep="")),outline=F,cex.lab=1.72,xaxt="n",cex.axis=1.25)
  axis(side = 1, at = seq(10,max(ii),by=20),cex=1.25)
  legend("topright",c("positive","negative"),fill=unique(col),cex=1.8)
  dev.off()
  print(paste("Saved visualization of significant PCs under",outF))
  
  return(list(nPCs=nSignificantComponents, pValues=pcSignificance, pValues_Bonferroni=pcSignificance_Adj, eigenVs=eigenVs, pcaU=pcaU))
}


interface_callPCA<-function(varargin){
  source("/mnt/ix2/nandor/Projects/code/RCode/scripts/do_pca.R")
  Y=varargin$Y;
  genesToIncl=varargin$genesToIncl;
  
  pca=.do_pca(Y,genesToIncl=genesToIncl)
  return(pca)
}
