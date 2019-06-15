analyzeFastProjectOutput<-function(fastProjectOutDir, projectionType="tSNE10", prefix="", lookAtRegex=c("_VS_","CELL"), maxConsistencyScore=-1,pca=F,folder="Fano", perReplicateTsne=F){
  library(plyr)
  library(matlab)
  HDTTYPES=c('PCA: 1,3', 'PCA: 1,2', 'tSNE10', '  PCA: 2,3', 'tSNE30', '  MDS', '     ISOMap', 'Spectral Embedding')
  
  ###########################
  ######Parameter settings###
  pcaString=""
  if(pca){
    pcaString="_PCA"
  }
  prefix=paste(prefix,pcaString,gsub(" ","",projectionType),"_",sep="")
  
  #################################
  ######Read FastProject output####
  projF=paste(fastProjectOutDir,"/Expression/",folder,"_Filter",pcaString,"/Projections.txt",sep="")
  proj=read.table(projF,sep='\t');
  pMatrix=read.table(gsub("Projections.txt","PMatrix.txt",projF,fixed = T),sep='\t',header = T,check.names = F);
  rownames(pMatrix)=pMatrix[,1]; pMatrix=pMatrix[,-1]
  colnames(proj)=c("DimensionalityReduction","barcodes","component1","component2")
  iK=grep(projectionType,proj$DimensionalityReduction)
  proj=proj[iK,]
  rownames(proj)=proj$barcodes
  
  ##########################################
  ######Plot components by sample origin####
  sampleOrigin=getSampleOriginVector(proj$barcodes)
  components=as.matrix(proj[,grep("^component",colnames(proj))])
  decam=.decamouflageSingleCells(singleCells=components,dm=as.numeric(sampleOrigin),figureName=paste(fastProjectOutDir,"/",prefix,"SampleOrigin.tiff",sep=""))
  
  ############################################
  ######Run tsne separately on each sample####
  if(perReplicateTsne){
    print("Recalculating projection separately per sample")
    weights=read.table("weights.txt",header=T,check.names = F)
    orig=as.matrix(read.table("../LPM018a12_PBMC13_hg19.anno.mtx",header=T,check.names = F))
    goi=intersect(rownames(orig),rownames(weights))
    orig=orig[goi,]+(1-weights[goi,])
    for(sI in unique(sampleOrigin)){
      print(paste("Calculating tSNE projection for sample",sI))
      this=rownames(sampleOrigin[sampleOrigin==sI,,drop=F])
      tsne=Rtsne(t(orig[,this]),pca=F)
      proj[this,paste("Sample",sI,"_component1",sep="")]=tsne$Y[,1]
      proj[this,paste("Sample",sI,"_component2",sep="")]=tsne$Y[,2]
    }
  }
  
  #########################################################################################################
  ######Annotate single cell clusters with pathways of interest: use signature-projection consistency######
  signatureScores=read.table("Expression/SignatureScores.txt",sep="\t",header=T,check.names = F)
  rownames(signatureScores)=signatureScores[,1];
  signatureScores=signatureScores[,colnames(signatureScores) %in% rownames(proj)]
  poi=grep(lookAtRegex[1],rownames(signatureScores),value = T)
  tmp=sort(pMatrix[poi,projectionType],index.return=T)
  poi=poi[tmp$ix]; 
  poi=poi[tmp$x<=maxConsistencyScore]; ##Keep only those terms  enriched above user defined threshold
  ##cell types gathering
  cellTypes=c()
  ##mean and std
  for(p in poi){
    ct12=strsplit(p,"_VS_"); ct2=ct12[[1]][2]; ##Cell type 1
    ct1=strsplit(ct12[[1]][1],"_"); 
    ct1=paste(ct1[[1]][2:length(ct1[[1]])],collapse = "_"); ##Cell type 2
    if(!isempty(grep("CELL",ct1)) && !isempty(grep("CELL",ct2))){
      cellTypes=c(cellTypes,ct1,ct2)        
    }
  }
  cellTypes=unique(cellTypes)
  
  cmap=matrix(NA,nrow(proj),length(cellTypes)); rownames(cmap)=rownames(proj);  colnames(cmap)=cellTypes
  ##mean and std
  for(p1 in cellTypes){
    for(p2 in cellTypes){
      ps=grep(paste('_',p1,"_VS_",p2,'$',sep=""),rownames(signatureScores),value = T)
      for(p in ps){
        m=0; std=sd(as.numeric(signatureScores[p,]))
        ct1B=names(signatureScores[p,signatureScores[p,]<=m-abs(std)]);
        ct2B=names(signatureScores[p,signatureScores[p,]>=m+abs(std)]);
        cmap[ct1B,p1]=apply(cbind( cmap[ct1B,p1], as.numeric(signatureScores[p,ct1B])),1,max,na.rm=T);
        cmap[ct2B,p2]=apply(cbind( cmap[ct2B,p2], as.numeric(signatureScores[p,ct2B])),1,max,na.rm=T);
      }
    }
  }
  cmap <- ifelse(cmap<0,NaN,cmap); 
  cmap=cmap[,-apply(is.na(cmap),2,all)]
  cellTypeI=as.numeric(apply(cmap,1,which.max))
  cmap=cmap[,unique(cellTypeI[!is.na(cellTypeI)])]
  
  for(x in c('^component',paste("Sample",unique(sampleOrigin),'_component',sep=""))){
    components=as.matrix(proj[,grep(x,colnames(proj))])
    ii=which(apply(!is.na(components),1,all));
    tmp=gsub('component','',gsub('^','',x,fixed=T))
    decam=.decamouflageSingleCells(singleCells=components[ii,],dm=cmap[ii,],cellres = 7,rev=T,figureName=paste(fastProjectOutDir,"/",prefix,tmp,"SignatureProjection.tiff",sep="")) 
  }
  return(list(projection=proj,signature=signatureScores))
}


