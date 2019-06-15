mergeSamples<-function(mtxFile1,mtxFile2,MATSUFF=".mtx", adjustRatio=T, NCELLS=NA){
  
  if(is.null(MATSUFF)){  
    tum_d=mtxFile1
    ctrl_d=mtxFile2
    mtxFile1=tum_d$mtxFile
    mtxFile2=ctrl_d$mtxFile
  }else if(MATSUFF==".mtx"){  
    tum_d=list(mat=t(readMM(mtxFile1)), genes=read_tsv(gsub(MATSUFF, "_genes.tsv",mtxFile1), col_names=F, col_types='cc')$X2, barcodes=read_tsv(gsub(".mtx", "_barcodes.tsv",mtxFile1), col_names=F, col_types='c')$X1) # Note: rows are barcodes, and columns are genes.
    ctrl_d=list(mat=t(readMM(mtxFile2)), genes=read_tsv(gsub(MATSUFF, "_genes.tsv",mtxFile2), col_names=F, col_types='cc')$X2, barcodes=read_tsv(gsub(MATSUFF, "_barcodes.tsv",mtxFile2), col_names=F, col_types='c')$X1)
  }else if (MATSUFF==".txt"){
    tum_d=read.table(mtxFile1)
    tum_d=list(mat=as.matrix(tum_d),genes=colnames(tum_d),barcodes=rownames(tum_d));
    ctrl_d=read.table(mtxFile2)
    ctrl_d=list(mat=as.matrix(ctrl_d),genes=colnames(ctrl_d),barcodes=rownames(ctrl_d));
  }
  
  ##Sample cells if NCELLS is set
  if(!is.na(NCELLS)){
    print(paste("Sampling ",NCELLS,"cells"))
    if(nrow(tum_d$mat)>NCELLS){
      ii=sort(sample(nrow(tum_d$mat),NCELLS))
      tum_d$mat=tum_d$mat[ii,]; tum_d$barcodes=tum_d$barcodes[ii]
    }
    if(nrow(ctrl_d$mat)>NCELLS){
      ii=sort(sample(nrow(ctrl_d$mat),NCELLS))
      ctrl_d$mat=ctrl_d$mat[ii,]; ctrl_d$barcodes=ctrl_d$barcodes[ii]
    }
  }
  
  #   ##Normalization
  #   tum<-filter_dataset(tum_d,NCELLS);    
  #   ctrl<-filter_dataset(ctrl_d,NCELLS);
  
  if (adjustRatio){
    # geneCov_T=rowSums(tum_d$mat>0);  geneCov_C=rowSums(ctrl_d$mat>0)
    geneCov_C=apply(ctrl_d$mat,1,mean,na.rm=T)
    geneCov_T=apply(tum_d$mat,1,mean,na.rm=T)
    ADJRATIO= median(geneCov_C)/median(geneCov_T)
    tum_d$mat=round(1*tum_d$mat*ADJRATIO)
    ctrl_d$mat=1*ctrl_d$mat
  }
  ##Normalization done
  
  ##Merge
  ctrl_d$barcodes=paste(ctrl_d$barcodes,'-0',sep='')
  tum_d$barcodes=paste(tum_d$barcodes,'-1',sep='')
  d=tum_d;
  d$barcodes=c(ctrl_d$barcodes,tum_d$barcodes)
  d$mat=rbind(ctrl_d$mat,tum_d$mat)
  
  outFmerged=paste(substr(fileparts(mtxFile2)$name,1,14),fileparts(mtxFile1)$name,MATSUFF,sep="")
  return(list(merged=d,s1=tum_d,s2=ctrl_d, s1file=mtxFile1, s2file=mtxFile2, mergedFile=outFmerged) )
}