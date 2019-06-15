removeFromEachListField<-function(X,rowIdx=NULL, colIdx=NULL, ignoreInconsistentDim=F){
  # X$mat=X$mat[iK,]; X$barcodes=X$barcodes[iK];  X$tsne$Y=X$tsne$Y[iK,];  X$gpc=X$gpc[iK];
  # dm_TREG$dm=dm_TREG$dm[iK,]; dm_TREG$P=dm_TREG$P[iK,];   dm_DMAP$dm=dm_DMAP$dm[iK,]; dm_DMAP$P=dm_DMAP$P[iK,]; 
  
  ##First remove from matrices & retrieve dimensions
  nrows=NULL; ncols=NULL;
  for(f in names(X)){
    x=X[[f]]
    print(paste("Changing list-field ",f))
    if(class(x)=="list"){
      ##Apply recursively:
      x=removeFromEachListField(x,rowIdx=rowIdx, colIdx=colIdx, ignoreInconsistentDim = ignoreInconsistentDim )
    }else if(is.null(dim(x))){
      next;
    }else{
      ##Record dimensions and check consistency
      if( (!is.null(rowIdx) && !is.null(nrows) && nrows!=nrow(x)) || (!is.null(colIdx) && !is.null(ncols) && ncols!=ncol(x)) ){
        if(ignoreInconsistentDim){
          print(paste("Field ignored due to inconsistent dimensions"))
          next
        }else{
          stop(paste("Field ",f,"contains matrix with inconsistent dimensions"))
        }
      }
      nrows=nrow(x); ncols=ncol(x);
      
      if(!is.null(rowIdx)){
        x=x[-rowIdx,,drop=F]
      }
      if(!is.null(colIdx)){
        x=x[,-colIdx,drop=F]
      }
    }
    X[[f]]=x
  }
  
  ##Next remove from vectors
  for(f in names(X)){
    x=X[[f]]
    if(!is.null(dim(x))){
      next;
    }
    
    if(!is.null(nrows) && length(x)==nrows && !is.null(rowIdx)){
      x=x[-rowIdx]
    }
    if(!is.null(ncols) && length(x)==ncols && !is.null(colIdx)){
      x=x[-colIdx]
    }
    
    X[[f]]=x
  }
  
  return(X)
}