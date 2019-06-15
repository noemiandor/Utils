getSampleOriginVector<-function(barcodes,prefix="Replicate",sampleNames=NULL, keepReplicateInfo=TRUE,batchID=NULL){
  ## @prefix = the prefix to append to sample Origins with pure numeric IDs
  
  sampleOrigin <- batch <- matrix(length(barcodes),1)
  
  for(i in 1:length(barcodes)){
    codeLength_10X=regexpr("-",barcodes[1])[[1]];##The length of a typical 10X barcode; ##TODO: what if the length varies?
    if(codeLength_10X>0){ ##This is indeed a barcode
      codeLength_10X=codeLength_10X+1;
    }
    idX=substr(barcodes[i],codeLength_10X+2,nchar(barcodes[i]))
    if(keepReplicateInfo){
      id=idX;
    }else{
      id=substr(idX,nchar(idX)-1,nchar(idX))
    }
    if (!is.na(as.numeric(id)) && as.numeric(id)>=0 && !is.null(idX) && nchar(idX)<=5){
      id=paste(prefix,"-",id,sep="")
    }else{
      sampleI=strsplit(id,"-")[[1]][2]
      if(!is.na(as.numeric(sampleI)) && !is.null(sampleNames) && nchar(idX)<=5){
        id=paste(sampleNames[as.numeric(sampleI)],"-",id,sep="");
      }
      if(!is.na(as.numeric(sampleI)) && !is.null(batchID)){
        batch[i]=batchID[as.numeric(sampleI)];
      }
    }
    sampleOrigin[i]=id
  }
  ##Matrix representation
  sampleOriginNum=matrix(0, length(barcodes),length(unique(sampleOrigin))); colnames(sampleOriginNum)=unique(sampleOrigin)
  rownames(sampleOriginNum)=substr(barcodes,1,codeLength_10X)
  for (sO in colnames(sampleOriginNum)){
    sampleOriginNum[sampleOrigin==sO,sO]=1
  }
  
  
#   tmp=unlist(strsplit(as.character(barcodes),"-")); 
#   sampleIDs=unique(tmp[seq(2,length(tmp),by=2)])
#   sampleOrigin=matrix(0,length(barcodes),1);  rownames(sampleOrigin)=barcodes
#   for(id in sampleIDs){
#     iX=grep(paste('-',id,'$',sep=""),barcodes); 
#     sampleOrigin[iX]=as.numeric(id);
#   }
  return(list(sampleOrigin=sampleOrigin,sampleOriginNum=sampleOriginNum,batch=batch))
}