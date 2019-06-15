parseSingleCellVCFoutput<-function(f,barcodesF=NULL,snvsToKeep=NULL, min_Cells_ExpressingGENE=1, min_Cells_ExpressingMUT=1){
  if(!is.null(snvsToKeep)){
    print(paste("Parsing ",f,"into matrix of B-allele expression and total expression. Keeping only SNVs that are also present in user-provided matrix: ",snvsToKeep))
  }
  if(is.null(barcodesF)){
    barcodesF=gsub("_3_","_3_barcodes",gsub("_2_","_2_barcodes",gsub("_1_","_1_barcodes",gsub(".rg","",fileparts(f)$name))));
    barcodesF=paste(fileparts(fileparts(f)$path)$path,filesep,barcodesF,".tsv",sep="")
    print(paste("Will read barcodes info from",barcodesF))
  }
  BALLELE_IDX=2; ##The index of the B-allele read count within the sample columns of a vcf file
  MINDEPTH=1; MINCELLS=10;
  outputF=gsub(".vcf",".Ballele.sparse.expression",f,fixed = T);
  outputFdepth=gsub(".vcf",".sparse.expression",f,fixed = T);
  
  tmp=read.table( f,nrows=900,comment.char="",fill=T,header=F,flush=T)
  skip= which(as.character(tmp[,1]) %in% '#CHROM')-1; ##This is the header line
  
  EMUT_raw=read.table( f,skip=skip,header=T,comment.char="",stringsAsFactors=FALSE,skipNul = F,flush = T)
  # EMUT_raw=EMUT_raw[1:38038,]; ##TODO: remove
  #ia=which(duplicated(EMUT_raw[,c("X.CHROM","POS")])); EMUT_raw=EMUT_raw[-ia,]
  realcellBarcodes=read.table(barcodesF)
  EMUT_raw[,"X.CHROM"]=gsub("chr","",EMUT_raw[,"X.CHROM"])
  iKeep=which(EMUT_raw[,"X.CHROM"] %in% as.character(seq(1:22))); EMUT_raw=EMUT_raw[iKeep,]
  EMUT_raw[,"X.CHROM"]=as.numeric(EMUT_raw[,"X.CHROM"])
  EMUT_loci=EMUT_raw[,c("X.CHROM","POS","REF","ALT")]; colnames(EMUT_loci)=c("chr","startpos","REF","ALT");
  
  iCells=match(realcellBarcodes[,1],paste(colnames(EMUT_raw),"-1",sep=""));   iCells=iCells[!is.na(iCells)];
  if (!is.null(snvsToKeep)){
    x=.intersectRows(snvsToKeep[,c("chr","startpos")],as.matrix(EMUT_loci[,c("chr","startpos")]))
    EMUT_af=EMUT_raw[x$ib,iCells]; ##keep only real cells at mutated loci
  }else{
    EMUT_af=EMUT_raw[,iCells];
  }
  EMUT_depth=EMUT_af
  for(i in 1:ncol(EMUT_af)){
    # GT1,GQ1,GL3,DP1,DV1,SP1,PL1
    if(mod(i,100)==1){
      print(paste("Reading mutation status for cell ",i ,'out of ',ncol(EMUT_af)))
    }
    tmp=strsplit(as.character(EMUT_af[,i]),":")
    DP=as.numeric(sapply(tmp, "[[", 2));  DV=as.numeric(sapply(tmp, "[[", 3))
    af=DV/DP;
    #af[DP<MINDEPTH]=NA;
    af[DP<MINDEPTH & DV<1]=NA; ##We reject concluding that the B-allele is not present if <MINDEPTH reads cover the locus
    #af[DP==0]=0; ##We reject concluding that the B-allele is not present if <MINDEPTH reads cover the locus
    EMUT_af[,i]=af;
    EMUT_depth[,i]=DP;
  }
  print(paste("Found ",ncol(EMUT_af),"cells.",length(which(apply(is.na(EMUT_af),2,all))),"cell(s) have inconclusive B-allele status for all",nrow( EMUT_af),"loci with mutations. Cause: B-allele not found at <",MINDEPTH ,"read depth") )
  EMUT_af=as.matrix(EMUT_af);#matrix(as.numeric(unlist(EMUT_af)),nrow=nrow(EMUT_af),dimnames = list(1:nrow(EMUT_af),colnames(EMUT_af)))
  EMUT_depth=as.matrix(EMUT_depth); #matrix(as.numeric(unlist(EMUT_depth)),nrow=nrow(EMUT_depth),dimnames = list(1:nrow(EMUT_depth),colnames(EMUT_depth)))
  colnames(EMUT_af)=paste(colnames(EMUT_af),"-1",sep=""); colnames(EMUT_depth)=colnames(EMUT_af);
  
  ##Keep only loci with minimum expression across cells:
  iKeep=which(apply(EMUT_depth>0,1,sum)>=min_Cells_ExpressingGENE & apply(EMUT_af>0,1,sum,na.rm=T)>=min_Cells_ExpressingMUT);
  EMUT_depth=EMUT_depth[iKeep,];    EMUT_af=EMUT_af[iKeep,];    EMUT_loci=EMUT_loci[iKeep,];
  
  # EMUT_af[EMUT_af>0]=1; ##Don't care what AF the mutation has as long as it's expressed
  if (!is.null(snvsToKeep)){
    EMUT_af=cbind(snvsToKeep[x$ia,c("SP","%maxP","AF_Tumor","PM_B","PM","chr","startpos")],EMUT_af)
    EMUT_depth=cbind(snvsToKeep[x$ia,c("SP","%maxP","AF_Tumor","PM_B","PM","chr","startpos")],EMUT_depth)
  }else{
    lociNames=apply(EMUT_loci[,c("chr","startpos")],1,paste,collapse="_")
    rownames(EMUT_af)=lociNames;    rownames(EMUT_depth)=lociNames;
  }
  #   write.table(EMUT_af, outputF,sep="\t",quote=F,row.names = F);
  #   write.table(EMUT_depth, outputFdepth, sep="\t",quote=F,row.names = F);
  #   print(paste("B-allele expression AF and total read depth per SNV saved under ",outputF,"and", outputFdepth))
  
  ##Done parsing vcf
  return(list(EMUT_af=EMUT_af, EMUT_depth=EMUT_depth, EMUT_loci=EMUT_loci ))
}