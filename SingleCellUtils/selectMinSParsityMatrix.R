selectMinSParsityMatrix<-function(EMUT_af,cellstep=1,genestep=2,MAXPERCENTMISSING=60){
  EMUT_af=t(EMUT_af)

  iCell=1:ncol(EMUT_af);
  ##Minimize matrix sparsity: --> keep only genes expressed in many cells & Cells expressing many genes
  iC=1:ncol(EMUT_af); iM=1:nrow(EMUT_af);
  percentMissing=Inf; tC=0; tM=0;
  while (percentMissing>MAXPERCENTMISSING && length(iC)>tC*2 && length(iM)>tM*2){
    tC=tC+cellstep; tM=tM+genestep;
    iM=which( apply(EMUT_af[,iC]>0,1,sum,na.rm=T) >= tC  );
    iC=which( apply(EMUT_af[iM,]>0,2,sum,na.rm=T) >= tM  );
    print(paste('Keeping only ',length(iM),' genes that are expressed in >= ',tC,' cells and ',length(iC),' cells that express >= ',tM,' genes'))
    percentMissing=100*sum(sum(EMUT_af[iM,iC]==0))/(length(iM)*length(iC));
    print(paste('Sparsity at: ',percentMissing,'% missing values'))
  }
  ## iCC=which.max(apply(!is.na(EMUT_af[iM,]),2,sum)); iM=iM[!is.na(EMUT_af[iM,iCC])]; ##Making sure there is at least 1 mutation expressed in all cells
   
  EMUT_af=EMUT_af[iM,iC];

  return(list(iC=iC,iG=iM))
}