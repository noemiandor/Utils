quantifyGenomicInstability<-function(X,outF, computeLOHscore=T, computeTAIscore=T){
  
  MAXADJDIST=2E6; ##Maximum distance between two segments for them to be considered adjacent
  pdf(paste0(outF,".pdf"))
  par(mfrow=c(2,2))
  
  gim=c("LOHscore","TAIscore","LSTscore","NCS_","wGII_","CNVburden_","DELburden_","AMPburden_","NCS_ps","wGII_ps","CNVburden_ps","DELburden_ps","AMPburden_ps","ploidy","bpc")
  M=matrix(NA,length(gim),length(unique(X$clones)))
  colnames(M)=unique(X$clones); rownames(M)=gim
  GIM=list()
  
  L=as.data.frame(parseGenomicLocus(rownames(X$eps)))
  L$seglength=1+L$endpos-L$startpos
  
  
  # i) the whole genome tumor loss of heterozygosity (LOH) score (55), <-- defined as the number of LOH regions >15 Mb, but less than a whole chromosome in length, within a tumour genome.
  if(computeLOHscore){
    L$chrlength=getChrLengths(as.character(L$chr))
    Y=X$eps[L$seglength>=15E6 & L$seglength<0.95*L$chrlength, names(X$clones), drop=F]
  }else{    Y = list()   }
  if(!isempty(Y)){
    GIM$LOHscore=t(as.matrix(apply(Y==1,2,sum)))
    M["LOHscore",]=grpstats(t(GIM$LOHscore),X$clones,"mean")$mean[colnames(M),]
    try(.addBoxplot(t(GIM$LOHscore),X$clones,"LOHscore"), silent = T)
  }else{
    M["LOHscore",]=NA  
  }
  
  
  # ii) the telomeric allelic imbalance (TAI) score (56) <-- Allelic imbalance was defined as any time the copy number of the two alleles were not equal, and at least one allele was present (Fig S1). 
  ## <-- Telomeric AI and telomeric CNA are defined as regions that extend to one of the sub-telomeres but do not cross the centromere. 
  ## Copy number of telomeric AI regions was defined as the mean copy number of the probes mapping to the region. 
  ## Copy loss was defined as a mean of less than 1.5 copies and copy gain was defined as a mean of greater than 2.5 copies.
  if(computeTAIscore){
    tac=getTelomeresAndCentromeres();
    tm=as.matrix(tac[tac$type=="telomere",c("chr","startpos","endpos","bin")])
    cm=as.matrix(tac[tac$type=="centromere",c("chr","startpos","endpos","bin")])
    L=assignToOverlappingSeg(as.matrix(L),tm,"bin"); colnames(L)=gsub("^bin$","telomere",colnames(L))
    L=assignToOverlappingSeg(L,cm,"bin"); colnames(L)=gsub("^bin$","centromere",colnames(L))
    L=as.data.frame(L)
    iK=which(!is.na(L[,"telomere"]) & is.na(L[,"centromere"]))
    if(!isempty(iK)){
      Y=X$eps[iK, names(X$clones),drop=F]
      GIM$TAIscore=t(as.matrix(apply(mod(Y,2),2,sum)))
      M["TAIscore",]=grpstats(t(GIM$TAIscore),X$clones,"mean")$mean[colnames(M),]
      try(.addBoxplot(t(GIM$TAIscore),X$clones,"TAIscore"), silent = T)
    }else{
      GIM$TAIscore=t(as.matrix(rep(NA,length(X$clones)))); colnames(GIM$TAIscore)=names(X$clones)
    }
  }
  
  # iii) the large-scale state transitions (LST) score which quantifies chromosomal breaks between adjacent regions of at least 10 Mb (57).
  ## <-- defined a state transition of the size S Mb if 2 adjacent chromosomal segments, each not less than S Mb in size, have different copy numbers and/or allelic contents.
  ## The number of state transitions in the tumor genomes displayed an approximately log-linear decay as a function of S (S = 3,…,20 Mb)
  Y=X$eps[, names(X$clones), drop=F]
  ix=sort(L[,'startpos']+1E9*L[,'chr'],index.return=T)$ix
  Y=Y[ix,,drop=F]; 
  GIM$LSTscore=matrix(0,1,ncol(Y)); colnames(GIM$LSTscore)=colnames(Y); ##LST score per cell
  for(i in 2:nrow(Y)){
    r=as.data.frame(parseGenomicLocus(rownames(Y)[(i-1):i]));   r$seglength=1+r$endpos-r$startpos
    c1=r$chr[1]==r$chr[2] && r$startpos[2]-r$endpos[1]<=MAXADJDIST; ##adjacent segment pair
    c2=all(r$seglength>10E6)                                        ##each segment has at least 10 Mb
    if(!c1 || !c2){
      break;
    }
    GIM$LSTscore=GIM$LSTscore+ (Y[(i-1),,drop=F]!=Y[i,,drop=F])
  }
  M["LSTscore",]=grpstats(t(GIM$LSTscore),X$clones,"mean")$mean[colnames(M),]
  try(.addBoxplot(t(GIM$LSTscore),X$clones,"LSTscore"), silent = T)
  
  # iv) ploidy
  Y=X$eps[, names(X$clones),drop=F]
  GIM$ppc=apply(Y,2,function(x) calcPloidy(cbind(L,x),cnColumn = "x")) #ploidy per cell
  GIM$ppc=t(as.matrix(GIM$ppc))
  M["ploidy",]=grpstats(t(GIM$ppc),X$clones,"mean")$mean[colnames(M),]
  try(.addBoxplot(t(GIM$ppc),X$clones,"ploidy"), silent = T)
  
  # v) CNV burden relative to diploid state & vi) CNV burden relative to ML ploidy,
  # vii) Numerical complexity score (NCS): the sum of all whole chromosome gains and losses (chromosomes with >75% of SNP copy number values higher or lower than the ploidy of the sample were counted as whole chromosome gains or losses respectively). 
  # @TODO: Multiple copy number events affecting the same chromosome were scored separately (e.g. −2 copies = 2 chromosome losses). 
  # @TODO: NCS scores were divided by 1.5 for triploid cell lines, and by 2 for tetraploid cell lines, to account for the increased likelihood of karyotypic abnormalities in polyploid genomes.
  # viii) Weighted genome instability index: for each of the 22 autosomal chromosomes, % gained and lost genomic material was calculated relative to the ploidy of the sample. The use of percentages eliminates the bias induced by differing chromosome sizes.
  # The wGII score of a sample is defined as the average of this percentage value over the 22 autosomal chromosomes.
  Y=X$eps[, names(X$clones),drop=F]
  cnv=GIM$LSTscore[,names(X$clones),drop=F]; cnv[T]=NA
  ncs <- wGII <- amp <- del <- cnv
  l=list("_"=2,"_ps"=round(median(GIM$ppc,na.rm=T)))
  cnvCoverage=cnvAbundance(L, 1:nrow(Y))
  chrlength=sapply(unique(L$chr), function(x) cnvAbundance(L[L$chr==x,,drop=F], 1:sum(L$chr==x)));      names(chrlength)=as.character(unique(L$chr))
  for(CN in names(l)){
    for(cell in names(X$clones)){
      x=rep(NA,22); names(x)=as.character(1:22)
      for(chr in unique(L$chr)){
        ii=which(L$chr==chr)
        x[as.character(chr)]=cnvAbundance(L[ii,,drop=F], which(Y[ii,cell,drop=F]!=l[[CN]]))
      }
      wGII[1,cell]=mean(x[names(chrlength)]/chrlength, na.rm=T)
      ncs[1,cell]=sum( x[names(chrlength)]/chrlength > 0.75)
      cnv[1,cell]=sum(x)/cnvCoverage
      amp[1,cell]=cnvAbundance(L, which(Y[,cell]>l[[CN]] ))/cnvCoverage
      del[1,cell]=cnvAbundance(L, which(Y[,cell]<l[[CN]] ))/cnvCoverage
    }
    GIM[[paste0("NCS",CN)]]=ncs
    M[paste0("NCS",CN),]=grpstats(t(ncs),X$clones,"mean")$mean[colnames(M),]
    try(.addBoxplot(t(ncs),X$clones,paste0("NCS",CN)), silent = T)
    GIM[[paste0("wGII",CN)]]=wGII
    M[paste0("wGII",CN),]=grpstats(t(wGII),X$clones,"mean")$mean[colnames(M),]
    try(.addBoxplot(t(wGII),X$clones,paste0("wGII",CN)), silent = T)
    GIM[[paste0("CNVburden",CN)]]=cnv; GIM[[paste0("AMPburden",CN)]]=amp; GIM[[paste0("DELburden",CN)]]=del
    M[paste0("CNVburden",CN),]=grpstats(t(cnv),X$clones,"mean")$mean[colnames(M),]
    try(.addBoxplot(t(cnv),X$clones,paste0("CNVburden",CN)), silent = T)
    M[paste0("DELburden",CN),]=grpstats(t(del),X$clones,"mean")$mean[colnames(M),]
    try(.addBoxplot(t(del),X$clones,paste0("DELburden",CN)), silent = T)
    M[paste0("AMPburden",CN),]=grpstats(t(amp),X$clones,"mean")$mean[colnames(M),]
    try(.addBoxplot(t(amp),X$clones,paste0("AMPburden",CN)), silent = T)
  }
  
  # ix) Number of breakpoints per cell
  ii = 2:nrow(X$eps)
  ##Variance across segments per cell:
  bp = sapply(colnames(X$eps), function(x) round(X$eps)[ii-1,x]!=round(X$eps)[ii,x]);
  GIM$bpc = rbind(rep(1,ncol(bp)), bp); ##First row is "start of genome" - all cells have this bp too
  rownames(GIM$bpc) = rownames(X$eps)
  gc=cnvAbundance(L,1:nrow(L))
  GIM$bpc=t(as.matrix(apply(GIM$bpc,2,sum,na.rm=T)/gc))
  M["bpc",]=grpstats(t(GIM$bpc[,names(X$clones),drop=F]),X$clones,"mean")$mean[colnames(M),]
  try(.addBoxplot(t(GIM$bpc),X$clones,"bpc"), silent = T)
  
  dev.off()
  colnames(M)=paste0("Clone",colnames(M));##TODO: move
  out=list(summary=M, GIM=GIM)
  return(out)
}

.addBoxplot<-function(x,g,feature){
  te=anova(lm(x~paste("Clone",g)))
  b=boxplot(x~paste("Clone",g),col="cyan",ylab=feature,main=paste0("Anova: P=",round(te$`Pr(>F)`[1],6)),las=3)
}