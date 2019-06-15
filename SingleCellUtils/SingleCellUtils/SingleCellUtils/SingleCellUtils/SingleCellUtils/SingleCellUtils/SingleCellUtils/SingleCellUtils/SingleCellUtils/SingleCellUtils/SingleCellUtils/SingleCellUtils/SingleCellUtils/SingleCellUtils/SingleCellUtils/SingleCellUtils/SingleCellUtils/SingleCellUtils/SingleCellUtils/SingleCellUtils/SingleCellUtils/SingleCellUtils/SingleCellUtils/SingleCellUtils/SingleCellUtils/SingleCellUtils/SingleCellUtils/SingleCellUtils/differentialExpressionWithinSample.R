differentialExpressionWithinSample<-function(d_cell,i1,i2,MINCELLS=1){
  allGI=which(apply(d_cell$mat[union(i1,i2),]>0,2,sum,na.rm=T)>=MINCELLS); ##any gene expressed in at least one cell of interest
  noNaI=which(apply(!is.na(d_cell$mat[union(i1,i2),]),2,all)); ##any gene with non-missing values across all cells of interest
  allG=intersect(allGI,noNaI);
  
  d_cell$mat=d_cell$mat[,allG]; d_cell$genes=d_cell$genes[allG]; 
  
  print(paste("# cells in group 1:",length(i1)));
  print(paste("# cells in group 2:",length(i2)));
  
  dE=matrix(NA,length(d_cell$genes),2); rownames(dE)=d_cell$genes; colnames(dE)=c("P-value","Ratio")
  cnt=0
  for (gI in 1:ncol(d_cell$mat)){
    cnt=cnt+1;
    test=try(wilcox.test(d_cell$mat[i1,gI], d_cell$mat[i2,gI]),silent=F)
    if(mod(cnt,500)==0){
      print(paste("Testing gene",cnt,"out of",length(allG),d_cell$genes[gI],"..."))
    }
    m1=mean(d_cell$mat[i1,gI]); m2=mean(d_cell$mat[i2,gI])
    dE[gI,]=c(test$p.value, (m1+0.1)/(m2+0.1))
  }  
  return(dE)
}
