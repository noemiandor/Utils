callCNVsFromSingleCells<-function(sc){
  source('/mnt/ix2/nandor/Projects/code/RCode/scripts/SingleCellUtils/annotateWithGeneLoc.R')
  GENEWINDOW=750
  MINP=1E-8; ##P-value below which differences between P- and Q-arms are considered significant
  
  ##Normalize single cell expression
  print("Normalize single cell expression...")
  sc$mat <- ifelse(sc$mat==0,NaN,sc$mat);
  sc$mat=log2(sc$mat+0.1)
  mpg=apply(sc$mat,2,mean,na.rm=T); #center gene across single cells
  sc$mat<- sweep(sc$mat,2,mpg)
  hkI=match(c("GAPDH","ACTB"),sc$genes); nohkI=setdiff(1:length(sc$genes),hkI)
  hkMean=apply(sc$mat[,hkI],1,mean,na.rm=T); # subtracting the average expression of housekeeping genes (GAPDH, ACTB) from all other genes
  sc$mat<- sweep(sc$mat[,nohkI],1,hkMean)
  print("Done")
  minv=min(min(sc$mat,na.rm=T))
  sc$mat=sc$mat-minv
  
  ##Prep gene coor
  sc=annotateWithGeneLoc(sc);
  
  ##Estimate CNV of each gene
  print(paste("Estimating copy number for each of ",ncol(sc$mat),"genes and ",nrow(sc$mat),"cells"))
  cnv=matrix(NA,nrow(sc$mat),ncol(sc$mat))
  p=(GENEWINDOW+1)
  cnv[,p]=apply(sc$mat[,(p-GENEWINDOW):(p+GENEWINDOW)],1,sum,na.rm=T)
  while(p<(ncol(sc$mat)-(GENEWINDOW+1))){
    cnv[,(p+1)]=cnv[,p]
    ii=which(!is.na(sc$mat[,(p-GENEWINDOW)]))
    if(!isempty(ii)){
      cnv[ii,(p+1)]=cnv[ii,p] - sc$mat[ii,(p-GENEWINDOW)]
    }
    ii=which(!is.na(sc$mat[,(p+(GENEWINDOW+1))]))
    if(!isempty(ii)){
      cnv[ii,(p+1)]=cnv[ii,p] + sc$mat[ii,(p+(GENEWINDOW+1))]
    }
    p=p+1;
  }
  cnv=cnv/(1+GENEWINDOW*2);
  
  ##Calculate whole-chromosome-arm CNVs per cell
  cen=read.table('/mnt/ix2/nandor/Projects/GenHetAcrossCancers/data/Annotation/hg19_Centromeres.txt',header = T); ##HG19 centromere info
  gpc=apply(sc$mat>0,1,sum,na.rm=T)
  chrArmCNVs=matrix(NA, nrow(sc$mat),22); ##For each chr and cell: 1 <=> significant difference between p- and q-arm; 0 <=> no significant difference
  cells=which(gpc>300); #1:nrow(cnv)
  # cells=sample(as.numeric(cells),500)
  for(cell in cells){
    if(mod(which(cells==cell),100)==0){
      print(paste("Calculated chr-arm CNVs for cell ",which(cells==cell),'out of',length(cells),"cells"))
    }
    noNa=which(!is.na(cnv[cell,]))
    for (chr in 1:22){
      CENB=cen[cen$chr==chr,]
      iP=intersect(noNa,which(sc$geneLoc[,"chr"]==chr & sc$geneLoc[,"endpos"]<CENB[1,"startpos"]))
      iQ=intersect(noNa,which(sc$geneLoc[,"chr"]==chr & sc$geneLoc[,"startpos"]>CENB[2,"endpos"]))
      if(length(iP)>=2 && length(iQ)>=2){
        # plot(remove_outliers(cnv[cell,iQ]),col="red",ylim=c(100,300),xlim=c(0,1200)); points(remove_outliers(cnv[cell,iP]),col="blue")
        t=try(t.test(remove_outliers(cnv[cell,iP]),remove_outliers(cnv[cell,c(iP,iQ)]) ))
        if(class(t)!="try-error"){
          effSize=abs(t$estimate[1]-t$estimate[2])
          minexpr=min(t$estimate[1],t$estimate[2])
          if(!is.na(t$p.value) && t$p.value<MINP && effSize>=minexpr*0.2){
            chrArmCNVs[cell,chr]=effSize/minexpr
          }
        }
      }
    }
  }
  
  ##Center the CNV vector of each cell at zero to avoid bias towards high values for high quality cells and low values for low-quality cells. 
  mpc=apply(cnv,1,mean,na.rm=T)
  cnv<- sweep(cnv,1,mpc)
  out=list(cnv=cnv,sc=sc, chrArmCNVs=chrArmCNVs)
  
  return(out)
}


remove_outliers <- function(x, na.rm = TRUE, ...) {
  qnt <- quantile(x, probs=c(.25, .75), na.rm = na.rm, ...)
  H <- 0.5 * IQR(x, na.rm = na.rm)
  y <- x
  y[x < (qnt[1] - H)] <- NA
  y[x > (qnt[2] + H)] <- NA
  y
}