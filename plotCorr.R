plotCorr<-function(vmr_dna,vmr_rna,xlab,ylab,main="",group=NULL,log='xy',textL=NULL, xlim=NULL,
                   cex=2,cex.text=0.95,col="cyan",legend="topleft",symm=F, cor.method="pearson"){
  LOGSHIFT=0.01
  te=try(cor.test(vmr_dna,vmr_rna,method=cor.method, na.rm=T))
  if(class(te)=="try-error"){ te=list(estimate=NA,p.value=NA) }
  
  vmr_dna=vmr_dna+LOGSHIFT
  vmr_rna=vmr_rna+LOGSHIFT
  
  if(is.null(xlim)){
    xlim=quantile(vmr_dna,c(0,1),na.rm=T)
  }
  ylim=quantile(vmr_rna,c(0,1),na.rm=T)
  if(symm){
    ylim<-xlim<-quantile(c(vmr_dna,vmr_rna),c(0,1), na.rm=T)
    ylim[1] <- xlim[1] <- ylim[1]-0.3*ylim[1]
    ylim[2] <- xlim[2] <- ylim[2]+0.3*ylim[1]
  }
  plot(vmr_dna,vmr_rna,log=log,pch=20,col=col,xlab=xlab,ylab=ylab,xlim=xlim,ylim=ylim,
       cex=cex,cex.lab=1.5,cex.axis=1.5,main=paste0(main,": r=",round(te$estimate,2),"; P=",round(te$p.value,4))); 
  #lines(c(0,1),c(0,1))
  
  if(!is.null(textL)){
    text(vmr_dna,vmr_rna,textL,col="black",cex=cex.text)
  }
  
  ##Color code by group membership
  if(!is.null(group)){
    ucols=brewer.pal(length(group),"Paired"); names(ucols)=group
    for(sName in group){
      i=grep(sName,names(vmr_dna));
      points(vmr_dna[i],vmr_rna[i],pch=20,col=ucols[sName],cex=cex)
    }
    if(!is.null(legend)){
      legend(legend,group,fill=ucols,cex=0.9,bty = "n")
    }
  }
  
  return(te)
}