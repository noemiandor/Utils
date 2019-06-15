erinsTSNEplot <- function(marker, singleCells,transform=NA, alias="",countID="", cex=0.5,cexMarker=1.4,main="", addConcaveHull=F){ #x=map[,1], y=map[,2]){
  index=rownames(singleCells$tsne@cell.embeddings)
  
  xdat <- singleCells$tsne@cell.embeddings[,1]
  xlabel <- "TSNE 1"
  
  ydat <- singleCells$tsne@cell.embeddings[,2]
  ylabel <- "TSNE 2"
  
  if (nchar(alias)>0){
    alias=paste(" (",alias,")",sep="")
  }
  
  # Set color scale for plots
  mypal <- rev(colorRampPalette(brewer.pal(11,"Spectral")[2:10])(256))
  
  # Load data for plot
  gI=match(marker,singleCells$genes)
  dat=singleCells$mat[index,gI]
  ##Merge to 1 if multiple genes are provided
  if(is.matrix(dat)){
    dat=apply(dat,1,sum,na.rm=T)
    marker=countID
  }
  
  if(!is.na(transform)){
    if(transform=="log"){
      dat=log(dat+0.01) 
    }else if(transform=="log10"){
      dat=log10(dat+0.01) 
    }
  }

  
  
  if(quantile(dat, probs=0.98) > 0){
    datmin <- quantile(dat,probs=0.01)
    datmax <- quantile(dat,probs=0.99)
    dat[dat<=datmin] <- datmin
    dat[dat>=datmax] <- datmax
  }  else {
    datmin <- min(dat)
    datmax <- max(dat)
  }
  
  plotdata <- cbind(xdat, ydat)
  
  mycolors <- mypal[as.numeric(cut(as.matrix(dat),breaks = 256))]
  
  par(mar=c(4.1, 3.9, 4.1, 2.1),bg="white",fg="black",col.axis="black",col.lab="black",col.main="black",col.sub="black")
  layout(matrix(1:2,ncol=2), width = c(5,1),height = c(1,1))
  plot(plotdata, col=mycolors, pch=16, cex=cex, cex.axis=2.65,cex.lab=1.72,xlab=xlabel, ylab=ylabel,main=main); #, main=paste(titl,sep="\n"), sub=paste("min = ",round(datmin,2),"; max = ",round(datmax,2),sep=""))
  mtext(at = c(1.39*max(plotdata[,1],na.rm=T)), marker,cex=cexMarker)
  if(addConcaveHull){    ##Adds ConcaveHull
    for(l_ in unique(singleCells$ident)){
      ii=which(singleCells$ident==l_)
      plot_ConcaveHull(xdat[ii], ydat[ii], zoomIn=0.1)
    }
  }
  par(mar=c(6,1,6,3))
  image(1,seq(round(datmin,2),round(datmax,2),len=256),matrix(1:256,nrow=1),col=mypal,axes=FALSE,xlab="",ylab="")
  axis(2,cex.axis=2)

  return(dat)
}


