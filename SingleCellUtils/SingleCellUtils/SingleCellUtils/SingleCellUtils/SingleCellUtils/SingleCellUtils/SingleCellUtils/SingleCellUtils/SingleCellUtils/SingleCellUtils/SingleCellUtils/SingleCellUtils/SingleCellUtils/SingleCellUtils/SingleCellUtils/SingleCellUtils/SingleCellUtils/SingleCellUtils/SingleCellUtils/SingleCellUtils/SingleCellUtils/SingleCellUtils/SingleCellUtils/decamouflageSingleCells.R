.decamouflageSingleCells<-function(singleCells,dm,level=2,figureName=NA, cellres=15,rev=F,maxTypeCount=NULL,minCellFreq=0,cexSymbols=0.6, cexLegend=1.4, pch=20, add=F, closeFig=T,colorScheme="Paired",bg="white"){
  #exclude low freq cell types
  if(minCellFreq>0){
    if(is.null(dim(dm))){
      dm=dm-min(dm)+1
      fr=count(dm)
    }else{
      fr=count(apply(dm,1,which.max))
    }
    fr$freq=fr$freq/sum(fr$freq)
    maxTypeCount=sum(fr$freq>=minCellFreq);
  }
  if(!is.null(maxTypeCount)){
    
    if(is.null(dim(dm))){
      dm=dm-min(dm)+1
      fr=count(dm)
      minCells=sort(fr$freq, decreasing=T)[min(maxTypeCount,nrow(fr))];
      dm[dm %in% fr$x[fr$freq<minCells]]=0
      ##Now rename clusters
      fr=count(dm)
      for(i in 1:nrow(fr)){
        dm[dm==fr$x[i]]=10000+i
      }
      dm=dm-10000
    }else{
      fr=count(apply(dm,1,which.max))
      minCells=sort(fr$freq, decreasing=T)[min(maxTypeCount,nrow(fr))];
      dm=dm[,fr$x[fr$freq>=minCells]];
    }
  }
  
  ##plot cell type similarities
  #dm is either a similarity matrix or the cluster indices
  library(RColorBrewer)
  if(is.null(dim(dm))){
    cellTypeIDs=dm
    cellTypes=as.character(cellTypeIDs)
    
    tmp=sort(cellTypeIDs,index.return=T)
    cellTypeI=tmp$ix;
    singleCells=singleCells[cellTypeI,,drop=F]
  }else{
    iR=setdiff(1:ncol(dm),unique(unlist(apply(dm,1,which.max))) ); 
    if(length(iR)>0){
      dm=dm[,-iR]; ##Remove cell types that are not encountered
    }
    tmp=colnames(dm)
    cellTypes=c();
    for (i in 1:length(tmp)){
      if(!rev){
        cellTypes[i]=substr(tmp[i],1,cellres)
      }else{
        cellTypes[i]=substr(tmp[i],nchar(tmp[i])-cellres,nchar(tmp[i]))
      }
    }
    
    cellTypeIDs=try(cellTypes2numIDs(cellTypes,level))
    dup=cellTypeIDs[duplicated(cellTypeIDs)]
    cellTypes[cellTypeIDs==dup]="Other"
    if(class(cellTypeIDs)=="try-error"){
      cellTypeIDs=charIDs2numIDs(cellTypes)
    }
    
    ii=which(!apply(is.na(dm),1,all)); dmx=dm[ii,]; 
    singleCells=singleCells[ii,]
    cellTypeI=apply(dmx,1,which.max) 
    #     tmp=sort(cellTypeI,index.return=T)
    #     cellTypeI=tmp$x
    #     singleCells=singleCells[tmp$ix,]
  }
  
  N=max((cellTypeIDs[cellTypeI]),na.rm=T)
  uCT=unique(cellTypes[cellTypeI]);  uCTid=unique(cellTypeIDs[cellTypeI])
  ##Deal with colors
  colmap=gsub("#FFFF99","cyan",brewer.pal(N, colorScheme));  
  if(any(isCellType(uCT))){
    for(ct in uCT){
      cellI=uCTid[getCellIndices(uCT,ct)]
      if(!isempty(cellI) && !is.na(cellI)){
        colmap[cellI]=getCellColors(length(cellI),ct) 
      }
    }
  }
  if(min(cellTypeIDs)>brewer.pal.info[colorScheme,"maxcolors"]){
    N_=N-brewer.pal.info[colorScheme,"maxcolors"]
    if(N_>brewer.pal.info[colorScheme,"maxcolors"]){
      colmap[(brewer.pal.info[colorScheme,"maxcolors"]+1):N]=rainbow(N_);
    }else{
      colmap[(brewer.pal.info[colorScheme,"maxcolors"]+1):N]=colmap[1:N_]
    }
  }else if(N>brewer.pal.info[colorScheme,"maxcolors"]){
    colmap[13:N]=gray.colors(100)[1:(N-brewer.pal.info[colorScheme,"maxcolors"])]
  }
  colmap=colmap[cellTypeIDs[cellTypeI]]
  ##Done with colors
  x_minmax=quantile(singleCells[,1],c(0,1));  x_minmax[1]=x_minmax[1]-abs(x_minmax[1])*0.25;
  if(!is.na(figureName) && !add){
    jpeg(filename = figureName, width=10.55, height=7.65, units="in", res=200)
  }
  mainText=paste(gsub(".tif","",figureName),": ",nrow(singleCells),"cells")
  if(add){
    points(singleCells,xlab="TSNE 1", ylab="TSNE 2", xlim=x_minmax,main=mainText,pch=pch,cex=cexSymbols,cex.axis=2.3, cex.lab =2,col=colmap); #slot(cl,"clusters")
    legend("bottomleft",unique(cellTypes[cellTypeI]),fill=unique(colmap),cex=cexLegend,pch=pch,box.lwd = 0,box.col = "white",bg = "white");
  }else{
    par(bg = bg)
    plot(singleCells,xlab="TSNE 1", ylab="TSNE 2", xlim=x_minmax,main=mainText,pch=pch,cex=cexSymbols,cex.axis=2.3, cex.lab =2,col=colmap); #slot(cl,"clusters")
    legend("topleft",unique(cellTypes[cellTypeI]),fill=unique(colmap),cex=cexLegend,pch=pch,box.lwd = 0,box.col = "white",bg = "white");
  }
  print(paste("# of Cell Types",length(unique(cellTypes[cellTypeI]))))
  print(paste("# of Cell Type IDs",length(unique(cellTypeIDs[cellTypeI]))))
  print(paste("# of Cell Type color IDs",length(unique(colmap))))
  if(!is.na(figureName) && closeFig){
    dev.off()
    print(paste("Saved decamouflaged cells under ",figureName))
  }
  return(list(cellType=cellTypes[cellTypeI],ids=cellTypeIDs,cls=cellTypeIDs[cellTypeI]))
}
