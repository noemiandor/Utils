color10XcellsByGeneExpression<-function(gene,singleCellClusters,elog=F,saveFig=T){
  grp=singleCellClusters$mat[,match(gene,singleCellClusters$genes)]
  if (elog){
    grp=log(grp+0.01)
    grp=round(grp-min(grp));
  }
  while(length(unique(grp))>20){
    grp=round(grp/20)*20
  }
  fn=NA;
  if(saveFig){
    fn=gsub(".mtx", paste("_",gene,".tif",sep=""),singleCellClusters$mtxFile);
  }
  .decamouflageSingleCells(singleCellClusters$tsne$Y,grp,figureName=fn)
}
