plotSingleWellCounts <- function(f){
  library(xlsx)
  MINCNT=20
  # f = "~/Downloads/Cell Sorter.xlsx"
  wb <- loadWorkbook(f)
  sheets <- names(getSheets(wb))
  sheets = grep("Yr20", sheets, value = T)
  
  o = list()
  for(sheet in sheets){
    o[[sheet]] = read.xlsx(f, sheetName = sheet,)
    rows = o[[sheet]][,1]
    o[[sheet]] = o[[sheet]][,-1]
    o[[sheet]] = apply(o[[sheet]], 2, function(x) as.numeric(gsub(">","",x)))
    rownames(o[[sheet]]) = rows
  }
  par(mfrow=c(2,3))
  sapply(sort(names(o)), function(sheet) hist(o[[sheet]], xlab="# cells",30,col="blue", border="white", main=sheet, xlim=c(0,50)))
  
  day = as.Date(gsub("Yr","",gsub("_","/",names(o))))
  cols = plyr::count(unlist(lapply(o, colnames)))
  cols = as.character(cols$x[cols$freq==max(cols$freq)])
  rows = plyr::count(unlist(lapply(o, rownames)))
  rows = as.character(rows$x[rows$freq==max(rows$freq)])
  
  par(mfrow=c(length(rows), length(cols)), mai=c(0.4,0.2,0.4,0.2))
  for(i in rows){
    G=sapply(cols, function(j) sapply(o, function(x) x[i,j]))
    clr = c("black","red")
    sapply(colnames(G), function(j) plot(day, G[,j], pch=20, col =clr[1+(max(G[,j],na.rm=T)>=MINCNT)], type = "b", main=paste(i,j), ylim=c(0,max(10,max(G[,j], na.rm=T)))))
  }
}