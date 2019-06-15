grpstats <- function(x,g,statscols,q1=0.5){
  allOut=list()
  o=matrix(NA, length(unique(g)),ncol(x) )
  rownames(o)=unique(g); colnames(o)=colnames(x);
  for(col in statscols){
    for(m in unique(g)){
      ii=which(g==m);
      v=NA;
      if(col=='mean'){
        v=apply(x[ii,,drop=F],2,mean,na.rm=T)
      }else if(col=='sum'){
        v=apply(x[ii,,drop=F],2,sum,na.rm=T)
      }else if(col=='var'){
        v=apply(x[ii,,drop=F],2,var)
      }else if(col=='max'){
        v=apply(x[ii,,drop=F],2,max,na.rm=T)
      }else if(col=='min'){
        v=apply(x[ii,,drop=F],2,min,na.rm=T)
      }else if(col=='quantile'){
        v=apply(x[ii,,drop=F],2,quantile,q1,na.rm=T)
      }else if(col=='median'){
        v=apply(x[ii,,drop=F],2,median,na.rm=T)
      }else if(col=='numel+'){##Count elements >0 
        v=apply(x[ii,,drop=F]>0,2,sum,na.rm=T)
      }else if(col=='numel_u'){##Count unique elements 
        v=apply(x[ii,,drop=F],2, function(k) length(unique(k)))
      }else if(col=='fraction+'){##Fraction of elements> 0 out of all finite elements
        v1=apply(!is.na(x[ii,,drop=F]),2,sum,na.rm=T)
        v=apply(x[ii,,drop=F]>0,2,sum,na.rm=T)/v1
      }else if(col=='maxcount'){##Most frequent value
        v1=plyr::count(x[ii])
        v=v1$x[which.max(v1$freq)]
      }else{
        v=get(col)(x[ii,,drop=F])
      }
      o[as.character(m),]=v;
    }
    allOut[[col]]=o;
  }
  return(allOut)
}