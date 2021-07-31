annotateFromBioMart<-function(dm,join_id,excludeSex=F,GRCh=38, GOregex=NULL, mart=NULL, otherCoi=NULL){ #GOregex="cellular_component:integral to membrane,^membrane,*plasma membrane*,^outer membrane*,*extracell*") {
  ## Given a dataframe `dm` with genes as rows and an arbitrary number of columns, annotates each gene with information from biomart. 
  ## `join_id` specifies what type of gene-identifiers the rows of `dm` are (e.g. HUGO symbols, Entrez IDs, Go terms, etc).
  
  coi=unique(c(join_id,"ensembl_gene_id","entrezgene","hgnc_symbol" ,"affy_hg_u133a", "chromosome_name","start_position","end_position", otherCoi))
  library(biomaRt)
  
  if(!is.null(GOregex)){
    coi=unique(c(coi,"namespace_1003", "name_1006", "go_id")); #"go_linkage_type", "definition_1006"
  }
  
  if(GRCh==36){
    ensembl=useMart("ENSEMBL_MART_ENSEMBL",dataset="hsapiens_gene_ensembl" ,host="may2009.archive.ensembl.org")
  }else if (GRCh==37){
    ensembl=try(useMart("ENSEMBL_MART_ENSEMBL",dataset="hsapiens_gene_ensembl" ,host="feb2014.archive.ensembl.org"))
    if(class(ensembl)=="try-error"){
      ensembl=useMart("ENSEMBL_MART_ENSEMBL",dataset="hsapiens_gene_ensembl" ,host="grch37.ensembl.org")
    }
  }else if(GRCh==38){   
    # ensembl = useEnsembl(biomart="ensembl",GRCh=NULL)
    ensembl=try(useMart("ENSEMBL_MART_ENSEMBL",dataset="hsapiens_gene_ensembl" ,host="dec2017.archive.ensembl.org"))
    if(class(ensembl)=="try-error"){
      ensembl=try(useMart("ENSEMBL_MART_ENSEMBL",dataset="hsapiens_gene_ensembl" ,host="aug2017.archive.ensembl.org"))
      if(class(ensembl)=="try-error"){
        ensembl=useMart("ENSEMBL_MART_ENSEMBL",dataset="hsapiens_gene_ensembl" ,host="grch38.ensembl.org")
      }
    }
  }else{
    print(paste("GRCh",GRCh,"not supported. Only supporting GRCh=36, GRCh=37 or GRCh=38. Abborting."))
    return()
  }
  mart = biomaRt::useDataset("hsapiens_gene_ensembl",mart=ensembl)
  
  genes<-biomaRt::getBM( coi, mart=mart)
  ##Deal with GO terms/ categories of interest
  if(!is.null(GOregex)){
    tmp=strsplit(GOregex,":")[[1]]
    GOcategory=tmp[1]; GOterms=strsplit(tmp[2],",")[[1]]
    genes=genes[genes$namespace_1003==GOcategory | genes$namespace_1003=="",]
    colnames(genes)=gsub("name_1006",GOcategory,colnames(genes))
    
    ii=c();
    for(GOterm in GOterms){
      ii=c(ii,grep(GOterm,genes[,GOcategory]))
    }
    genes=genes[sort(unique(ii)),]
  }
  ##Done with GO
  x=intersect_MatlabV(rownames(dm),genes[,join_id] );  ##All chromosomes
  dm=dm[x$ia,,drop=F];  genes=genes[x$ib,,drop=F]
  if(excludeSex){
    genes[,"chromosome_name"]=as.numeric(gsub("Y","24",gsub("X","23",genes[,"chromosome_name"])))
    iK=which(!is.na(genes[,"chromosome_name"])); dm=dm[iK,,drop=F]; genes=genes[iK,,drop=F]
  }
  colnames(genes)=gsub("end_position","endpos",gsub("start_position","startpos",gsub("chromosome_name","chr",colnames(genes))))
  return(cbind(genes,dm))
  
}


getGenesInvolvedIn<-function(pathway, dbs="ReactomePA"){
  library(ReactomePA)
  library(igraph)
  goi <- ReactomePA::viewPathway(pathway)$data$name
  goi=gsub("^SYMBOL:","",goi, fixed = F)
  return(goi)
}


getAllPathways <- function(include_genes=F){
  library(ReactomePA)
  library(matlab)
  xx <- as.list(reactome.db::reactomePATHID2NAME)
  tmp=sapply(xx,function(x) strsplit(x,": ")[[1]])
  tmp=tmp[sapply(tmp,"[[",1)=="Homo sapiens"]; 
  pathways=sapply(tmp,"[[",2)  
  pathways=pathways[!duplicated(pathways)]
  
  if(include_genes){
    names=pathways[T]
    pathways=sapply(pathways, getGenesInvolvedIn)
    names(pathways)=names
  }
  return(pathways)
}


exportNetwork2CytoscapeExample<-function(){
  ## Cytoscape needs to be installed: https://cytoscape.org/
  library("RCy3")
  library("WGCNA");
  library(Matrix)
  targetP = c("TP53", "EGFR","PI3K", "MYC", "LGL1")
  M=matrix(0,length(targetP),length(targetP)); rownames(M)<- colnames(M)<-targetP
  M[sample(length(M),length(M)/3)]=sample((1:100)/100, length(M)/3)
  ## Diagonale
  M[row(M)==col(M)]=1
  ## Symmetry
  M = forceSymmetric(M)
  ##Visualize network
  net=exportNetworkToCytoscape(M,threshold = 0.1)
  colnames(net$edgeData)[1:4]=c("source","target","weight","interaction")
  colnames(net$nodeData)[1]=c("id")
  # net$nodeData$group=as.character(net$nodeData$id %in% getGeneAliases(immReceptors))
  net$edgeData=as.data.frame(as.matrix(net$edgeData), stringsAsFactors=FALSE)
  net$edgeData$score=as.numeric(net$edgeData$weight)
  net$nodeData=as.data.frame(as.matrix(net$nodeData), stringsAsFactors=FALSE)
  
  cw=createNetworkFromDataFrames(net$nodeData,net$edgeData, title="network", collection="DataFrame Example")
}


loadGEOdataset <- function(geo_id, loadRaw=F, header = T,sep="\t"){
  library(Biobase)
  library(GEOquery)
  olddir=getwd()
  # load series and platform data from GEO
  gset <- getGEO(geo_id, GSEMatrix =F, getGPL=T)
  samples=lapply(gset@gsms, function(x) paste(x@header$characteristics_ch1, x@header$description))
  if(!loadRaw){
    dat=sapply(gset@gsms, function(x) x@dataTable@table[,"VALUE"])
    rownames(dat)=gset@gsms[[1]]@dataTable@table$ID_REF
    anno=gset@gpls[[1]]@dataTable@table
    rownames(anno)=anno$ID
  }else{
    anno=NULL
    ##Load raw data
    # sf=gset@gsms[[1]]@header$supplementary_file_1    
    # sf=grep("^ftp", sf, value = T)
    o=getGEOSuppFiles(geo_id)
    setwd(geo_id)
    if(!isempty(grep(".tar",rownames(o)))){
      untar(grep(".tar",rownames(o),value=T))
      files=untar(grep(".tar",rownames(o),value=T), list = T)
      sapply(files,function(x) try(gunzip(x)))
      files=gsub(".gz$","",files)
    }else{
      try(sapply(grep(".gz",rownames(o),value=T),gunzip))
      files=gsub(".gz$","",rownames(o))
    }
    dat=lapply(files, function(x) read.table(x,sep=sep,check.names = F, stringsAsFactors = F, header = header))
    names(dat)=files
    ##Column names
    tmp=sapply(gset@gsms, function(x) grep("Supplementary_files_format_and_content", x@header$data_processing, value=T))
    if(!isempty(grep("columns", tmp))){
      tmp=sapply(strsplit(tmp, "columns:"),"[[",2)
      tmp=sapply(tmp, function(x) strsplit(x,"\""))
      tmp=sapply(tmp, function(x) trimws(x))
      tmp=sapply(tmp, function(x) x[x!=""])
      columns=tmp[[1]]
      if(any(sapply(dat,ncol)!=length(columns))){
        dat=dat[sapply(dat,ncol)==length(columns)]
      }
      for(i in 1:length(dat)){ colnames(dat[[i]])=columns }
    }
    # ##Sample info
    # stype=sapply(samples,function(x) x[1])
    # ab= sapply(samples,function(x) grep(sinfo,x,value=T))
    # ab= trimws(sapply(strsplit(ab,":"), "[[",2))
    # ii=1:length(stype)
    # if(!is.null(grep_sample)){ ##subset of samples are of interest
    #   ii=grep(grep_sample, stype)
    # }
    # stype=trimws(sapply(strsplit(stype,":"), "[[",2))
    # alias=paste(stype[ii], ab[ii]); names(alias)=names(stype)[ii]
    # if(any(duplicated(alias))){
    #   warning("Alias is not unique. Aborting", immediate. = T);
    #   return();
    # }
    # ##Rename data
    # idx=sapply(names(alias), function(x) grep(x,names(dat)))
    # idx=idx[sapply(idx,length)>0]
    # dat=dat[unlist(idx)]; names(dat)=alias[names(idx)]
    # if(length(dat)==1){
    #   dat = dat[[1]]
    # }else{
    #   si = sapply(samples, paste, collapse=" ++ ")
    #   names(dat)[sapply(names(si), grep, names(dat))] = si
    # }
  }
  setwd(olddir)
  return(list(samples=samples, data=dat, anno=anno))
}