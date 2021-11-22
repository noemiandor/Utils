## Given a dataframe `dm` with genes as rows and an arbitrary number of columns, annotates each gene with information from biomart. 
## `join_id` specifies what type of gene-identifiers the rows of `dm` are (e.g. HUGO symbols, Entrez IDs, etc).
annotateFromBioMart<-function(dm,join_id,excludeSex=F,GRCh=38, dataset="hsapiens_gene_ensembl", GOregex=NULL, mart=NULL, otherCoi=NULL){ #GOregex="cellular_component:integral to membrane,^membrane,*plasma membrane*,^outer membrane*,*extracell*") {

  coi=unique(c(join_id,"ensembl_gene_id","entrezgene","hgnc_symbol" ,"chromosome_name","start_position","end_position", otherCoi)); #"affy_hg_u133a", 
  library(biomaRt)
  
  if(!is.null(GOregex)){
    coi=unique(c(coi,"namespace_1003", "name_1006", "go_id")); #"go_linkage_type", "definition_1006"
  }
  
  if(GRCh==36){
    ensembl=useMart("ENSEMBL_MART_ENSEMBL",dataset=dataset ,host="may2009.archive.ensembl.org")
  }else if (GRCh==37){
    ensembl=try(useMart("ENSEMBL_MART_ENSEMBL",dataset=dataset ,host="feb2014.archive.ensembl.org"))
    if(class(ensembl)=="try-error"){
      ensembl=useMart("ENSEMBL_MART_ENSEMBL",dataset=dataset ,host="grch37.ensembl.org")
    }
  }else if(GRCh==38){   
    # ensembl = useEnsembl(biomart="ensembl",GRCh=NULL)
    ensembl=try(useMart("ENSEMBL_MART_ENSEMBL",dataset=dataset ,host="dec2017.archive.ensembl.org"))
    if(class(ensembl)=="try-error"){
      ensembl=try(useMart("ENSEMBL_MART_ENSEMBL",dataset=dataset ,host="aug2017.archive.ensembl.org"))
      if(class(ensembl)=="try-error"){
        ensembl=useMart("ENSEMBL_MART_ENSEMBL",dataset=dataset ,host="grch38.ensembl.org")
      }
    }
  }else{
    print(paste("GRCh",GRCh,"not supported. Only supporting GRCh=36, GRCh=37 or GRCh=38. Abborting."))
    return()
  }
  mart = biomaRt::useDataset(dataset,mart=ensembl)
  
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
  x=intersect_MatlabV(toupper(rownames(dm)),toupper(genes[,join_id]) );  ##All chromosomes
  dm=dm[x$ia,,drop=F];  genes=genes[x$ib,,drop=F]
  if(excludeSex){
    genes[,"chromosome_name"]=as.numeric(gsub("Y","24",gsub("X","23",genes[,"chromosome_name"])))
    iK=which(!is.na(genes[,"chromosome_name"])); dm=dm[iK,,drop=F]; genes=genes[iK,,drop=F]
  }
  colnames(genes)=gsub("end_position","endpos",gsub("start_position","startpos",gsub("chromosome_name","chr",colnames(genes))))
  return(cbind(genes,dm))
  
}
