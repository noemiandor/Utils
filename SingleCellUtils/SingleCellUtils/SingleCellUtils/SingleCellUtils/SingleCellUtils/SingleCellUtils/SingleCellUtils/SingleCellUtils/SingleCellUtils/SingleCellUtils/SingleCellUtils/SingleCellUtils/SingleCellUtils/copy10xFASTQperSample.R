copy10xFASTQperSample<-function(){
library(plyr)
library(xlsx)
library(matlab)
seqD="/mnt/ix1/Seq_Runs"
# targetD="/mnt/noemi_temp/raw"
targetD="/mnt/ix2/nandor/noemi_temp/raw"
dm=read.xlsx("/mnt/ix2/nandor/SingleCellRnaSeq_10X/BeforeSubmission2/data/scRNAseq/SeqRunsInfo_FLproject_2017March27.xlsx",sheetName = "Indices_ATGC",row.names=T,check.names=F,stringsAsFactors=F)

for(sName in setdiff(colnames(dm),c("LPM011_1","LPM011_2")) ){
  patient=dm["Patient",sName]
  replicate=dm["Replicate",sName]
  outD=paste0(targetD,filesep,patient,filesep,patient,"_",replicate)
  dir.create(outD)
  ##
  d=list.files(seqD,pattern=paste0("*",dm["Run",sName]),full.names = T);
  d=paste0(d,filesep,dm["Flowcell",sName],"/outs/fastq_path/")
  for(index in c("Index_A",	"Index_C",	"Index_G","Index_T")){
    print(paste("Copying",index,"for",sName,"to",outD))
    f=list.files(d,pattern=dm[index,sName],full.names = T)
    ##Copy files
    for(x in f){
      y=paste0(outD,filesep,fileparts(x)$name,fileparts(x)$ext)
      file.copy(x,y)
    }
  }
}
}
