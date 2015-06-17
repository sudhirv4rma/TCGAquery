read.cnvlowpassdnaseq.lvl3=function(data.dir){
  file.filter="(---|_)Segment\\.tsv$"
  files=dir(data.dir, pattern=file.filter)
  X=NULL
  if(length(files)>0){
    for(i in 1:length(files)){
      d=read.delim(file.path(data.dir, files[i]), na=c("null", "NULL"))
      d=data.frame(Sample=sample.name, d, check.names=F, stringsAsFactors=F)
      X=rbind(X, d)
    }    
  }else{
    warning("No data files matching pattern ", file.filter, " found in directory ", data.dir)
  }
  return(list(X=X))
}

read.cnvsnparray.lvl3=function(data.dir){
  #Select the hg19 segmentation results after germline CNVs have been removed
  file.filter="\\.nocnv_hg19\\.seg\\.txt$"
  files=dir(data.dir, pattern=file.filter)
  X=NULL
  if(length(files)>0){
    for(i in 1:length(files)){
      d=read.delim(file.path(data.dir, files[i]), na=c("null", "NULL"))
      X=rbind(X, d)
    }    
  }else{
    warning("No data files matching pattern ", file.filter, " found in directory ", data.dir)
  }
  return(list(X=X))
}
