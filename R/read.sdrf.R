read.sdrf=function(sdrf.file){
  require(plyr)
  sdrf=read.delim(sdrf.file, check.names=F, na="->")
  
  #Get details for each file in the sdrf
  protocol.cols=which(colnames(sdrf)=="Protocol REF")
  file.cols=which(colnames(sdrf) %in% c("Array Data File", "Derived Array Data File", "Array Data Matrix File", "Derived Array Data Matrix File", "Derived Data File"))
  if(length(file.cols)==0)
    stop("Could not find filenames associated with analysis results in sdrf")
  geno.cols=which(colnames(sdrf)=="Comment [Genome reference]")
  level.cols=which(colnames(sdrf) %in% c("Comment [TCGA Data Level]"))
  data.type.cols=which(colnames(sdrf) %in% c("Comment [TCGA Data Type]"))
  if(length(data.type.cols)==0)
    stop("Could not find TCGA data type associated with analysis results in sdrf")
  
  #Replicate the sdrf rows for each data file
  new.sdrf=NULL
  file.name=NULL
  data.level=NULL
  final.protocol.name=NULL
  genomic.reference=NULL
  tcga.data.type=NULL
  n=length(protocol.cols)
  for(i in file.cols){
    #Find out the range of columns that refer to this output file
    if(i>protocol.cols[n]){
      col.range=protocol.cols[n]:ncol(sdrf)
    }else{
      q2=which(protocol.cols[-1]>i & protocol.cols[-n]<i)
      col.range=protocol.cols[q2]:(protocol.cols[q2+1]-1)
    }
    level.col=level.cols[which(level.cols %in% col.range)]
    protocol.col=protocol.cols[which(protocol.cols %in% col.range)]
    geno.col=NA
    if(length(geno.cols)>0)
      geno.col=geno.cols[which(geno.cols %in% col.range)]
    data.type.col=data.type.cols[which(data.type.cols %in% col.range)]
    qw=which(!is.na(sdrf[,i]))
    temp=sdrf[qw,]
    file.name=c(file.name, temp[,i])
    data.level=c(data.level, temp[,level.col])
    final.protocol.name=c(final.protocol.name, temp[,protocol.col])
    if(length(geno.col)>0){
      if(!is.na(geno.col)){
        genomic.reference=c(genomic.reference, temp[,geno.col])
      }else{
        genomic.reference=c(genomic.reference, rep(NA, nrow(temp)))
      }
    }else{
      genomic.reference=c(genomic.reference, rep(NA, nrow(temp)))
    }
    
    tcga.data.type=c(tcga.data.type, temp[,data.type.col])
    if(max(col.range)<ncol(temp))
      temp[,max(col.range):ncol(temp)]=NA
    new.sdrf=rbind(new.sdrf, temp)
  }
  sdrf=new.sdrf
  sdrf$file.name=file.name
  sdrf$data.level=data.level
  sdrf$final.protocol.name=final.protocol.name
  #The format of the data files depend on the genomic reference (see load.formats())
  if(!all(is.na(genomic.reference)))
    sdrf$genomic.reference=genomic.reference
  sdrf$tcga.data.type=tcga.data.type
  
  #For each file, get the route of sample analysis
  analysis.route=ddply(sdrf[,c(protocol.cols, which(colnames(sdrf)=="file.name"))], "file.name", 
                       function(x){
                         apply(x, 2, function(y){return(paste(sort(unique(y)), collapse="|"))})
                       })
  qm=match(sdrf$file.name, analysis.route$file.name)
  analysis.route=apply(analysis.route[,-ncol(analysis.route)], 1, function(x){return(paste(x[!(is.na(x)|x=="")], collapse="\t"))})
  sdrf$analysis.route=analysis.route[qm]
    
  return(sdrf)
}