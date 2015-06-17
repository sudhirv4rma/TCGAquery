#Create the data objects required for the library
TCGA.data.types=read.delim("TCGA.data.types.txt", sep="\t", stringsAsFactors=F)
save(TCGA.data.types, file="data/TCGA.data.types.rda")

TCGA.diseases=read.delim("TCGA.diseases.txt", sep="\t", stringsAsFactors=F)
save(TCGA.diseases, file="data/TCGA.diseases.rda")

load.formats=function(format.file){
  #Loads a file containing formatting information for various types of TCGA data files
  options(stringsAsFactors=FALSE)
  #Read a structured file that contains the formats to expect for each type of TCGA data file
  require(gsubfn)
  d=read.delim(format.file, sep="\t", header=F, stringsAsFactors=FALSE, comment.char="#")
  q1=which(d[,1]==">terms")
  q2=which(d[,1]==">formats")
  terms=d[(q1+1):(q2-1),1:2]
  colnames(terms)=c("term", "description")
  formats=d[(q2+2):nrow(d),1:8]
  colnames(formats)=c("data.type", "format", "skip", "col.classes", "combine.type", "cbind.using", "data.col", "annot.cols")
  skip=as.integer(formats$skip)
  col.classes=strsplit(formats$col.classes, split=",")
  combine.type=formats$combine.type
  cbind.using=formats$cbind.using
  data.col=as.integer(formats$data.col)
  annot.cols=sapply(strsplit(formats$annot.cols, split=","), as.integer)
  
  #Parse the formats
  columns=strsplit(formats$format, split="\\\\t")
  data.pattern=sapply(columns, function(x){
    #Replace special characters with escaped versions
    for(ch in c(".", "|","(", ")", "[", "{", "^", "$", "*", "+", "?"))
      x=gsub(x, pattern=paste("\\", ch, sep=""), replacement=paste("\\\\", ch, sep=""))
    return(gsub(x, pattern=paste(terms[,1], collapse="|"), replacement="(.+)"))
  })
  split.pattern=sapply(columns, function(x){
    #Replace special characters with escaped versions
    for(ch in c(".", "|","(", ")", "[", "{", "^", "$", "*", "+", "?"))
      x=gsub(x, pattern=paste("\\", ch, sep=""), replacement=paste("\\\\", ch, sep=""))
    
    splits=strsplit(x, split=paste(terms[,1], collapse="|"))
    splits=sapply(splits, function(y){y=unique(y)
                                      y=y[y!=""]
                                      return(paste(y, collapse="|"))})
    return(splits)
  })
  names(split.pattern)=formats$data.type
  
  #Number of final columns to expect after splitting each column in the original data
  split.column.nums=sapply(columns, function(x){
    #Replace special characters with escaped versions
    for(ch in c(".", "|","(", ")", "[", "{", "^", "$", "*", "+", "?"))
      x=gsub(x, pattern=paste("\\", ch, sep=""), replacement=paste("\\\\", ch, sep=""))
    
    splits=strsplit(x, split=paste(terms[,1], collapse="|"))
    
    return(sapply(splits, length))
  })
  names(split.column.nums)=formats$data.type
  
  data.names=list()
  for(i in 1:length(columns)){
    temp=NULL
    for(j in 1:length(columns[[i]])){
      temp=c(temp, strapply(columns[[i]][j], pattern=data.pattern[[i]][j], function(...){return(c(...))})[[1]])
    }
    data.names[[i]]=temp
  }
  names(data.names)=formats$data.type
  names(skip)=formats$data.type
  names(col.classes)=formats$data.type
  names(combine.type)=formats$data.type
  names(cbind.using)=formats$data.type
  names(data.col)=formats$data.type
  names(annot.cols)=formats$data.type
  return(list(pattern=split.pattern, column.nums=split.column.nums, names=data.names, 
              skip=skip, col.classes=col.classes,
              combine.type=combine.type, cbind.using=cbind.using,
              data.col=data.col, annot.cols=annot.cols,
              description=terms))
}
data.file.formats=load.formats("TCGA.data.file.format.txt")
save(data.file.formats, file="data/data.file.formats.rda")

