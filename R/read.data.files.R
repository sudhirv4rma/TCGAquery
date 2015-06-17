#Functions to read the files containing the array data
read.array.data=function(data.dir, sdrf, out.file.prefix, file.sample.match, 
                         clin.data, out.dir, file.type=NULL, 
                         tcga.data.type=NULL, 
                         return.available.data.types=FALSE){
  if("genomic.reference" %in% colnames(sdrf)){
    file.info=data.frame(analysis.route=sdrf$analysis.route, 
                         genomic.reference=sdrf$genomic.reference, 
                         file.type=sdrf$final.protocol.name, 
                         tcga.data.type=sdrf$tcga.data.type,
                         file.name=sdrf$file.name, stringsAsFactors=FALSE)
  }else{
    file.info=data.frame(analysis.route=sdrf$analysis.route, 
                         file.type=sdrf$final.protocol.name, 
                         tcga.data.type=sdrf$tcga.data.type, 
                         file.name=sdrf$file.name, stringsAsFactors=FALSE)
  }
  
  file.info=unique(file.info)
  #Remove file names that are not in the data directory
  qw=which(!(file.info$file.name %in% dir(data.dir)))
  if(length(qw)>0){
    file.info=file.info[-qw,]
    warning(length(qw), " data files were not found in the specified data directory")
  }
  
  #Select only the requested data types
  if(!is.null(file.type) | !is.null(tcga.data.type)){
    selected=rep(FALSE, nrow(file.info))
    if(!is.null(file.type))
      selected[file.info$file.type %in% file.type]=TRUE
    if(!is.null(tcga.data.type))
      selected[file.info$tcga.data.type %in% tcga.data.type]=TRUE
    if(sum(selected)==0){
      warning("No data files of specified type found")
      return(NULL)
    }
    file.info=file.info[which(selected), ]
  }
  
  #Group the files according to the data type
  data.types=unique(file.info[,-ncol(file.info)])
  if(return.available.data.types)
    return(data.types)
  
  data.type.id=match(apply(file.info[,-ncol(file.info)], 1, paste, collapse="|"), 
                     apply(data.types, 1, paste, collapse="|"))  
 
  data.types$out.file.name=create.filenames(out.file.prefix, data.types)
  
  #compiled.data=list()
  for(i in 1:nrow(data.types)){
    sel=which(data.type.id==i)
    #tcga.data.type=unique(file.info$tcga.data.type[sel])
    #genomic.reference=NA
    if("genomic.reference" %in% colnames(file.info))
      genomic.reference=unique(file.info$genomic.reference[sel])
    file.info.sel=file.info[sel,]
    res=read.data.files(data.dir, file.info.sel$file.name, data.types$file.type[i])
    #file.type=unique(file.info$file.type[sel])
    #tcga.data.type=tcga.data.type
    #if(!is.na(genomic.reference))
    #  res$genomic.reference=genomic.reference
    
    X=res$X
    annot=res$annot
    probes=res$probes
    #Match up dat with file.sample.match, sdrf and clin.data2
    qm=match(colnames(X), file.sample.match$filename)
    file.sample.match2=file.sample.match[qm,]
    qm=match(colnames(X), sdrf$file.name)
    sdrf2=sdrf[qm,]
    if(!is.null(clin.data)){
      biospecimen=clin.data$biospecimen[qm,]
      qm=match(biospecimen$patient.uuid, clin.data$patient$patient.uuid)
      patient=clin.data2$patient[qm,]
      expd=cbind(patient, biospecimen, sdrf2, file.sample.match2)
    }else{
      expd=cbind(sdrf2, file.sample.match2)
    }
    #analysis.route=names(x)[k]
    #file.type=x[[k]]$file.type
    data.info=as.list(data.types[i,])
    save(X, annot, expd, probes, data.info,
         file=file.path(out.dir, data.types$out.file.name[i]))
    
#     compiled.data[[i]]=res
  }
  #names(compiled.data)=names(file.groups)
  return(data.types)
}

create.filenames=function(out.file.prefix, data.types){
  out.file.name=paste(out.file.prefix, data.types$tcga.data.type, sep=".")
  if(any(duplicated(out.file.name))){
    #Try adding the genomic reference to distinguish between the duplicates
    if("genomic.reference" %in% names(data.types)){
      out.file.name=paste(out.file.name, data.types$genomic.reference, sep=".")
    }
    if(any(duplicated(out.file.name))){
      #Add 1, 2, 3 etc to distinguish between the duplicates
      new.file.names=c()
      for(i in 1:length(out.file.name))
        new.file.names[i]=paste(out.file.name, sum(out.file.name[1:i]==out.file.name[i]), sep="-")
      out.file.name=new.file.names
    }
  }
  out.file.name=paste(out.file.name, "Rdata", sep=".")
  return(out.file.name)
}


read.data.files=function(data.dir, file.names, file.type){
  options(stringsAsFactors=FALSE)
  
  data(data.file.formats)
  #Use the format specifications to load and check the data in file
  qw=which(file.type==names(data.file.formats$pattern))
  if(length(qw)==0)
    stop("Format for ", file.type, " not found (e.g. ", file.names[1], ")\n")
  
  pattern=data.file.formats$pattern[[qw]]
  column.nums=data.file.formats$column.nums[[qw]]
  data.names=data.file.formats$names[[qw]]
  skip=data.file.formats$skip[qw]
  col.classes=data.file.formats$col.classes[[qw]]
  combine.type=data.file.formats$combine.type[qw]
  data.col=data.file.formats$data.col[qw]
  annot.cols=data.file.formats$annot.cols[[qw]]
  
  if(combine.type=="cbind")
    cbind.using=data.file.formats$cbind.using[qw]
  
  X=NULL
  probes=NULL
  annot=NULL
  for(file.name in file.names){
    d=read.delim(file.path(data.dir, file.name), skip=skip, colClasses="character")
    #Parse the columns of d
    parsed=NULL
    for(i in 1:ncol(d)){
      if(pattern[i]==""){
        parsed=cbind(parsed, d[,i])
      }else{
        dat=split2matrix(d[,i], split=pattern[i])
        #If the split column has fewer columns than expected, pad with missing data
        if(ncol(dat)<column.nums[i])
          dat=cbind(dat, array(data=NA, dim=c(nrow(dat), column.nums[i]-ncol(dat))))
        parsed=cbind(parsed, dat)
      }
    }
    parsed=as.data.frame(parsed)
    qw=which(data.names=="ignore")
    if(length(qw)>0){
      parsed=parsed[,-qw]
      colnames(parsed)=data.names[-qw]
    }else{
      colnames(parsed)=data.names
    }
    
    for(i in 1:ncol(parsed))
      parsed[,i]=as(parsed[,i], col.classes[i])
    
    if(combine.type=="rbind"){
      parsed$file.name=file.name
      X=rbind(X, parsed)
    }else if(combine.type=="cbind"){
      if(grepl(cbind.using, pattern="^col[[:digit:]]$")){
        unique.id=d[,as.integer(gsub(cbind.using, pattern="^col", replacement=""))]
      }else{
        unique.id=parsed[,colnames(parsed)==cbind.using]
      }
      
      if(is.null(probes)){
        probes=unique.id
        annot=parsed[,annot.cols]
        X=cbind(X, parsed[,data.col])
      }else{
        qm=match(probes, unique.id)
        parsed1=parsed[qm,]
        unmatched=setdiff(1:length(unique.id), qm)
        if(length(unmatched)>0){
          parsed2=parsed[unmatched,]
          parsed1=rbind(parsed1, parsed2)
          probes2=c(probes, unique.id[unmatched])
          annot=rbind(annot, parsed2[,annot.cols])
          qm=match(probes2, probes)
          X=X[qm,]
          probes=probes2
        }
        X=cbind(X, parsed1[,data.col])
      }
    }
    cat(which(file.names==file.name), "of", length(file.names), "data files read\n")
  }
  if(combine.type=="cbind")
    colnames(X)=file.names
  
  return(list(probes=probes, annot=annot, X=X))
}
