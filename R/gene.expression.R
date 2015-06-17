read.expression.genes.lvl3=function(data.dir, derived.data.files, analysis.route){
  file.filter="gene_expression_analysis\\.txt$"
  files=dir(data.dir, pattern=file.filter)
  X=annot=expd=NULL
  if(length(files)>0){
    for(i in 1:length(files)){
      d=read.delim(file.path(data.dir, files[i]), na=c("null", "NULL"))
      if(i==1)
      {
        genes=d[,2]
      }else
      {
        qm=match(genes, d[,2])
        d=d[qm,]
      }
      X=cbind(X, as.numeric(d[,3]))
      expd=c(expd, d$barcode[1])
    }
    expd=data.frame(sample=expd, stringsAsFactors=FALSE)
    annot=data.frame(gene=genes, stringsAsFactors=FALSE)
  }else{
    warning("No data files matching pattern ", file.filter, " found in directory ", data.dir)
  }
  return(list(annot=annot, expd=expd, X=X))
}


read.rnaseqv2.lvl3=function(data.dir, derived.data.files, analysis.route){
  file.filter="rsem.gene.normalized"
  files=dir(data.dir, pattern=file.filter)
  X=annot=NULL
  probes=NULL
  if(length(files)>0){
    for(i in 1:length(files)){
      d=read.delim(file.path(data.dir, files[i]), na=c("null", "NULL"))
      new.probes=d[,1]
      x=as.numeric(d[,4])
      q1=which(!(new.probes %in% probes))
      if(length(q1)>0){
        probes=c(probes, new.probes[q1])
        if(!is.null(X))
          X=rbind(X, array(data=NA, dim=c(length(q1), ncol(X))))
      }
      qm=match(probes, new.probes)
      X=cbind(X, x[qm])
    }
    annot=strsplit(probes, split="\\|")
    annot=data.frame(gene=sapply(annot, function(x){return(x[1])}), 
                     entrez.id=as.integer(gsub(sapply(annot, function(x){return(x[2])}), pattern="_calculated", replacement="")), 
                     stringsAsFactors=FALSE)
    qw=which(annot$gene=="?")
    if(length(qw)>0)
      annot$gene[qw]=NA
    
  }else{
    warning("No data files matching pattern ", file.filter, " found in directory ", data.dir)
  }
  return(list(annot=annot, files=files, X=X))
}

read.rnaseq.lvl3=function(data.dir, derived.data.files, analysis.route){
  file.filter="gene.quantification\\.txt$"
  files=dir(data.dir, pattern=file.filter)
  X=annot=NULL
  probes=NULL
  if(length(files)>0){
    for(i in 1:length(files)){
      d=read.delim(file.path(data.dir, files[i]), na=c("null", "NULL"))
      new.probes=d[,1]
      x=as.numeric(d[,4])
      q1=which(!(new.probes %in% probes))
      if(length(q1)>0){
        probes=c(probes, new.probes[q1])
        if(!is.null(X))
          X=rbind(X, array(data=NA, dim=c(length(q1), ncol(X))))
      }
      qm=match(probes, new.probes)
      X=cbind(X, x[qm])
    }
    annot=strsplit(probes, split="\\|")
    annot=data.frame(gene=sapply(annot, function(x){return(x[1])}), 
                     entrez.id=as.integer(gsub(sapply(annot, function(x){return(x[2])}), pattern="_calculated", replacement="")), 
                     stringsAsFactors=FALSE)
    qw=which(annot$gene=="?")
    if(length(qw)>0)
      annot$gene[qw]=NA
    
  }else{
    warning("No data files matching pattern ", file.filter, " found in directory ", data.dir)
  }
  return(list(annot=annot, files=files, X=X))
}
