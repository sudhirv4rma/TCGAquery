read.somatic.mutations.lvl2=function(data.dir){
  annot.cols=c("Hugo_Symbol", "Entrez_Gene_Id", "Center",
               "Ncbi_Build", "Chrom", "Start_Position", 
               "End_Position", "Strand", "Variant_Classification", 
               "Variant_Type", "Reference_Allele",
               "Tumor_Seq_Allele1", "Tumor_Seq_Allele2", 
               "Dbsnp_Rs", "Dbsnp_Val_Status", "Match_Norm_Seq_Allele1",
               "Match_Norm_Seq_Allele2", "Tumor_Validation_Allele1", 
               "Tumor_Validation_Allele2", "Match_Norm_Validation_Allele1",
               "Match_Norm_Validation_Allele2")
  expd.cols=c("Tumor_Sample_Barcode", "Matched_Norm_Sample_Barcode", 
              "Tumor_Sample_UUID", "Matched_Norm_Sample_UUID")
  
  file.filter="\\.maf$"
  files=dir(data.dir, pattern=file.filter)
  X=annot=expd=NULL
  if(length(files)>0){
    d=c()
    for(j in 1:length(files)){    
      d=rbind(d, read.delim(file.path(data.dir, files[j]), check.names=F, na="", comment.char="#"))
    }      
    if(any(!annot.cols %in% colnames(d)) | any(!expd.cols %in% colnames(d))){
      cat("Some columns could not be found")
      browser()
    }
    annot=unique(d[,annot.cols])
    annot.p=apply(d[,annot.cols], 1, paste, collapse=";")
    expd=unique(d[,expd.cols])
    expd.p=apply(d[,expd.cols], 1, paste, collapse=";")
    X=as.matrix(table(annot.p, expd.p))
    qm=match(apply(annot, 1, paste, collapse=";"), rownames(X))
    X=X[qm,]
    qm=match(apply(expd, 1, paste, collapse=";"), colnames(X))
    X=X[,qm]
    rownames(X)=colnames(X)=NULL
  }else{
    warning("No data files matching pattern ", file.filter, " found in directory ", data.dir)
  }
  return(list(annot=annot, expd=expd, X=X))
}
