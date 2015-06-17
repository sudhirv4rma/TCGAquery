empty=function(x, na.test=TRUE, null.){
  #Tests whether the matrix, data.frame or vector x is empty or all NULL or NA
  if(is.null(dim(x))){
    #Empty or a vector
    if(length(x)==0){
      return(TRUE)
    }else if(all(is.na(x))|all(is.null(x))){
      return(TRUE)
    }
  }else{
    if(nrow(x)==0 | ncol(x)==0){
      return(TRUE)
    }else if(all(is.na(x))|all(is.null(x))){
      return(TRUE)
    }
  }
  return(FALSE)
}

barcode.size=function(barcode){
  s=strsplit(barcode, split="-")
  return(sapply(s, length))
}

barcode.trim=function(barcode, nterms){
  #Trims everything but the first nterms number of terms from the vector of barcodes
  #E.g
  #barcode="TCGA-A6-6137-F18036", nterms=2, result is "TCGA-A6"
  pattern=paste("^(([[:alnum:]]+-){", nterms-1, "}([[:alnum:]]+))(.+)", sep="")
  return(gsub(barcode, pattern=pattern, replacement="\\1"))
}

match.barcode=function(barcode1, barcode2){
  n1=unique(barcode.size(barcode1))
  n2=unique(barcode.size(barcode2))
  if(length(n1)>1|length(n2)>1)
    stop("Barcode has different numbers of terms")
  if(n1>n2){
    barcode1=barcode.trim(barcode1, n2)
  }else if(n2>n1){
    barcode2=barcode.trim(barcode2, n1)
  }
  return(match(barcode1, barcode2))
}

match.id=function(x, y, use.any=FALSE){
  #Match dataframes x and y using TCGA barcodes
  
  #I check for matches between x and y first using the aliquot barcode, then using the sample barcode
  #and finally with the patient barcode
  
  if("bcr_aliquot_barcode" %in% colnames(x)){
    x.barcode=x$bcr_aliquot_barcode
    x.barcode.type="aliquot"
  }else if("bcr_sample_barcode" %in% colnames(x)){
    x.barcode=x$bcr_sample_barcode
    x.barcode.type="sample"
  }else if("bcr_patient_barcode" %in% colnames(x)){
    x.barcode=x$bcr_patient_barcode
    x.barcode.type="patient"
  }else{
    if(!use.any){
      stop("No columns to match by in x")
    }else{
      barcode.cols=grep(colnames(x), pattern="barcode")
      if(length(barcode.cols)==0)
        stop("No columns to match by in x")
      #Choose the column with the longest barcodes
      barcode.nterms=NULL
      for(i in 1:length(barcode.cols)){
        n=unique(barcode.size(x[,barcode.cols[i]]))
        if(length(n)==1)
          barcode.nterms[i]=n
      }
      qw=which.max(barcode.nterms)
      x.barcode=x[,barcode.cols[qw]]
      x.barcode.type=colnames(x)[barcode.cols[qw]]
    }
  }
  
  if("bcr_aliquot_barcode" %in% colnames(y)){
    y.barcode=y$bcr_aliquot_barcode
    y.barcode.type="aliquot"
  }else if("bcr_sample_barcode" %in% colnames(y)){
    y.barcode=y$bcr_sample_barcode
    y.barcode.type="sample"
  }else if("bcr_patient_barcode" %in% colnames(y)){
    y.barcode=y$bcr_patient_barcode
    y.barcode.type="patient"
  }else{
    if(!use.any){
      stop("No columns to match by in y")
    }else{
      barcode.cols=grep(colnames(y), pattern="barcode")
      if(length(barcode.cols)==0)
        stop("No columns to match by in y")
  
      #Choose the column with the longest barcodes
      barcode.nterms=NULL
      for(i in 1:length(barcode.cols)){
        n=unique(barcode.size(x[,barcode.cols[i]]))
        if(length(n)==1)
          barcode.nterms[i]=n
      }
      qw=which.max(barcode.nterms)
      y.barcode=y[,barcode.cols[qw]]
      y.barcode.type=colnames(y)[barcode.cols[qw]]
    }
  }
  #Remove rows with missing barcodes
  qx=which(!is.na(x.barcode))
  x.barcode=x.barcode[qx]
  qy=which(!is.na(y.barcode))
  y.barcode=y.barcode[qy]
  
  x.barcode.nterms=unique(barcode.size(x.barcode))
  if(length(x.barcode.nterms)>1){
    warning("x barcode has different lengths")
    browser()
  }
  y.barcode.nterms=unique(barcode.size(y.barcode))
  if(length(y.barcode.nterms)>1){
    warning("y barcode has different lengths")
    browser()
  }
  #Trim the barcode with the longer number of characters to match the number of characters in the shorter barcode
  match.type=y.barcode.type
  if(x.barcode.nterms>y.barcode.nterms){
    x.barcode=barcode.trim(x.barcode, y.barcode.nterms)
  }else if(y.barcode.nterms>x.barcode.nterms){
    y.barcode=barcode.trim(y.barcode, x.barcode.nterms)
    match.type=x.barcode.type
  }
  common=intersect(x.barcode, y.barcode)
  matches=NULL
  if(length(common)>0){
    qm=match(x.barcode, y.barcode)
    matches=rbind(matches, data.frame(x.index=1:length(x.barcode), y.index=qm, match.type=match.type))
    qm=match(y.barcode, x.barcode)
    matches=rbind(matches, data.frame(x.index=qm, y.index=1:length(y.barcode), match.type=match.type))
    qw=which(apply(!is.na(matches), 1, all))
    matches=matches[qw,]
    matches=unique(matches)
  }
  matches$x.index=qx[matches$x.index]
  matches$y.index=qy[matches$y.index]
  return(matches)
}

combine.data=function(x,y,matches, all.x=FALSE, all.y=FALSE){
  combined.info=cbind(x[matches$x.index,], y[matches$y.index,])
  
  if(all.x){
    qw=setdiff(1:nrow(x), matches$x.index)
    if(length(qw)>0){
      null.row=y[1,]
      for(i in 1:ncol(null.row))
        null.row[,i]=NA
      combined.info=rbind(combined.info, cbind(x[qw,], null.row[rep(1, length(qw)),]))
    }
  }
  
  if(all.y){
    qw=setdiff(1:nrow(y), matches$y.index)
    if(length(qw)>0){
      null.row=x[1,]
      for(i in 1:ncol(null.row))
        null.row[,i]=NA
      combined.info=rbind(combined.info, cbind(null.row[rep(1, length(qw)),], y[qw,]))
    }
  }
  return(combined.info)
}

combine.same.cols=function(x, check.colname=TRUE){
  #Combine columns in x that contain the same information
  #If check.colname==TRUE, the columns are combined after confirming that the column names are the same
  #Remove a period followed by one or more numbers at the end of the column names
  #They are added by R to make the column names unique
  coln=colnames(x)
  coln=gsub(coln, pattern="\\.[[:digit:]]+$", replacement="")
  repeated.unmatched.colnames=NULL#Contains list of columns that are repeated in x, but whose data don't match
  groups=list()
  for(i in 1:(ncol(x)-1)){
    for(j in (i+1):ncol(x)){
      if((check.colname & coln[i]==coln[j]) | !check.colname){
        if(sum(x[,i]!=x[,j], na.rm=T)==0){
           if(length(groups)>0){
            for(k in 1:length(groups)){
              if(i%in% groups[[k]] | j %in% groups[[k]])
                groups[[k]]=unique(c(groups[[k]]), i, j)
            }
          }else{
            groups[[1]]=c(i,j)
          }
        }else{
          if(check.colname)
            repeated.unmatched.colnames=unique(c(repeated.unmatched.colnames, coln[i]))
        }
      }
    }
  }
  if(length(repeated.unmatched.colnames)>0)
    warning("Matching column names with mismatched data: \"", paste(repeated.unmatched.colnames, collapse="\", \""), "\"")
  
  if(length(groups)>0){
    to.remove=NULL
    for(k in 1:length(groups)){
      indx=groups[[k]]
      dat=x[,indx[1]]
      for(l in indx[-1]){
        qw=which(is.na(dat) & !is.na(x[,l]))
        if(length(qw)>0)
          dat[qw]=x[qw,l]
      }
      x[,indx[1]]=dat
      to.remove=c(to.remove, indx[-1])
    }
    x=x[,-to.remove]
  }
  return(x)
}

combine.rows=function(x, combine.by=c("patient", "sample", "aliquot")){
  if(combine.by=="aliquot" & "bcr_aliquot_barcode" %in% colnames(x)){
    barcode=x$bcr_aliquot_barcode
  }else if(combine.by=="sample" & "bcr_sample_barcode" %in% colnames(x)){
    barcode=x$bcr_sample_barcode
  }else if(combine.by=="patient" & "bcr_patient_barcode" %in% colnames(x)){
    barcode=x$bcr_patient_barcode
  }else{
    stop("Could not find column to combine")
  }
  
  require(plyr)
  
  for(i in 1:ncol(x))
    x[,i]=as.character(x[,i])
  
  res=ddply(x, "barcode",
               function(z){
                 if(nrow(z)==1){
                   return(z)
                 }else{
                   z=unique(z)
                   y=z[1,]
                   for(i in 1:ncol(z)){
                     if(length(unique(z[,i]))==1){
                       y[1,i]=unique(z[,i])
                     }else{
                       y[1,i]=paste(z[,i], collapse=";")
                     }
                   }
                   return(y)
                 }
               }
               
  )
}