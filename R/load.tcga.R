#Load file information from TCGA downloads
get.data.types=function(tcga.download.dir){
  #Return the types and levels of data that are available in a tcga download directory
  #The download directory is the directory within which a TCGA archive has been uncompressed
  options(stringsAsFactors=FALSE)  
  res=NULL
  
  data.type=dir(tcga.download.dir)
  if(length(data.type)>0){
    for(j in 1:length(data.type)){
      if(data.type[j]=="Clinical"){
        format=dir(file.path(tcga.download.dir, data.type[j]))
        res=rbind(res, data.frame(directory=tcga.download.dir, data.type=data.type[j], platform=NA, level=NA, format=format))
      }else if(data.type[j]=="METADATA"){
        platform=dir(file.path(tcga.download.dir, data.type[j]))
        res=rbind(res, data.frame(directory=tcga.download.dir, data.type=data.type[j], platform=platform, level=NA, format=NA))
      }else if(data.type[j]=="file_manifest.txt"){
        res=rbind(res, data.frame(directory=tcga.download.dir, data.type=data.type[j], platform=NA, level=NA, format=NA))
      }else if(data.type[j]=="FILE_SAMPLE_MAP.txt"){
        res=rbind(res, data.frame(directory=tcga.download.dir, data.type=data.type[j], platform=NA, level=NA, format=NA))
      }else{
        #Check if data.type[j] is a directory
        info=file.info(file.path(tcga.download.dir, data.type[j]))
        if(info$isdir){
          platform=dir(file.path(tcga.download.dir, data.type[j]))
          if(length(platform)>0){
            for(k in 1:length(platform)){
              level=dir(file.path(tcga.download.dir, data.type[j], platform[k]), pattern="^Level_")
              if(length(level)>0)
                res=rbind(res, data.frame(directory=tcga.download.dir, data.type=data.type[j], platform=platform[k], level=level, format=NA))
            }
          }
        }
      }
    }
  }
  return(res)
}

get.data=function(tcga.download.dir, output.dir, clin.data=NULL, file.type=NULL, tcga.data.type=NULL, return.available.data.types=FALSE){
  #Load all the data and returns a list of names of Rdata files containing the data
  options(stringsAsFactors=FALSE)
  results=list()
  
  data.info=get.data.types(tcga.download.dir)

  if(nrow(data.info)==0)
    stop("No data directories found in ", tcga.download.dir)
  
    clin.data.new=NULL
    qw=which(data.info$data.type=="Clinical" & data.info$format=="XML")
    if(length(qw)==1){
      clin.directory=file.path(tcga.download.dir, data.info$data.type[qw], data.info$format[qw])
      clin.data.new=clinical.data.xml(clin.directory)
    }else{
      #Check for a Biotab directory
      #NOT FINISHED
      #       qw=which(data.info$data.type=="Clinical" & data.info$format=="Biotab")
      #       clin.directory=file.path(tcga.download.dir, data.info$data.type[qw], data.info$format[qw])
      #       clin.data.new=clinical.data.biotab(clin.directory)      
    }
    #If clinical data is already provided, combine it with the new clinical data from the download directory
    if(is.null(clin.data)){
      clin.data=clin.data.new
    }else{
      clin.data$biospecimen=unique(rbind(clin.data$biospecimen, clin.data.new$biospecimen))
      clin.data$patient=unique(rbind(clin.data$patient, clin.data.new$patient))      
    }
    if(is.null(clin.data)){
      warning("Clinical data not found")
    }
    
    #Check if there is a file sample match file, i.e. a file that matches the sample names to the raw data files
    file.sample.match=NULL
    qw=which(data.info$data.type=="FILE_SAMPLE_MAP.txt")
    if(length(qw)==1){
      file.sample.match=read.delim(file.path(tcga.download.dir, data.info$data.type[qw]))
      colnames(file.sample.match)=c("filename", "barcode")
    }
    #Add the data from the file manifest to the file.sample.match
    qw=which(data.info$data.type=="file_manifest.txt")
    if(length(qw)==1){
      d=read.delim(file.path(tcga.download.dir, data.info$data.type[qw]), na=c("n/a", "N/A"))
      qw=which(!is.na(d$Barcode))
      d=d[qw,]
      qw=which(colnames(d)=="File.Name")
      colnames(d)[qw]="filename"
      qw=which(colnames(d)=="Barcode")
      colnames(d)[qw]="barcode"
      if(is.null(file.sample.match)){
        file.sample.match=d
      }else{
        file.sample.match=merge(file.sample.match, d, by="filename", suffixes=c(".map", ".manifest"))
      }
    }
   
    #Find out if there are any SDRF (Sample and Data Relationship Format) files available 
    qw=which(!(data.info$data.type %in% c("Clinical", "FILE_SAMPLE_MAP.txt", "file_manifest.txt", "METADATA")))    
    data.types=data.info[qw, c("data.type", "platform", "level")]
    #Check if METADATA directory exists
    data.types$sdrf=NA
    for(j in 1:nrow(data.types)){
      qw=which(data.info$data.type=="METADATA" & data.info$platform==data.types$platform[j])
      if(length(qw)==1){
        #Get the names of the sdrf files for each row in data.types
        sdrf.dir=file.path(tcga.download.dir, data.info$data.type[qw], data.types$platform[j])
        file.name=dir(sdrf.dir, pattern="\\.sdrf\\.txt")
        if(length(file.name)==1){
          data.types$sdrf[j]=file.name
        }
      }
    }
    
    qw=which(is.na(data.types$sdrf))
    if(length(qw)>0){
      for(j in qw)
        warning(paste("Cannot load", tcga.download.dir, data.types$data.type[j], data.types$platform[j], data.types$level[j], "since sdrf file is missing"))
      data.types=data.types[-qw,]
    }
    
    qw=which(data.types$level!="Level_3")
    if(length(qw)>0){
      for(j in qw)
        warning(paste("Cannot load", tcga.download.dir, data.types$data.type[j], data.types$platform[j], data.types$level[j], "since only Level 3 data is supported"))
      data.types=data.types[-qw,]
    }
    
    
    #Read the data for each data type
    if(nrow(data.types)==0)
      next

    for(j in 1:nrow(data.types))
    {
      #sdrf file
      sdrf.file=file.path(tcga.download.dir, "METADATA", data.types$platform[j], data.types$sdrf[j])
      sdrf=read.sdrf(sdrf.file)      
      
      #Select the subset of clinical data that corresponds to the aliquots used in the current sdrf file
      qm=match(sdrf[,"Extract Name"], clin.data$biospecimen$aliquot.uuid)
      if(any(is.na(qm)))
        warning("Biospecimen data not available for all experiments in the sdrf")
      
      if(!is.null(clin.data)){
        clin.data2=clin.data
        clin.data2$biospecimen=clin.data2$biospecimen[qm,]
        #Keep information for only those patients whose aliquots are used in the current sdrf file
        qw=which(clin.data2$patient$patient.uuid %in% clin.data2$biospecimen$patient.uuid)
        clin.data2$patient=clin.data2$patient[qw,]
      }else{
        clin.data2=NULL
      }
      
      cat(tcga.download.dir, data.types$data.type[j], data.types$platform[j], data.types$level[j], length(unique(clin.data2$biospecimen$aliquot.uuid)), "\n")
  
      if(!is.na(data.types$level[j])){
        data.dir=file.path(tcga.download.dir, data.types$data.type[j], data.types$platform[j], data.types$level[j])
      }else{
        stop("Check!!")
      }
      data.type=data.types$data.type[j]
      platform=data.types$platform[j]
      level=data.types$level[j]

      if(length(setdiff(dir(data.dir), sdrf$file.name))>0)
        stop("More files than specified in sdrf file")
      qm=match(dir(data.dir), sdrf$file.name)
      if(any(is.na(qm)|duplicated(qm)))
        stop("Unspecified match for file in sdrf")
      sdrf=sdrf[qm,]
      
      out.file.prefix=paste(basename(tcga.download.dir), platform, level, sep=".")
      #Read in the data from the data files, keeping data from different analysis routes separate
      results[[j]]=read.array.data(data.dir, sdrf, out.file.prefix, file.sample.match, clin.data2, output.dir, file.type, tcga.data.type, return.available.data.types)
         
#       if(length(x)>0){
#         out.file.names=NULL
#         for(k in 1:length(x)){
#           file.name=paste(tcga.download.dir, platform, level, sep=".")
#           file.name=paste(c(file.name, x[[k]]$tcga.data.type), collapse=".")
#           if("genomic.reference" %in% names(x[[k]]))
#             file.name=paste(c(file.name, x[[k]]$genomic.reference), collapse=".")
#           
#           out.file.names[k]=file.name
#         }
#         if(any(duplicated(out.file.names))){
#           dups=unique(out.file.names[duplicated(out.file.names)])
#           for(k in 1:length(dups)){
#             qw=which(out.file.names==dups[k])
#             for(l in 1:length(qw))
#               out.file.names[qw[l]]=paste(out.file.names[qw[l]], l, sep=".")
#           }   
#         }
#         #Check again if there are duplicated file names
#         if(any(duplicated(out.file.names))){
#           warning("Duplicate output files")
#           browser()
#         }
#         out.file.names=paste(out.file.names, "Rdata", sep=".")
# 
#         for(k in 1:length(x)){
#           X=x[[k]]$X
#           annot=x[[k]]$annot
#           probes=x[[k]]$probes
#           #Match up dat with file.sample.match, sdrf and clin.data2
#           qm=match(colnames(X), file.sample.match$filename)
#           file.sample.match2=file.sample.match[qm,]
#           qm=match(colnames(X), sdrf$file.name)
#           sdrf2=sdrf[qm,]
#           if(!is.null(clin.data2)){
#             biospecimen=clin.data2$biospecimen[qm,]
#             qm=match(biospecimen$patient.uuid, clin.data2$patient$patient.uuid)
#             patient=clin.data2$patient[qm,]
#             expd=cbind(patient, biospecimen, sdrf2, file.sample.match2)
#           }else{
#             expd=cbind(sdrf2, file.sample.match2)
#           }
#           analysis.route=names(x)[k]
#           file.type=x[[k]]$file.type
#           save(X, annot, expd, probes, analysis.route, file.type,
#                file=file.path(output.dir, out.file.names[k]))
#           results=c(results, out.file.names[k])
#         }
#       }
      
    }
  
  return(results)
}

extract.files=function(tcga.download.dir){
  files=dir(tcga.download.dir, pattern=".tar.gz")
  tissue=gsub(files, pattern=".tar.gz", replacement="")
  for(i in 1:length(tissue)){
    tissue.dir=file.path(tcga.download.dir, tissue[i])
    gz.file=file.path(tcga.download.dir, files[i])
    tar.file=file.path(tcga.download.dir, gsub(files[i], pattern="\\.gz", replacement=""))
  
    dir.create(tissue.dir)
    system(paste("gzip -d \"", gz.file, "\"", sep=""))
    system(paste("tar -xvf \"", tar.file , "\" -C \"", tissue.dir, "\"", sep=""))
  }
}

# clinical.sdrf.match=function(clin.data, sdrf){
#   d=NULL
#   clin.empty=empty(clin.data)
#   sdrf.empty=empty(sdrf)
#   if(!clin.empty){
#     if(!sdrf.empty){
#       d=merge(clin.data, sdrf, by.x="bcr_aliquot_uuid", by.y="Extract Name")
#       #Sometimes the "Extract Name" column in clinical data contains barcodes instead of uuids
#       if(nrow(d)==0)
#         d=merge(clin.data, sdrf, by.x="bcr_aliquot_barcode", by.y="Extract Name")
#       if(nrow(d)==0)
#         warning("Could not find any matches between clinical data and sdrf")
#     }else{
#       warning("Sdrf data is empty")
#       d=clin.data
#     }
#   }else{
#     warning("Clinical data is empty")
#     if(!sdrf.empty){
#       d=sdrf
#     }else{
#       warning("Sdrf data is empty")
#     }
#   }
#   return(d)
# }
