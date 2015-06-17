#Download data from TCGA with specified parameters
#and returns a directory with the unzipped files
download.tcga=function(center=NULL, disease=NULL, platform=NULL,
                       platformType=NULL, level=NULL, sampleList=NULL,
                       tumorNormal=NULL,#T for unmatched tumor, TN for tumor with matched normal, NT for vice versa and N for unmatched normal
                       download.dir=NULL, force.delete=FALSE, verbose=FALSE){
  require(RCurl)
  require(XML)
  if(is.null(download.dir)){
    #Create a temporary directory to download into
    download.dir=file.path(tempdir(), disease)
  }
  
  if(file.exists(download.dir)){
    if(!force.delete){
      if(readline("Extract directory already exists. Delete all files? (y/n)")=="y")
        unlink(dir(download.dir, full.names=T), recursive=T, force=T)
    }else{
      unlink(dir(download.dir, full.names=T), recursive=T, force=T)
    }
  }else{
    dir.create(download.dir)
  }
  
  status.check.url=send.request(center, disease, platform, platformType, level, sampleList, tumorNormal, verbose)
  if(is.null(status.check.url)){
    warning("Could not send data request")
    return(NULL)
  }
  
  if(verbose)
    cat("Status check URL:", status.check.url, "\n")
  
  #Check status until an archive url is returned
  archive.url=check.status(status.check.url, verbose)

  #Download archive file
  archive.file=download.archive(archive.url, download.dir, verbose)
  if(is.null(archive.file)){
    warning("Could not download archive file")
    return(NULL)
  }
  
  uncompress.archive(archive.file, download.dir, verbose)
  return(download.dir)
}

send.request=function(center=NULL, disease=NULL, platform=NULL,
                      platformType=NULL, level=NULL, sampleList=NULL,
                      tumorNormal=NULL,#T for unmatched tumor, TN for tumor with matched normal, NT for normals with matched tumor and N for unmatched normal
                      verbose=FALSE){
  #Send the request for data to TCGA
  url=paste("http://tcga-data.nci.nih.gov/tcga/damws/jobprocess/xml?",
            ifelse(is.null(center), "", paste("center=", center, "&", sep="")),
            ifelse(is.null(disease), "", paste("disease=", disease, "&", sep="")),
            ifelse(is.null(platform), "", paste("platform=", platform, "&", sep="")),
            ifelse(is.null(platformType), "", paste("platformType=", platformType, "&", sep="")),
            ifelse(is.null(level), "", paste("level=", level, "&", sep="")),
            ifelse(is.null(tumorNormal), "", paste("tumorNormal=", tumorNormal, "&", sep="")),
            ifelse(is.null(sampleList), "", paste("sampleList=", paste(sampleList, collapse=","), "&", sep="")),            
            sep="")
  #Remove final "&"
  url=gsub(url, pattern="&$", replace="")
  if(verbose)
    cat("Querying TCGA web service with URL", url, "...")
  
  res=getURL(url)
  if(grepl(res, pattern="HTTP STATUS 204", ignore.case=T)){
    warning("No data for specified disease/platform\n", url, "\n", res, "\n")
    return(NULL)
  }else if(grepl(res, pattern="HTTP STATUS 412", ignore.case=T)){
    warning("Error specifying database parameters\n", url, "\n", res, "\n")
    return(NULL)
  }else if(grepl(res, pattern="HTTP STATUS 413", ignore.case=T)){
    warning("Connection limit reached (1 connection every 10 seconds). Try again later\n", url, "\n", res, "\n")
    return(NULL)
  }else if(grepl(res, pattern="302 Found|500 Internal Server Error", ignore.case=T)){
    warning("TCGA website is down\n", url, "\n", res, "\n")
    return(NULL)
  }
  if(verbose)
    cat("done\n")
  doc = xmlInternalTreeParse(file=res, asText=TRUE)
  status.check.url=xmlValue(getNodeSet(doc, path="/job-process/status-check-url")[[1]])
  return(status.check.url)
}

check.status=function(status.check.url, verbose=FALSE){
  #Check status of request until an archive url is available
  status.code="202"
  if(verbose)
    cat("Waiting for archive to be ready.")
  status.message=""
  while(status.code!="200"){
    Sys.sleep(30)
    res=getURL(status.check.url)
    doc=xmlInternalTreeParse(file=res, asText=TRUE)
    status.code=xmlValue(getNodeSet(doc, path="/job-process/job-status/status-code")[[1]])
    if(verbose){
      new.message=xmlValue(getNodeSet(doc, path="/job-process/job-status/status-message")[[1]])
      if(new.message!=status.message){
        status.message=new.message
        cat("\n", format(Sys.time(), "%a %b %d %H:%M:%S"), "Status:", status.message)
      }else{
        cat(".")
      }
    }
  }
  if(verbose)
    cat("done\n")
  archive.url=xmlValue(getNodeSet(doc, path="/job-process/job-status/archive-url")[[1]])
  return(archive.url)
}

download.archive=function(archive.url, download.dir, verbose=FALSE){
  if(verbose)
    cat("Archive URL:", archive.url, "\nDownloading archive...")

  archive.file=file.path(download.dir, basename(archive.url))
  quiet=!verbose
  
  #download.file(archive.url, destfile=file.path(download.dir, basename(archive.url)), mode="wb", method=download.method, quiet=quiet)
  
  #Using function "download" from package "downloader" instead of "download.file" for handling https
  if(download(archive.url, destfile=archive.file, mode="wb", quiet=quiet,
              extra=c("--no-check-certificate", "-k"))!=0){#"--no-check-certificate" for wget and "-k" for curl
    warning("Could not download archive file: ", archive.url)
    return(NULL)
  }
  
  if(verbose)
    cat("done\n")
  return(archive.file)
}

uncompress.archive=function(archive.file, destination.dir, verbose=FALSE){
  if(verbose)
    cat("Uncompressing archive...")
  archive.file=path.expand(archive.file)
  destination.dir=path.expand(destination.dir)
  
  if(grepl(archive.file, pattern="\\.tar\\.gz$")){
    untar(archive.file, exdir=destination.dir, compressed="gzip")
  }else if(grepl(archive.file, pattern="\\.tar$")){
    untar(archive.file, exdir=destination.dir)
  }
  if(verbose)
    cat("done\n")
}