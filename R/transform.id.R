barcode2uuid=function(barcode){
  #Convert barcode to uuid
  if(length(barcode)>500)
    stop("A maximum of 500 barcodes can be converted at a time")
  require(RCurl)
  require(XML)
  resp=getURL("https://tcga-data.nci.nih.gov/uuid/uuidws/mapping/xml/barcode/batch", 
              customrequest="POST", httpheader=c("Content-Type: text/plain"), 
              postfields=paste(barcode, collapse=","), ssl.verifypeer = FALSE)
  doc=xmlInternalTreeParse(resp, asText=T)
  nodes=getNodeSet(doc, path="/uuidBarcodeMappings/uuidMapping")
  if(length(nodes)==0)
    stop("Some of the barcodes do not match any UUIDs")
  uuid=c()
  uuid.names=c()
  for(i in 1:length(nodes)){
    children=xmlChildren(nodes[[i]])
    uuid=c(uuid, xmlValue(children[["uuid"]]))
    uuid.names=c(uuid.names, xmlValue(children[["barcode"]]))
    rm(children)
  }
  qm=match(barcode, uuid.names)
  uuid=uuid[qm]
  names(uuid)=barcode
  
  rm(nodes)
  gc()
  free(doc)
  rm(doc)
  gc()
  
  return(uuid)
}

uuid2barcode=function(uuid){
  #Convert uuid to barcode
  if(length(uuid)>500)
    stop("A maximum of 500 uuids can be converted at a time")
  require(RCurl)
  require(XML)
  resp=getURL("https://tcga-data.nci.nih.gov/uuid/uuidws/mapping/xml/uuid/batch", 
              customrequest="POST", httpheader=c("Content-Type: text/plain"), 
              postfields=paste(uuid, collapse=","), ssl.verifypeer = FALSE)
  doc=xmlInternalTreeParse(resp, asText=T)
  nodes=getNodeSet(doc, path="/uuidBarcodeMappings/uuidMapping")
  if(length(nodes)==0)
    stop("Some of the barcodes do not match any UUIDs")
  barcode=c()
  barcode.names=c()
  for(i in 1:length(nodes)){
    children=xmlChildren(nodes[[i]])
    barcode=c(barcode, xmlValue(children[["barcode"]]))
    barcode.names=c(barcode.names, xmlValue(children[["uuid"]]))
    rm(children)
  }
  qm=match(uuid, barcode.names)
  barcode=barcode[qm]
  names(barcode)=uuid
  
  rm(nodes)
  gc()
  free(doc)
  rm(doc)
  gc()
  
  return(barcode)
}
# uuid=c("a2f14743-13f1-4480-859b-6e22435420a4",
#        "128a3673-6dfc-4a62-a14a-4a7635d8ded4",
#        "0c57bfeb-3dc5-4b0b-b27a-7c031803f4ad",
#        "06e73696-510e-48fe-975b-93908e8ca545",
#        "04ffbe69-badc-4636-ba81-d0efa3c2196d",
#        "c08b110a-b36e-4042-b599-fa9e088b4ae6",
#        "67faec27-840a-491c-bd18-4e3a43f36899",
#        "209db446-9fda-40bc-879b-5615409908bd",
#        "8606e511-48f7-4880-9da8-644d2f700403",
#        "8e195b3c-4512-4887-9191-faaac56c1c65")
# barcode=uuid2barcode(uuid)
# barcode=c("TCGA-02-0071-01A-01R-0202-01",
#           "TCGA-02-0086-01A-01R-0202-01",
#           "TCGA-02-0115-01A-01R-0202-01",
#           "TCGA-02-0075-01A-01R-0202-01",
#           "TCGA-02-0080-01A-01R-0202-01",
#           "TCGA-02-0089-01A-01R-0202-01")
# uuid=barcode2uuid(barcode)
