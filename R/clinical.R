accumulate.df=function(df, rows, stringsAsFactors = default.stringsAsFactors()){
  #Adds new rows to an existing dataframe "df"
  #The function is similar to "rbind", however
  #the difference is that columns in rows that
  #are not present in df are padded onto df
  #with all NA values before adding the rows.
  #Also the columns of rows are re-arranged
  #to match the column order in df (with possible
  #NAs for columns present in df but not in rows)
  if(is.null(df))
    return(rows)
  
  if(is.null(rows))
    return(df)
  
  #Match the new rows to the accumulated data
  q1=match(colnames(df), colnames(rows))
  qw=which(is.na(q1))
  if(length(qw)>0){
    #Add all-NA columns to match up with df
    temp.cols=df[1,qw, drop=F]
    for(j in 1:ncol(temp.cols))
      temp.cols[,j]=NA
    rownms=rownames(rows)#Save the rownames to apply to the cbind-ed result
    rownames(rows)=NULL
    rownames(temp.cols)=NULL
    rows=cbind(rows, temp.cols, stringsAsFactors=stringsAsFactors)
    rownames(rows)=rownms
    
    q1=match(colnames(df), colnames(rows))
  }
  df=rbind(df, rows[,q1])
  
  #Add any new data fields from the new rows to the accumulated data
  q2=which(!(colnames(rows) %in% colnames(df)))
  if(length(q2)>0){
    #Pad the patient data with the new columns
    rownms=rownames(df)#Save the rownames to apply to the cbind-ed result
    rownames(df)=NULL
    rownames(rows)=NULL
    df=cbind(df, rows[rep(1, nrow(df)), q2, drop=F], stringsAsFactors=stringsAsFactors)
    rownames(df)=rownms
    #Set all but the newly added rows as missing
    df[1:(nrow(df)-nrow(rows)), (ncol(df)-length(q2)+1):ncol(df)]=NA
  }
  return(df)
}

read.clinical.xml=function(file){
  doc = xmlInternalTreeParse(file)
  
  #Get the namespace name
  ns=names(xmlNamespace(xmlRoot(doc)))
  
  patient=getNodeSet(doc, path=paste("/", ns, ":tcga_bcr/", ns,":patient", sep=""))
  
  row=xmlToDataFrame(nodes=patient, stringsAsFactors = F)
  
  patient.fields=xmlChildren(patient[[1]])
  #For any fields that have subfields, concatenate the subfields
  for(j in 1:length(patient.fields)){
    children=xmlChildren(patient.fields[[j]])
    if(!any(names(children)=="text") & length(children)>0){
      to.df=FALSE
      for(k in 1:length(children)){
        if(length(xmlChildren(children[[k]]))>0){
          to.df=TRUE
          break
        }
      }
      if(!to.df)
        next
      sub.data=xmlToDataFrame(children, stringsAsFactors = F)
      n=rowSums(!is.na(sub.data))
      sub.data=sub.data[which(n>0),1:ncol(sub.data), drop=F]
      if(ncol(sub.data)==1 & nrow(sub.data)==length(children)){  
        sub.data=as.data.frame(t(sub.data))
        colnames(sub.data)=names(children)
      }
      dat.array=c()
      for(k in 1:nrow(sub.data))
        dat.array[k]=paste(paste(colnames(sub.data), unlist(sub.data[k,]), sep="="), collapse=";")
      
      row[,j]=paste(dat.array, collapse="|")
    }
    rm(children)
  }
  #       patient.barcode=xmlValue(patient.fields[["bcr_patient_barcode"]])
  #       row=data.frame(patient.barcode=patient.barcode)
  #       row$patient.uuid=xmlValue(patient.fields[["bcr_patient_uuid"]])
  #       row$patient.id=xmlValue(patient.fields[["patient_id"]])
  #       row$tumor.tissue.site=xmlValue(patient.fields[["tumor_tissue_site"]])
  #       row$histological.type=xmlValue(patient.fields[["histological_type"]])
  #       row$gender=xmlValue(patient.fields[["gender"]])
  #       row$vital.status=xmlValue(patient.fields[["vital_status"]])
  #       row$days.to.birth=xmlValue(patient.fields[["days_to_birth"]])
  #       row$days.to.death=xmlValue(patient.fields[["days_to_death"]])
  #       row$days.to.last.followup=xmlValue(patient.fields[["days_to_last_followup"]])
  #       row$race=xmlValue(patient.fields[["race"]])
  #       row$tissue.source.site=xmlValue(patient.fields[["tissue_source_site"]])
  #       row$history.of.neoadjuvant.treatment=xmlValue(patient.fields[["history_of_neoadjuvant_treatment"]])
  #       row$icd.o.3.site=xmlValue(patient.fields[["icd_o_3_site"]])
  #       row$icd.o.3.histology=xmlValue(patient.fields[["icd_o_3_histology"]])
  #       row$icd.10=xmlValue(patient.fields[["icd_10"]])
  #       row$days.to.initial.pathologic.diagnosis=xmlValue(patient.fields[["days_to_initial_pathologic_diagnosis"]])
  #       row$age.at.initial.pathologic.diagnosis=xmlValue(patient.fields[["age_at_initial_pathologic_diagnosis"]])
  #       row$year.of.initial.pathologic.diagnosis=xmlValue(patient.fields[["year_of_initial_pathologic_diagnosis"]])
  #       row$person.neoplasm.cancer.status=xmlValue(patient.fields[["person_neoplasm_cancer_status"]])
  #       row$ethnicity=xmlValue(patient.fields[["ethnicity"]])
  #       row$karnofsky.performance.score=xmlValue(patient.fields[["karnofsky_performance_score"]])
  #       row$eastern.cancer.oncology.group=xmlValue(patient.fields[["eastern_cancer_oncology_group"]])
  #       row$performance.status.scale.timing=xmlValue(patient.fields[["performance_status_scale_timing"]])
  # 
  #       row$stage.event=xmlValue((xmlChildren(patient.fields[["stage_event"]]))[["tnm_categories"]])
  #       
  #       row$neoplasm.histologic.grade=xmlValue(patient.fields[["neoplasm_histologic_grade"]])
  #       row$residual.tumor=xmlValue(patient.fields[["residual_tumor"]])
  #       row$tumor.residual.disease=xmlValue(patient.fields[["tumor_residual_disease"]])
  #       row$jewish.origin=xmlValue(patient.fields[["jewish_origin"]])
  #       row$anatomic.neoplasm.subdivision=xmlValue(patient.fields[["anatomic_neoplasm_subdivision"]])
  #       row$initial.pathologic.diagnosis.method=xmlValue(patient.fields[["initial_pathologic_diagnosis_method"]])
  #       row$venous.invasion=xmlValue(patient.fields[["venous_invasion"]])
  #       row$lymphatic.invasion=xmlValue(patient.fields[["lymphatic_invasion"]])
  #       
  #       drugs=xmlChildren(patient.fields[["drugs"]])
  #       row$drug=""
  #       if(length(drugs)>0){
  #         full.drug.data=NULL
  #         for(j in 1:length(drugs)){
  #           drug.fields=xmlChildren(drugs[[j]])
  #           drug.data=c()
  #           for(k in 1:length(drug.fields))
  #             drug.data[k]=paste(names(drug.fields)[k], xmlValue(drug.fields[[k]]), sep="=")
  #           full.drug.data[j]=paste(drug.data, collapse=";")
  #           rm(drug.fields)
  #           gc()
  #         }
  #         row$drug=paste(full.drug.data, collapse="|") 
  #       }
  #       
  #       radiations=xmlChildren(patient.fields[["radiations"]])
  #       row$radiation=""
  #       if(length(radiations)>0){
  #         full.radiation.data=c()
  #         for(j in 1:length(radiations)){
  #           radiation.fields=xmlChildren(radiations[[j]])
  #           radiation.data=c()
  #           for(k in 1:length(radiation.fields))
  #             radiation.data[k]=paste(names(radiation.fields)[k], xmlValue(radiation.fields[[k]]), sep="=")
  #           full.radiation.data[j]=paste(radiation.data, collapse=";")
  #           rm(radiation.fields)
  #           gc()
  #         }
  #         row$radiation=paste(full.radiation.data, collapse="|") 
  #       }
  #       
  #       follow.ups=xmlChildren(patient.fields[["follow_ups"]])
  #       row$follow.up=""
  #       if(length(follow.ups)>0){
  #         full.follow.up.data=c()
  #         for(j in 1:length(follow.ups)){
  #           follow.up.fields=xmlChildren(follow.ups[[j]])
  #           follow.up.data=c()
  #           for(k in 1:length(follow.up.fields))
  #             follow.up.data[k]=paste(names(follow.up.fields)[k], xmlValue(follow.up.fields[[k]]), sep="=")
  #           full.follow.up.data[j]=paste(follow.up.data, collapse=";")
  #           rm(follow.up.fields)
  #           gc()
  #         }
  #         row$follow.up=paste(full.follow.up.data, collapse="|") 
  #       }
  
  rm(patient.fields)
  rm(patient)
  free(doc)
  rm(doc)
  return(row)
}

read.biospecimen.xml=function(file, simplify=T){
  #Read an xml file describing a biospecimen and return a data.frame
  #(if simplify=T) or a list of data frames (if simplify=F) of data
  #for the biospecimen, sample, portion, analyte and aliquot.
  #Only biospecimens from which aliquots have been made are returned
  require(XML)
  
  biospecimen.result=c()
  sample.result=c()
  portion.result=c()
  analyte.result=c()
  aliquot.result=c()
  
  doc = xmlInternalTreeParse(file)
  biospecimens=getNodeSet(doc, path="/bio:tcga_bcr/bio:patient")
  biospecimen.data=xmlToDataFrame(biospecimens, stringsAsFactors = F)
  for(i in 1:length(biospecimens)){
    #Split the samples node
    biospecimen.children=xmlChildren(biospecimens[[i]])
    if(length(biospecimen.children)==0)
      next
    if("samples" %in% names(biospecimen.children)){
      samples=xmlChildren(biospecimen.children[["samples"]])
      if(length(samples)==0)
        next
      sample.data=xmlToDataFrame(samples, stringsAsFactors = F)
      for(j in 1:length(samples)){
        sample.children=xmlChildren(samples[[j]])
        if(length(sample.children)==0)
          next
        if("portions" %in% names(sample.children)){
          portions=xmlChildren(sample.children[["portions"]])
          if(length(portions)==0)
            next
          portion.data=xmlToDataFrame(portions, stringsAsFactors = F)
          for(k in 1:length(portions)){
            portion.children=xmlChildren(portions[[k]])
            if(length(portion.children)==0)
              next
            if("analytes" %in% names(portion.children)){
              analytes=xmlChildren(portion.children[["analytes"]])
              if(length(analytes)==0)
                next
              #Not using xmlToDataFrame here since it gives warnings
              analyte.data=NULL
              for(l in 1:length(analytes))
                analyte.data=accumulate.df(analyte.data, xmlToDataFrame(analytes[l], stringsAsFactors = F), stringsAsFactors = F)

              for(l in 1:length(analytes)){
                analyte.children=xmlChildren(analytes[[l]])
                if(length(analyte.children)==0)
                  next
                if("aliquots" %in% names(analyte.children)){
                  aliquots=xmlChildren(analyte.children[["aliquots"]])
                  if(length(aliquots)==0)
                    next
                  aliquot.data=xmlToDataFrame(aliquots, stringsAsFactors = F)
                  n=nrow(aliquot.data)
                  biospecimen.result=accumulate.df(biospecimen.result, biospecimen.data[rep(i, n),])
                  sample.result=accumulate.df(sample.result, sample.data[rep(j, n),])
                  portion.result=accumulate.df(portion.result, portion.data[rep(k, n),])
                  analyte.result=accumulate.df(analyte.result, analyte.data[rep(l, n),])
                  aliquot.result=accumulate.df(aliquot.result, aliquot.data)
                  rm(aliquots)
                }
                rm(analyte.children)
              }
              rm(analytes)
            }
            rm(portion.children)
          }
          rm(portions)
        }
        rm(sample.children)
      }
      rm(samples)
    }
    rm(biospecimen.children)
  }
  rm(biospecimens)

  free(doc)
  rm(doc)
  
  qw=which(colnames(biospecimen.result)!="samples")
  biospecimen.result=biospecimen.result[,qw] 

  qw=which(colnames(sample.result)!="portions")
  sample.result=sample.result[,qw]
  
  qw=which(colnames(portion.result)!="analytes")
  portion.result=portion.result[,qw]
  
  qw=which(colnames(analyte.result)!="aliquots")
  analyte.result=analyte.result[,qw]
  
  if(simplify){
    #Combine the five data.frames together and simplify the output
    #by removing columns that are not commonly used
    qm=match(c("bcr_patient_barcode", "bcr_patient_uuid", "patient_id", "gender", "tissue_source_site"), colnames(biospecimen.result))
    biospecimen.result=biospecimen.result[,qm[!is.na(qm)]]
    qm=match(c("bcr_sample_barcode", "bcr_sample_uuid", "sample_type", "tumor_pathology", "is_ffpe"), colnames(sample.result))
    sample.result=sample.result[,qm[!is.na(qm)]]
    qm=match(c("bcr_portion_barcode", "bcr_portion_uuid", "portion_number"), colnames(portion.result))
    portion.result=portion.result[,qm[!is.na(qm)]]
    qm=match(c("bcr_analyte_barcode", "bcr_analyte_uuid", "analyte_type", "a260_a280_ratio", "protocols"), colnames(analyte.result))
    analyte.result=analyte.result[,qm[!is.na(qm)]]
    qm=match(c("bcr_aliquot_barcode", "bcr_aliquot_uuid"), colnames(aliquot.result))
    aliquot.result=aliquot.result[,qm[!is.na(qm)]]
    return(cbind(biospecimen.result, sample.result, portion.result, analyte.result, aliquot.result))
  }else{
    return(list(biospecimen=biospecimen.result, sample=sample.result, portion=portion.result, analyte=analyte.result,aliquot=aliquot.result))
  }
   
#   patient=getNodeSet(doc, path="/bio:tcga_bcr/bio:patient")
#   patient.fields=xmlChildren(patient[[1]])
#   bcr_patient_barcode=xmlValue(patient.fields[["bcr_patient_barcode"]])
#   patient.uuid=xmlValue(patient.fields[["bcr_patient_uuid"]])
#   patient.id=xmlValue(patient.fields[["patient_id"]])
#   tissue.source.site=xmlValue(patient.fields[["tissue_source_site"]])
#   
#   if("samples" %in% names(patient.fields)){
#     samples=xmlChildren(patient.fields[["samples"]])
#     if(length(samples)>0){
#       for(j in 1:length(samples)){
#         sample.fields=xmlChildren(samples[[j]])
#         sample.barcode=xmlValue(sample.fields[["bcr_sample_barcode"]])
#         sample.uuid=xmlValue(sample.fields[["bcr_sample_uuid"]])
#         sample.type=xmlValue(sample.fields[["sample_type"]])
#         tumor.pathology=xmlValue(sample.fields[["tumor_pathology"]])
#         if("portions" %in% names(sample.fields)){
#           portions=xmlChildren(sample.fields[["portions"]])
#           if(length(portions)>0){
#             for(k in 1:length(portions)){
#               portion.fields=xmlChildren(portions[[k]])
#               portion.barcode=xmlValue(portion.fields[["bcr_portion_barcode"]])
#               portion.uuid=xmlValue(portion.fields[["bcr_portion_uuid"]])
#               portion.number=xmlValue(portion.fields[["portion_number"]])
#               portion.is.ffpe=xmlValue(portion.fields[["is_ffpe"]])
#               if("analytes" %in% names(portion.fields)){
#                 analytes=xmlChildren(portion.fields[["analytes"]])
#                 if(length(analytes)>0){
#                   for(l in 1:length(analytes)){
#                     analyte.fields=xmlChildren(analytes[[l]])
#                     analyte.barcode=xmlValue(analyte.fields[["bcr_analyte_barcode"]])
#                     analyte.uuid=xmlValue(analyte.fields[["bcr_analyte_uuid"]])
#                     analyte.type=xmlValue(analyte.fields[["analyte_type"]])
#                     if("aliquots" %in% names(analyte.fields)){
#                       aliquots=xmlChildren(analyte.fields[["aliquots"]])
#                       if(length(aliquots)>0){
#                         for(m in 1:length(aliquots)){
#                           aliquot.fields=xmlChildren(aliquots[[m]])
#                           aliquot.barcode=xmlValue(aliquot.fields[["bcr_aliquot_barcode"]])
#                           aliquot.uuid=xmlValue(aliquot.fields[["bcr_aliquot_uuid"]])
#                           biospecimen.data=rbind(biospecimen.data, data.frame(bcr_patient_barcode=bcr_patient_barcode,
#                                                                               patient.uuid=patient.uuid,
#                                                                               patient.id=patient.id,
#                                                                               tissue.source.site=tissue.source.site,
#                                                                               sample.barcode=sample.barcode,
#                                                                               sample.uuid=sample.uuid,
#                                                                               sample.type=sample.type,
#                                                                               tumor.pathology=tumor.pathology,
#                                                                               portion.barcode=portion.barcode,
#                                                                               portion.uuid=portion.uuid,
#                                                                               portion.number=portion.number,
#                                                                               analyte.barcode=analyte.barcode,
#                                                                               analyte.uuid=analyte.uuid,
#                                                                               analyte.type=analyte.type,
#                                                                               aliquot.barcode=aliquot.barcode,
#                                                                               aliquot.uuid=aliquot.uuid,
#                                                                               stringsAsFactors=FALSE))
#                           
#                           rm(aliquot.fields)
#                         }
#                       }
#                       rm(aliquots)
#                     }
#                     rm(analyte.fields)
#                   }
#                 }
#                 rm(analytes)
#               }
#               rm(portion.fields)
#             }
#           }
#           rm(portions)
#         }
#         rm(sample.fields)
#       }
#     }
#     rm(samples)
#   }
#   rm(patient.fields)
#   rm(patient)
#   gc()
#   free(doc)
#   rm(doc)
#   gc()
#   return(biospecimen.data)
}

clinical.data.xml=function(clin.directory, simplify.biospecimen=T, verbose=F){
  require(XML)
  require(gsubfn)
  
  files=dir(clin.directory)
  biospecimen.files=files[grep(files, pattern="biospecimen")]
  if(length(biospecimen.files)>0){
    biospecimen.patients=strapply(biospecimen.files, pattern="TCGA-[[:alnum:]]+-[[:alnum:]]+")
    n=sapply(biospecimen.patients, length)
    if(all(n==1)){
      biospecimen.patients=unlist(biospecimen.patients)
    }else{
      warning("File without TCGA patient id")
      browser()
    }
  }
  
  clinical.files=files[grep(files, pattern="clinical")]
  if(length(clinical.files)>0){
    clinical.patients=strapply(clinical.files, pattern="TCGA-[[:alnum:]]+-[[:alnum:]]+")
    n=sapply(clinical.patients, length)
    if(all(n==1)){
      clinical.patients=unlist(clinical.patients)
    }else{
      warning("File without TCGA patient id")
      browser()
    }
  }
  #First load the biospecimen files
  if(simplify.biospecimen){
    biospecimen.data=NULL
  }else{
    biospecimen.data=list(biospecimen=c(), 
                          sample=c(),
                          portion=c(),
                          analyte=c(),
                          aliquot=c())
  }
  if(length(biospecimen.files)>0)
    for(i in 1:length(biospecimen.files)){
      bio.dat=read.biospecimen.xml(file.path(clin.directory, biospecimen.files[i]), simplify=simplify.biospecimen)
      if(simplify.biospecimen){
        biospecimen.data=accumulate.df(biospecimen.data, bio.dat)
      }else{
        for(field in c("biospecimen", "sample", "portion", "analyte", "aliquot"))
          biospecimen.data[[field]]=accumulate.df(biospecimen.data[[field]], bio.dat[[field]])
      }
      if(verbose)
        cat(i, "of", length(biospecimen.files), "biospecimen files read\n")
    }
  
  patient.data=NULL
  if(length(clinical.files)>0)
    for(i in 1:length(clinical.files)){
      row=read.clinical.xml(file.path(clin.directory, clinical.files[i]))
      patient.data=accumulate.df(patient.data, row)
      if(verbose)
        cat(i, "of", length(clinical.files), "clinical files read\n")
    }
  
  return(list(biospecimen=biospecimen.data, patient=patient.data))
}

#clin.directory="C:/Users/varmas/Documents/Work/Common data/TCGA data/Gene expression/TCGA download 2014-02-25/OV/Clinical/XML"
#res=clinica.data(clin.directory)

# biospecimen.fields=data.frame(
#   type=c("patient",
#          "patient",
#          "patient",
#          "patient",
#          "sample",
#          "sample",
#          "sample",
#          "sample",
#          "portion",
#          "portion",
#          "portion",
#          "portion",
#          "analyte",
#          "analyte",
#          "analyte"),
#   
#   field=c("bcr_patient_barcode",
#           "bcr_patient_uuid",
#           "patient_id",
#           "tissue_source_site",
#           "bcr_sample_barcode",
#           "bcr_sample_uuid",
#           "sample_type",
#           "tumor_pathology",
#           "bcr_portion_barcode",
#           "bcr_portion_uuid",
#           "portion_number",
#           "is_ffpe",
#           "bcr_analyte_barcode",
#           "bcr_analyte_uuid",
#           "analyte_type"),
#   
#   name=c("bcr_patient_barcode",
#          "patient.uuid",
#          "patient.id",
#          "tissue.source.site",
#          "sample.barcode",
#          "sample.uuid",
#          "sample.type",
#          "tumor.pathology",
#          "portion.barcode",
#          "portion.uuid",
#          "portion.number",
#          "portion.is.ffpe",
#          "analyte.barcode",
#          "analyte.uuid",
#          "analyte.type"),
#   stringsAsFactors=F)
# 
# read.file=function(filename){
#   d=read.delim(filename, check.names=F, stringsAsFactors=F, na=c("[Not Applicable]", "[Not Available]", "null", ""), header=F)
#   for(i in 1:10){
#     if(any(grepl(as.vector(d[i,]), pattern="^bcr"))){
#       header=d[i,]
#       break
#     }
#   }
#   for(i in 1:10){
#     if(any(grepl(as.vector(d[i,]), pattern="^TCGA"))){
#       data.row=i
#       break
#     }
#   }
#   d=d[data.row:nrow(d),]
#   colnames(d)=header
#   return(d)
#   
# }
# read.clinical=function(clin.directory, data.type=c("patient", "drug", "radiation", "follow_up")){
#   file.pattern=paste("clinical", data.type, sep="_")
#   file=dir(clin.directory, full.names=T, pattern=file.pattern)
#   if(length(file)==1){
#     d=read.file(file)
#    }else{
#     warning(paste("Could not find",file.pattern, "file"))
#     return(NULL)
#   }
# }
# patient.data=function(clin.directory){
#   #Read the patient information from various files and match them up together
#   patient.data=read.clinical(clin.directory, "patient")
#   
#   drug.data=read.clinical(clin.directory, "drug")
#   if(!is.null(patient.data) & !is.null(drug.data)){
#     qm=match.id(patient.data, drug.data)
#     patient.data=combine.data(patient.data, drug.data, qm, all.x=TRUE, all.y=TRUE)
#     patient.data=combine.same.cols(patient.data)
#   }else if(!is.null(drug.data)){
#     patient.data=drug.data
#   }
#   
#   radiation.data=read.clinical(clin.directory, "radiation")
#   if(!is.null(patient.data) & !is.null(radiation.data)){
#     qm=match.id(patient.data, radiation.data)
#     patient.data=combine.data(patient.data, radiation.data, qm, all.x=TRUE, all.y=TRUE)
#     patient.data=combine.same.cols(patient.data)
#   }else if(!is.null(radiation.data)){
#     patient.data=radiation.data
#   }
#   
#   followup.data=read.clinical(clin.directory, "follow_up")
#   if(!is.null(patient.data) & !is.null(followup.data)){
#     qm=match.id(patient.data, followup.data, use.any=TRUE)
#     patient.data=combine.data(patient.data, followup.data, qm, all.x=TRUE, all.y=TRUE)
#     patient.data=combine.same.cols(patient.data)
#   }else if(!is.null(followup.data)){
#     patient.data=followup.data
#   }
#   
#   return(patient.data)
# }
# 
# read.biospecimen=function(clin.directory, data.type){
#   file.pattern=paste("biospecimen", data.type, sep="_")
#   file=dir(clin.directory, full.names=T, pattern=file.pattern)
#   if(length(file)==1){
#     d=read.file(file)
#     return(d)
#   }else{
#     warning(paste("Could not find",file.pattern, "file"))
#     return(NULL)
#   }
# }
# 
# 
# biospecimen.data=function(clin.directory){
#   #Read the biospecimen information from various files and match them up together
#   
#   #####NOT FINISHED!!#############
#   browser()
#   normal.patient=read.biospecimen(clin.directory, "normal_control")
#   tumor.patient=read.biospecimen(clin.directory, "tumor_sample")
#   
#   sample.data=read.biospecimen(clin.directory, "sample")
#   qm=match(c("bcr_sample_barcode", "bcr_sample_uuid", "sample_type", "tumor_pathology"), colnames(sample.data))
#   qw=which(!is.na(qm))
#   sample.data=sample.data[,qm[qw]]
#   colnames(sample.data)=c("sample.barcode", "sample.uuid", "sample.type", "tumor.pathology")[qw]
#   
#   portion.data=read.biospecimen(clin.directory, "portion")
#   qm=match(c("bcr_sample_barcode", "bcr_portion_barcode", "bcr_portion_uuid", "portion_number", "is_ffpe"), colnames(portion.data))
#   qw=which(!is.na(qm))
#   portion.data=portion.data[,qm[qw]]
#   colnames(portion.data)=c("sample.barcode", "portion.barcode", "portion.uuid", "portion.number","portion.is.ffpe")[qw]
#   biospecimen.data=merge(sample.data, portion.data, by="sample.barcode")
#   
#   analyte.data=read.biospecimen(clin.directory, "analyte")
#   
#   analyte.barcode=xmlValue(analyte.fields[["bcr_analyte_barcode"]])
#   analyte.uuid=xmlValue(analyte.fields[["bcr_analyte_uuid"]])
#   analyte.type=xmlValue(analyte.fields[["analyte_type"]])
#   
#   if(!is.null(biospecimen.data) & !is.null(aliquot.data)){
#     qm=match.id(biospecimen.data, aliquot.data)
#     biospecimen.data=combine.data(biospecimen.data, aliquot.data, qm, all.x=TRUE, all.y=TRUE)
#     biospecimen.data=combine.same.cols(biospecimen.data)
#   }else if(!is.null(aliquot.data)){
#     biospecimen.data=aliquot.data
#   }
#   cqcf.data=read.biospecimen(clin.directory, "cqcf")
#   if(!is.null(biospecimen.data) & !is.null(cqcf.data)){
#     qm=match.id(biospecimen.data, cqcf.data)
#     biospecimen.data=combine.data(biospecimen.data, cqcf.data, qm, all.x=TRUE, all.y=TRUE)
#     biospecimen.data=combine.same.cols(biospecimen.data)
#   }else if(!is.null(cqcf.data)){
#     biospecimen.data=cqcf.data
#   }
#   
#   qm=match(c("bcr_patient_barcode",
#              "bcr_patient_uuid",
#              "patient_id",
#              "tissue_source_site",
#              "bcr_sample_barcode",
#              "bcr_sample_uuid",
#              "sample_type",
#              "tumor_pathology",
#              "bcr_portion_barcode",
#              "bcr_portion_uuid",
#              "portion_number",
#              "bcr_analyte_barcode",
#              "bcr_analyte_uuid",
#              "analyte_type"), colnames(biospecimen.data))
#   browser()
#   
#   biospecimen.data=biospecimen.data[,qm]
#   colnames(biospecimen.data)=c("bcr_patient_barcode",
#                                "patient.uuid",
#                                "patient.id",
#                                "tissue.source.site",
#                                "sample.barcode",
#                                "sample.uuid",
#                                "sample.type",
#                                "tumor.pathology",
#                                "portion.barcode",
#                                "portion.uuid",
#                                "portion.number",
#                                "analyte.barcode",
#                                "analyte.uuid",
#                                "analyte.type")
#   return(biospecimen.data)
# }
# 
# clinical.data.biotab=function(clin.directory){
#   patient.data=patient.data(clin.directory)
#   biospecimen.data=biospecimen.data(clin.directory)
# #   if(!is.null(clin.data)){
# #     qm=match.id(clin.data, bio.data)
# #     clin.data=combine.data(clin.data, bio.data, qm, all.x=TRUE, all.y=TRUE)
# #   }
#   return(list(biospecimen=biospecimen.data, patient=patient.data))
# }
