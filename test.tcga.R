#Test the TCGA scripts
library(TCGAquery)

tcga.download.dir="Test-dir"
if(!file.exists(tcga.download.dir))
  dir.create(tcga.download.dir)
#Set the data level
level=3

#Select the data type to download
data(TCGA.data.types)
data.type=TCGA.data.types[8,]

#Download data for all diseases
data(TCGA.diseases)
for(i in 5){#1:nrow(TCGA.diseases)){
  disease=TCGA.diseases[i,]
  disease.abbr=disease$Study.Abbreviation
  
  #Download and read assay data
  center=data.type$Short.Name
  platform=data.type$Platform.Alias
  disease.download.dir=file.path(tcga.download.dir, disease$Study.Abbreviation)
  if(!file.exists(disease.download.dir))
    dir.create(disease.download.dir)
  data.download.dir=file.path(disease.download.dir, "Data")
  if(!file.exists(data.download.dir))
    dir.create(data.download.dir)
  result=download.tcga(center=center, disease=disease$Study.Abbreviation, 
                       platform=platform, level=level,
                       download.dir=data.download.dir,
                       force.delete=TRUE, verbose=FALSE)
  
  #data.files=get.data(data.download.dir, data.download.dir, tcga.data.type="RSEM_genes_normalized")
  data.files=get.data(data.download.dir, data.download.dir, file.type="broad.mit.edu:segmented_scna_hg19:Genome_Wide_SNP_6:01")
  
  #Download and read the clinical data
  platformType="C"
  #Download the clinical data into a temporary directory
  clinical.download.dir=file.path(disease.download.dir, "Clinical")
  if(!file.exists(clinical.download.dir))
    dir.create(clinical.download.dir)
  download.tcga(disease=disease$Study.Abbreviation,
                platformType=platformType,
                download.dir=clinical.download.dir,
                force.delete=TRUE, verbose=FALSE)
  
  #Parse and load the clinical data
  clin.data=clinical.data.xml(clin.directory=file.path(clinical.download.dir, "Clinical", "XML"))
  
  save(clin.data, file=file.path(clinical.download.dir, paste(disease.abbr,"clinical.Rdata", sep=".")))
  
  #Match the clinical data to the methylation data
  load(file.path(data.download.dir, data.files[[1]]$out.file.name[1]))
  
  qm=match(toupper(expd[,"Extract Name"]), toupper(clin.data$biospecimen$aliquot.uuid))
  #   qw=which(!is.na(qm))
  #   expd=expd[qw,]
  #   X=X[,qw]
  #   biospecimen.data=clin.data$biospecimen[qm[qw],]
  biospecimen.data=clin.data$biospecimen[qm,]
  
  qm=match(biospecimen.data$patient.barcode, clin.data$patient$patient.barcode)
  #   qw=which(!is.na(qm))
  #   expd=expd[qw,]
  #   X=X[,qw]
  #   biospecimen.data=biospecimen.data[qw,]
  #   patient.data=clin.data$patient[qm[qw],]
  patient.data=clin.data$patient[qm,]
  save(expd, X, biospecimen.data, patient.data, annot,file=paste(disease.abbr,"2", "matched.Rdata", sep="."))
}