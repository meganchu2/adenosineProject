## ----setup, include=FALSE------------------------------------------------
knitr::opts_chunk$set(echo = TRUE)
knitr::opts_knit$set(progress = FALSE)

## ----message=FALSE, warning=FALS/E, include=FALSE-------------------------
library(TCGAbiolinks)
library(SummarizedExperiment)
library(dplyr)
library(DT)
library(data.table)
library(regexPipes)

  getclinical <- function(proj){
      message(proj)
      while(1){
          result = tryCatch({
              query <- GDCquery(project = proj, data.category = "Clinical",file.type = "xml")
              GDCdownload(query)
              clinical <- GDCprepare_clinic(query, clinical.info = "patient")
              for(i in c("admin","radiation","follow_up","drug","new_tumor_event")){
                  message(i)
                  aux <- GDCprepare_clinic(query, clinical.info = i)
                  if(is.null(aux) || nrow(aux) == 0) next
                  # add suffix manually if it already exists
                  replicated <- which(grep("bcr_patient_barcode",colnames(aux), value = T,invert = T) %in% colnames(clinical))
                  colnames(aux)[replicated] <- paste0(colnames(aux)[replicated],".",i)
                  if(!is.null(aux)) clinical <- merge(clinical,aux,by = "bcr_patient_barcode", all = TRUE)
              }
              readr::write_csv(clinical,path = paste0(proj,"_clinical_from_XML.csv")) # Save the clinical data into a csv file
              return(clinical)
          }, error = function(e) {
              message(paste0("Error clinical: ", proj))
          })
      }
  }