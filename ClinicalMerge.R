#######################################################################################
# merge all duplicated columns (clinical variables) ###################################

raw <- read.csv("TCGA-LUAD_clinical_from_XML.csv", colClasses = "character", na.strings = c("NA",""), header = FALSE)
colnames(raw) <- t(raw[1,])
raw <- raw[-1,]

mat <- matrix(nrow = nrow(raw))
newMat <- 1

matNames <- c()

count <- 1

while(1) {
  thisCol <- colnames(raw)[count]
  len <- nchar(thisCol)
  if(substr(thisCol,len-1,len) %in% c(".1",".2")) {
    thisCol <- substr(thisCol,1,len-2)
  }

  len <- nchar(thisCol)
  if(substr(thisCol,len-1,len) %in% c(".x",".y")) {
    thisCol <- substr(thisCol,1,len-2)
  }


  if(thisCol %in% matNames) {
    j <- which(matNames == thisCol)
    for (i in 1:nrow(raw)) {
      if(!is.na(raw[i,count])) { # get rows with values
        mat[i,j] <- raw[i,count]
      }
    }
  } else { # thisCol not in matrix yet
    matNames <- c(matNames, thisCol) 
    if(newMat == 1) {
      newMat <- 0
      mat <- raw[,count]
    } else {
      mat <- cbind(mat,raw[,count])
    }
  }
  
  count <- count+1
  if(count > ncol(raw)){
    break
  }
}

colnames(mat) <- matNames
write.table(mat, file="LUAD_clinical_merged_variable.csv", sep=",", quote=FALSE, row.names = FALSE, col.names = TRUE, na="")


#######################################################################################
# merge all duplicated rows (cases/patients) ##########################################

count <- nrow(mat)
numVariables <- ncol(mat)
uuid <- which(matNames == "bcr_patient_uuid")

while(1) {
  thisCase <- as.character(mat[count,uuid])
  if(thisCase == as.character(mat[count-1,uuid])) {
    for (j in 1:numVariables) {
      if(!(is.na(mat[count,j]))) { # get cols with values
        mat[count-1,j] <- mat[count,j]
      }
    }
    mat <- mat[-count,]
  }

  count <- count-1
  if(count <= 1){
    break
  }
}

write.table(mat, file="LUAD_clinical_merged_variable+case.csv", sep=",", quote=FALSE, row.names = FALSE, col.names = TRUE, na="")


#######################################################################################
# split stage tnm categories ##########################################################
library(stringr)

tnm <- which(matNames == "stage_event_tnm_categories")
mat <- cbind(mat,mat[,tnm]) # col for T
colT <- ncol(mat)
mat <- cbind(mat,mat[,tnm]) # col for N
colN <- ncol(mat)
mat <- cbind(mat,mat[,tnm]) # col for M
colM <- ncol(mat)
matNames <- c(matNames, "pathologic_T", "pathologic_N", "pathologic_M")
for (i in 1:nrow(mat)) {
  indT <- str_locate(mat[i,tnm], "T")[1,1]
  indN <- str_locate(mat[i,tnm], "N")[1,1]
  indM <- str_locate(mat[i,tnm], "M")[1,1]
  if (!is.na(indT) && !is.na(indN) && !is.na(indM)) {
    mat[i,colT] <- substr(mat[i,tnm],indT,indN-1)
    mat[i,colN] <- substr(mat[i,tnm],indN,indM-1)
    mat[i,colM] <- substr(mat[i,tnm],indM,nchar(mat[i,tnm]))
  }
  if (!is.na(indT) && !is.na(indN) && is.na(indM)) { # M DNE
    mat[i,colT] <- substr(mat[i,tnm],indT,indN-1)
    mat[i,colN] <- substr(mat[i,tnm],indN,nchar(mat[i,tnm]))
    mat[i,colM] <- ""
  }
  if (!is.na(indT) && is.na(indN) && !is.na(indM)) { # N DNE
    mat[i,colT] <- substr(mat[i,tnm],indT,indM-1)
    mat[i,colN] <- ""
    mat[i,colM] <- substr(mat[i,tnm],indM,nchar(mat[i,tnm]))
  }
  if (!is.na(indT) && is.na(indN) && is.na(indM)) { # T DNE
    mat[i,colT] <- ""
    mat[i,colN] <- substr(mat[i,tnm],indN,indM-1)
    mat[i,colM] <- substr(mat[i,tnm],indM,nchar(mat[i,tnm]))
  }
  if (is.na(indT) && is.na(indN) && !is.na(indM)) { # T DNE
    mat[i,colT] <- ""
    mat[i,colN] <- ""
    mat[i,colM] <- substr(mat[i,tnm],indM,nchar(mat[i,tnm]))
  }
  if (is.na(indT) && !is.na(indN) && is.na(indM)) { # T DNE
    mat[i,colT] <- ""
    mat[i,colN] <- substr(mat[i,tnm],indN,nchar(mat[i,tnm]))
    mat[i,colM] <- ""
  }
  if (!is.na(indT) && is.na(indN) && is.na(indM)) { # T DNE
    mat[i,colT] <- substr(mat[i,tnm],indT,nchar(mat[i,tnm]))
    mat[i,colN] <- ""
    mat[i,colM] <- ""
  }
  if (is.na(indT) && is.na(indN) && is.na(indM)) {
    mat[i,colT] <- ""
    mat[i,colN] <- ""
    mat[i,colM] <- ""
  }
}

matNames <- matNames[-tnm]
mat <- mat[,-tnm]

colnames(mat) <- matNames
write.table(mat, file="LUAD_clinical_tnm_split.csv", sep=",", quote=FALSE, row.names = FALSE, col.names = TRUE, na="")
