targets <- read.delim("Targets.txt", header = TRUE)
cpm<-read.delim("HNSC_CPM", header = FALSE, colClasses = "character")
corr <- read.delim("adenosine correlated genes.txt", header = TRUE)
colNames <- corr[,1]
counter <- 0

mat <- data.frame(row.names = targets$files)

for (i in 1:nrow(corr)) {
  x <- which(substr(cpm[,1],1,15) == corr[i,2])
  if (length(x) > 0) {
    temp <- t(cpm[x[1],])
    temp <- temp[2:length(temp)]
    mat <- cbind(mat,temp)
  } else {
    colNames <- colNames[-(i-counter)]
    counter <- (counter + 1)
  }
}

colnames(mat) <- colNames
write.table(mat, file="SpearmanInput", sep="\t", quote=FALSE, row.names = TRUE, col.names = TRUE)
