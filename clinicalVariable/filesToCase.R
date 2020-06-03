library("rjson")
data <- fromJSON(file="HNSCnorm.json")

output <- data.frame()

files <- list()
for ( i in 1:length(data)) {
  files <- c(files, substr(data[[i]]$file_name,1,49))
}

case_id <- list()
for (i in 1:length(data)) {
  case_id <- c(case_id, data[[i]]$cases[[1]]$case_id)
}

output <- cbind(files, case_id)
colnames(output) = c("file_name", "case_uuid")

write.table(output, file="HNSCnorm_case_uuid", sep = "\t", quote=FALSE, row.names = FALSE, col.names = TRUE)

