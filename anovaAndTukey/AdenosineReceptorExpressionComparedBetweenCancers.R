#This takes around 5 min to run per gene

# you can change the list of genes you want to look at
genes = c("ENSG00000163485",#"ADORA1"
			"ENSG00000128271", #"ADORA2A"
			"ENSG00000170425") #"ADORA2B"
geneName = c("ADORA1","ADORA2A","ADORA2B") # make sure this matches with order in genes variable


cancers <- c("BRCA","COAD","HNSC","LUAD","LUSC","PRAD") 

for(i in 1:length(genes)){
mat <- data.frame()
for(c in cancers){
	#cat(c,end="\n")
	targets <- read.delim(paste("Targets",c,".txt",sep=""),header=TRUE)
	#cat(paste("Total Samples: ",nrow(targets)),end="\n")
	cancerSample <- which(targets$group == "Cancer") # get which expression values are cancer ones
	#cat(paste("Num Cancer Samp: ", length(cancerSample)),end="\n")
	
	cpm <- read.delim(paste(c,"_CPM",sep=""), header = FALSE, colClasses = "character")
	x <- which(substr(cpm[,1],1,15) == genes[i]) # get expression values for this gene
	temp <- t(cpm[x[1],]) # get only expression values of this gene
	expression0 <- temp[2:length(temp)]
	#cat(paste("Num of expression data: ",length(expression0)),end="\n")
	expression <- expression0[cancerSample] # get only cancer expression values
	#cat(paste("Num of cancer expression data: ",length(expression)),end="\n\n")
	cancer <- rep(c, times=length(expression))
	toadd <- cbind(expression,cancer)
	if(ncol(mat) == 0 & nrow(mat) == 0){ mat <- toadd }
	else { mat <- rbind(mat,toadd) }
	
} 
mat2 <- mat
mat <-as.data.frame(mat)
mat$expression <- as.numeric(mat$expression)

pdf(paste(geneName[i]," expression vs. cancer .pdf"));boxplot(mat$expression~mat$cancer,main=paste(geneName[i]," expression vs. cancer"),xlab="cancer", ylab = paste(geneName[i]," Expression"),pch=16);dev.off()

#perform anova test
res.aov <- aov(expression ~ cancer, data = mat)
cat(paste("ANOVA for ",geneName[i]),end="\n")
print(summary(res.aov))
cat("\n")

#perform Tukey test
cat(paste("TUKEY TEST for ", geneName[i]),end="\n")
print(TukeyHSD(res.aov))
cat("\n")

}



