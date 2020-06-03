# you can change the list of genes you want to look at
genes = c("ADORA1","ADORA2A","ADORA2B") # ,"AC011462.1","ENTPD1","ESR1","GP5","IGHM","IL4R","NFKB1")


# you can change group1 and group2 cancers to the ones that you want to compare
group1 <- c("BRCA","COAD","HNSC","LUSC")  
group2 <- c("PRAD") 


# you can change these cutoffs
sig <- .25 # can make this value smaller if you want smaller FDR cutoff or more significance (.25 is standard)
insig <- .25 # can make this value bigger if you want higher FDR curoff for more insignificance


print(paste("group1: ",paste(group1,collapse=", ")))
print(paste("group2: ",paste(group2,collapse=", ")))


samePos = c()
sameNeg = c()	

for(j in 1:length(genes)){#genes loop
print(genes[j])

g1upenriched <- c()
g1dnenriched <- c()
g2upenriched <- c()
g2dnenriched <- c()

for (i in 1:length(group1)){#cancers loop 1
#print(cancers[i])
pos <- read.delim(paste("gsea_report_for_",genes[j],"_pos_",group1[i],".xls",sep=""))
neg <- read.delim(paste("gsea_report_for_",genes[j],"_neg_",group1[i],".xls",sep=""))

subPos <- pos[pos$FDR.q.val < sig,] 
subNeg <- neg[neg$FDR.q.val < sig,]

sortPos <- subPos[order(subPos$FDR.q.val),]
sortNeg <- subNeg[order(subNeg$FDR.q.val),]


if (length(samePos)==0 & length(sameNeg)==0 & i==1){g1upenriched = as.vector(sortPos$NAME);g1dnenriched=as.vector(sortNeg$NAME)}
else {
  both <- intersect(sortPos$NAME,g1upenriched)
  g1upenriched <- g1upenriched[g1upenriched %in% both]
  both <- intersect(sortNeg$NAME,g1dnenriched)
  g1dnenriched <- g1dnenriched[g1dnenriched %in% both]
}

}# end cancers loop 1
#print(paste(length(g1upenriched)," ",length(g1dnenriched)))


for (i in 1:length(group2)){#check if unenriched in group 2
#print(cancers[i])
pos <- read.delim(paste("gsea_report_for_",genes[j],"_pos_",group2[i],".xls",sep=""))
neg <- read.delim(paste("gsea_report_for_",genes[j],"_neg_",group2[i],".xls",sep=""))

subPos <- pos[pos$FDR.q.val < insig,] 
subNeg <- neg[neg$FDR.q.val < insig,]

sortPos <- subPos[order(subPos$FDR.q.val),]
sortNeg <- subNeg[order(subNeg$FDR.q.val),]


both <- intersect(sortPos$NAME,g1upenriched)
g1upenriched <- g1upenriched[!(g1upenriched %in% both)]

both <- intersect(sortNeg$NAME,g1upenriched)
g1upenriched <- g1upenriched[!(g1upenriched %in% both)]

both <- intersect(sortPos$NAME,g1dnenriched)
g1dnenriched <- g1dnenriched[!(g1dnenriched %in% both)]

both <- intersect(sortNeg$NAME,g1dnenriched)
g1dnenriched <- g1dnenriched[!(g1dnenriched %in% both)]
#print(paste(length(g1upenriched)," ",length(g1dnenriched)))

}# end of check
print(paste("g1upenriched: ", length(g1upenriched)))
print(g1upenriched)
print(paste("g1dnenriched: ",length(g1dnenriched)))
print(g1dnenriched)



for (i in 1:length(group2)){#cancers loop 2
#print(cancers[i])
pos <- read.delim(paste("gsea_report_for_",genes[j],"_pos_",group2[i],".xls",sep=""))
neg <- read.delim(paste("gsea_report_for_",genes[j],"_neg_",group2[i],".xls",sep=""))

subPos <- pos[pos$FDR.q.val < sig,] 
subNeg <- neg[neg$FDR.q.val < sig,]

sortPos <- subPos[order(subPos$FDR.q.val),]
sortNeg <- subNeg[order(subNeg$FDR.q.val),]


if (length(samePos)==0 & length(sameNeg)==0 & i==1){g2upenriched = as.vector(sortPos$NAME);g2dnenriched=as.vector(sortNeg$NAME)}
else {
  both <- intersect(sortPos$NAME,g2upenriched)
  g2upenriched <- g2upenriched[g2upenriched %in% both]
  both <- intersect(sortNeg$NAME,g2dnenriched)
  g2dnenriched <- g2dnenriched[g2dnenriched %in% both]
}

}# end cancers loop 2


#print("group2")
#print(paste(length(g2upenriched)," ",length(g2dnenriched)))
for (i in 1:length(group1)){#check if unenriched in group 1
#print(cancers[i])
pos <- read.delim(paste("gsea_report_for_",genes[j],"_pos_",group1[i],".xls",sep=""))
neg <- read.delim(paste("gsea_report_for_",genes[j],"_neg_",group1[i],".xls",sep=""))

subPos <- pos[pos$FDR.q.val < insig,] 
subNeg <- neg[neg$FDR.q.val < insig,]

sortPos <- subPos[order(subPos$FDR.q.val),]
sortNeg <- subNeg[order(subNeg$FDR.q.val),]


both <- intersect(sortPos$NAME,g2upenriched)
g2upenriched <- g2upenriched[!(g2upenriched %in% both)]

both <- intersect(sortNeg$NAME,g2upenriched)
g2upenriched <- g2upenriched[!(g2upenriched %in% both)]

both <- intersect(sortPos$NAME,g2dnenriched)
g2dnenriched <- g2dnenriched[!(g2dnenriched %in% both)]

both <- intersect(sortNeg$NAME,g2dnenriched)
g2dnenriched <- g2dnenriched[!(g2dnenriched %in% both)]
#print(paste(length(g2upenriched)," ",length(g2dnenriched)))
}# end of check
print(paste("g2upenriched: ", length(g2upenriched)))
print(g2upenriched)
print(paste("g2dnenriched: ",length(g2dnenriched)))
print(g2dnenriched)

}# end of genes loop




