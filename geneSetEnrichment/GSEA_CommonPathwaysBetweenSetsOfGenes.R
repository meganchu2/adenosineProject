genes = c("ADORA1","ADORA2A","ADORA2B","AC011462.1","ENTPD1","ESR1","GP5","IGHM","IL4R","NFKB1")
cancers = c("BRCA","COAD","HNSC","LUAD","LUSC","PRAD")
samePos = c()
sameNeg = c()

# you can change this list of vectors to compare any set of genes you want, just make sure the
# numbers in your vector are listed lowest to highest
# Ex: if I want to compare "ENTPD1" and "ESR1" (which are number 5 and 6 in the 'genes' vector 
# from line 1) I want to put c(5,6) in the list below and not c(6,5).
for(k in list(c(1,2,3), #A1 vs. A2A vs A2B
				c(1,2), #A1 vs. A2A
				c(1,3), #A1 vs. A2B
				c(2,3), #A2A vs. A2B
				c(2,5), #A2A vs. ENTPD1
				c(2,4), #A2A vs. TGFB1(AC011462.1)
				c(3,6), #A2B vs. ESR1
				c(1,9), #A1 vs. IL4R
				c(2,8), #A2A vs. IGHM
				c(3,7), #A2B vs. GP5
				c(1,10) #A1 vs. NFKB1
				)){
samePos = c()
sameNeg = c()	
for(j in k){#genes loop
print(genes[j])
for (i in 1:6){#cancers loop 1
#print(cancers[i])
pos <- read.delim(paste("gsea_report_for_",genes[j],"_pos_",cancers[i],".xls",sep=""))
neg <- read.delim(paste("gsea_report_for_",genes[j],"_neg_",cancers[i],".xls",sep=""))

subPos <- pos[pos$FDR.q.val < .25,] # I used .25 for cancer and canonical pathways, but .15 for immune pathways 
subNeg <- neg[neg$FDR.q.val < .25,] # because there were too many pathways in common

sortPos <- subPos[order(subPos$FDR.q.val),]
sortNeg <- subNeg[order(subNeg$FDR.q.val),]

if (length(samePos)==0 & length(sameNeg)==0 & i==1 & j==k[1]){samePos = as.vector(sortPos$NAME);sameNeg=as.vector(sortNeg$NAME)}
else {
  both <- intersect(sortPos$NAME,samePos)
  samePos <- samePos[samePos %in% both]
  both <- intersect(sortNeg$NAME,sameNeg)
  sameNeg <- sameNeg[sameNeg %in% both]
}

}# end cancers loop 1
#print(paste(length(samePos),length(sameNeg))) # checkpoint to make sure code is working correctly
}
print(paste("(+): ",length(samePos)))
if(length(samePos)>0){ print(samePos) }
print(paste("(-): ",length(sameNeg)))
if(length(sameNeg)>0){ print(sameNeg) }
}

