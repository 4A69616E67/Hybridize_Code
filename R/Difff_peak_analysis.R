library(DiffBind)
min_log2_fold_change <- log2(1.5)
max_pvalue <- 0.05
prefix <- "example.prefix"
raw_dba <- dba(sampleSheet = "SampleSheet.example.csv", minOverlap = 2)
dbObj <- dba.count(raw_dba, bUseSummarizeOverlaps = T)
dba.plotPCA(dbObj, attributes = DBA_CONDITION, label = DBA_ID)
plot(dbObj)

dbObj <- dba.contrast(dbObj, categories = DBA_CONDITION, minMembers = 2)
dbObj <- dba.analyze(dbObj, method = DBA_ALL_METHODS)
dba.report(dbObj)
dba.show(dbObj)

dba.plotVenn(dbObj, contrast = 1, method = DBA_ALL_METHODS)
# dba.plotVolcano(dbObj)
# dba.plotBox(dbObj)

comp1.edgeR <- as.data.frame(dba.report(dbObj, method = DBA_EDGER, contrast = 1, th = 1))

comp1.deseq <- as.data.frame(dba.report(dbObj, method = DBA_DESEQ2, contrast = 1, th = 1))

# distinguish significance
comp1.deseq$significant <- "normal"
comp1.deseq$significant[which(comp1.deseq$Fold > min_log2_fold_change & comp1.deseq$FDR < max_pvalue)] <- "up"
comp1.deseq$significant[which(comp1.deseq$Fold < -min_log2_fold_change & comp1.deseq$FDR < max_pvalue)] <- "down"

######### write the result
openxlsx::write.xlsx(list("DESeq"=comp1.deseq,"edgR"=comp1.edgeR), file = paste(prefix, "DiffBind.xlsx", sep = "."))

write.table(comp1.deseq, paste(prefix, "DiffBind.deseq.tsv", sep = "."), sep = "\t", row.names = F, quote = F)
write.table(comp1.deseq, paste(prefix, "DiffBind.deseq.tsv", sep = "."), sep = "\t", row.names = F, quote = F)