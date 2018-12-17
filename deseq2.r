library("pasilla")
library("DESeq2")
setwd("d:\\RNA3\\A7L")
fn<-c("gene")

#for (fn in fileNamebase){
alldata <- read.csv(paste(fn,"_count_matrix.csv",sep="") ,row.names=1)
coldata <- read.table("sample-group.txt", row.names=1,header = T)
cts<-alldata[,c(1:34)]

all(rownames(coldata) %in% colnames(cts))
dds <- DESeqDataSetFromMatrix(countData = cts,
                              colData = coldata,
                              design = ~ group)
dds <- DESeq(dds)
degs<-resultsNames(dds) # lists the coefficients

for(i in 2:length(degs)){
  print(degs[i])
  res <- results(dds, name=degs[i])
# or to shrink log fold changes association with condition:
#res <- lfcShrink(dds, coef=degs[i], type="apeglm")
plotMA(res, ylim=c(-2,2))
write.csv(as.data.frame(res), 
          file=paste(degs[i],paste(fn,"_allexp.csv",sep="")))
diff_gene_deseq2 <- subset(res,padj < 0.05 & (log2FoldChange >1 | log2FoldChange < -1)) 
write.csv(as.data.frame(diff_gene_deseq2), 
          file=paste(degs[i],paste(fn,"_deffexp.csv",sep="")))
}


}
