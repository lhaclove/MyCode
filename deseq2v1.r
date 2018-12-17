library("pasilla")
library("DESeq2")
library("pheatmap")
library("amap")
library("ggplot2")
library("BiocParallel")
register(MulticoreParam(4))
library("gplots")
library("RColorBrewer")
setwd("d:\\RNA3\\tmp")
fn<-c("gene")

#for (fn in fileNamebase){


##define file and file name
alldata <- read.csv(paste(fn,"_count_matrix.csv",sep="") ,row.names=1)

alldata <- read.csv("count-a7-a7l.csv" ,row.names=1)


#coldata <- read.table("sample-group.txt", row.names=1,header = T)

cts<-alldata[,c(5,7:9,13,23,27,28)]
##因为对照重复性不好，因此可能会因为对照的选择影响表达量，因此使用两组对照做分析？
###group1:ZJMA,ZJMD,2611Y4:19,22,26,
###group2:2611Y7,2611Y5,2503Y2:23,27,28
###group3:ZJMB.ZJMC,2503Y4,2503Y5: 改组可能不好
###A7OE:15,17,18
###A7RNAi:10,11 A7RNAi2:12,13

19,22,26,

19,22,24,25,20:21,23,26:28
group<-c(rep("test8_RNAi_group1",4),rep("CK",3))

{
aa<-data.frame(group)
rownames(aa)<-colnames(cts)

all(rownames(aa) %in% colnames(cts))
dds <- DESeqDataSetFromMatrix(countData = cts,
                              colData = aa,
                              design = ~ group)
dds <- DESeq(dds)
degs<-resultsNames(dds) # lists the coefficients
for(i in 2:length(degs)){
  print(degs[i])
  res <- results(dds, name=degs[i])
# or to shrink log fold changes association with condition:
#res <- lfcShrink(dds, coef=degs[i], type="apeglm")

  pdf(file=paste(degs[i],paste(fn,"_MA.pdf",sep="")), pointsize=10)
  plotMA(res, ylim=c(-2,2))
dev.off()
  


###normaliz counts
normalized_counts <- counts(dds, normalized=TRUE)

normalized_counts_mad <- apply(normalized_counts, 1, mad)
normalized_counts <- normalized_counts[order(normalized_counts_mad, decreasing=T), ]

{aa<-as.data.frame(res)
  Gene<-rownames(aa)
  da<-data.frame(Gene,aa)
  bb<-as.data.frame(normalized_counts)
  Gene<-rownames(bb)
  db<-data.frame(Gene,bb)
  aadd<-merge(db,da,by="Gene")}

write.csv(as.data.frame(aadd), 
          file=paste(degs[i],paste(fn,"_all_count_exp.csv",sep="")),row.names = FALSE)

diff_gene_deseq2 <- subset(aadd,padj < 0.05 & (log2FoldChange >1 | log2FoldChange < -1)) 
row.names(diff_gene_deseq2)<-diff_gene_deseq2[,1]
write.csv(as.data.frame(diff_gene_deseq2), file=paste(degs[i],paste(fn,"_deff_count_exp.csv",sep="")),row.names = FALSE)

rld <- rlog(dds, blind=FALSE)
rlogMat <- assay(rld)
rlogMat <- rlogMat[order(normalized_counts_mad, decreasing=T), ]
hmcol <- colorRampPalette(brewer.pal(9, "GnBu"))(100)
pearson_cor <- as.matrix(cor(rlogMat, method="pearson"))
hc <- hcluster(t(rlogMat), method="pearson")
pdf(file=paste(degs[i],paste(fn,"_cluster.pdf",sep="")), pointsize=10)
heatmap.2(pearson_cor, Rowv=as.dendrogram(hc), symm=T, trace="none",
          col=hmcol, margins=c(11,11), main="The pearson correlation of each
sample")
dev.off()

pca_data <- plotPCA(rld, intgroup=c("group"), returnData=T, ntop=50000)
pdf(file=paste(degs[i],paste(fn,"_pca.pdf",sep="")), pointsize=10)
plot(pca_data[,1:2],pch=19,col=c("red","red","red","blue","blue","blue"))
text(pca_data[,1],pca_data[,2]+1,row.names(pca_data),cex=1 )
dev.off()
}

}

