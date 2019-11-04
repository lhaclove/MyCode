library("pasilla")
library("DESeq2")
setwd("d:\\RNA3\\")
fn<-c("B2B3gene")

#for (fn in fileNamebase){
alldata <- read.csv(paste(fn,"_count_matrix.csv",sep="") ,row.names=1)
alldata <- read.csv("count-a7-a7l.csv" ,row.names=1)
coldata <- read.table("sample-group.txt", row.names=1,header = T)
cts<-alldata[,c(1:40)]

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
  
  pdf(file=paste(degs[i],paste(fn,"_MA.pdf",sep="")), pointsize=10)
  plotMA(res, ylim=c(-2,2))
  dev.off()
  
  
  
  ###normaliz counts
  normalized_counts <- counts(dds, normalized=TRUE)
  ###????出现什么问题了？？？？
  normalized_counts_mad <- apply(normalized_counts, 1, mad)
  #normalized_counts <- normalized_counts[order(normalized_counts_mad, decreasing=T), ]
  
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
  deg<-merge(diff_gene_deseq2,subset(funaa , funaa$Gene%in%diff_gene_deseq2[,1]),by="Gene",all=T)
  
  write.csv(as.data.frame(deg), file=paste(degs[i],paste(fn,"_deff_count_exp.csv",sep="")),row.names = FALSE)
  write.table(as.data.frame(diff_gene_deseq2[1]), file=paste(degs[i],paste(fn,"_DEG_ID.txt",sep="")),quote=F,row.names = F)
  
  
  
  
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
  plot(pca_data[,1:2],pch=19)
  text(pca_data[,1],pca_data[,2]+1,row.names(pca_data),cex=0.5 )
  dev.off()
}

}
}
