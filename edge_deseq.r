library("pasilla")
library("DESeq2")
library("edgeR")

fn<-c("B2B3gene")
alldata <- read.csv(paste(fn,"_count_matrix.csv",sep="") ,row.names=1)
coldata <- read.table("sample-group.txt", row.names=1,header = T)
cts<-alldata[,c(1:46)]

all(rownames(coldata) %in% colnames(cts))
dds <- DESeqDataSetFromMatrix(countData = cts,
                              colData = coldata,
                              design = ~ group + condition)
dds <- DESeq(dds)
degs<-resultsNames(dds) # lists the coefficients

normalized_counts <- counts(dds, normalized=TRUE)
###????鍑虹幇浠€涔堥棶棰樹簡锛燂紵锛燂紵
normalized_counts_mad <- apply(normalized_counts, 1, mad)
#normalized_counts <- normalized_counts[order(normalized_counts_mad, decreasing=T), ]


rld <- rlog(dds, blind=FALSE)
rlogMat <- assay(rld)
rlogMat <- rlogMat[order(normalized_counts_mad, decreasing=T), ]
hmcol <- colorRampPalette(brewer.pal(9, "GnBu"))(100)
pearson_cor <- as.matrix(cor(rlogMat, method="pearson"))
hc <- hcluster(t(rlogMat), method="pearson")
pdf(file=paste(paste(fn,"_cluster.pdf",sep="")), pointsize=10)
heatmap.2(pearson_cor, Rowv=as.dendrogram(hc), symm=T, trace="none",
          col=hmcol, margins=c(11,11), main="The pearson correlation of each
            sample")
dev.off()

pca_data <- plotPCA(rld, intgroup=c("group"), returnData=T, ntop=50000)
pdf(file=paste(paste(fn,"_pca.pdf",sep="")), pointsize=10)
plot(pca_data[,1:2],pch=19)
text(pca_data[,1],pca_data[,2]+1,row.names(pca_data),cex=0.5 )
dev.off()


design <- model.matrix(~ coldata$group+coldata$condition)

y <- estimateDisp(cts, design, robust=TRUE)

y <- DGEList(counts=cts,group=coldata$group)

y <- calcNormFactors(y)
#y <- estimateGLMCommonDisp(y,design,verbose=TRUE)
pdf(paste(fn,"-MDS.PDF",sep=""))
plotMDS(y)# show the outlier
dev.off()
