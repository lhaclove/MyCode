setwd("D:\\New folder1\\404")
library(edgeR)
#read count data

filename=("404")
alldata <- read.csv("D:\\New folder1\\all.csv",row.names=1)
#alldata <- read.delim("allcountB2.txt",row.names=1)
rnaidata<- read.csv("307V4count.csv",row.names=1)

x<-alldata[c(1:12)]
#define the data and group
#targets<-data.frame(TRAN=x[1:3],NT=x[4:6])


#defane the grorup and design
group<-factor(c(1,1,1,2,2,2,3,3,3,3,3,3))

y <- DGEList(counts=x,group=group)

y <- calcNormFactors(y)
#y <- estimateGLMCommonDisp(y,design,verbose=TRUE)
pdf(paste(filename,"-MDS.PDF",sep=""))
plotMDS(y)# show the outlier
dev.off()
y <- estimateCommonDisp(y,verbose=TRUE)
y <- estimateTagwiseDisp(y)
#plotBCV(y);
et <- exactTest(y)
de<-topTags(et,n=50000)



#order pvalues and fdr with same order of CPM
pval<-de$table[order(rownames(de$table)),]
pval<-pval[,-2]
Gene<-NULL
Gene<-rownames(pval)
da<-data.frame(Gene,cpm(x))
db<-data.frame(Gene,pval)
allexp<-merge(da,db,by="Gene")
write.csv(allexp, file = paste(filename,"-AllExP.csv",sep=""), row.names = F)
# Summary the DEGs, default filtered by FDR=0.05 and draw valcano map
detags<-rownames(de$table)
summary(de<-decideTestsDGE(et)) 
head(de)



detags <- rownames(y)[as.logical(de)]
pdf(paste(filename,"-DEG.PDF",sep=""))
plotSmear(et, de.tags=detags);
abline(h=c(-1, 1), col="blue");
dev.off()
#papre DEGs data
degs<-cpm(x)[as.logical(de),]
degpval<-pval[as.logical(de),]
keep <- rowSums(degs>1) >= 3& ((degpval[,1]>=1 )| (degpval[,1]<=-1))
degs<-degs[keep,]
degpval<-degpval[keep,]
dim(degs)
Gene<-rownames(degs)
da<-data.frame(Gene,degs)
db<-data.frame(Gene,degpval)
alldegs<-merge(da,db,by="Gene")
write.csv(alldegs, file = paste(filename,"-DEG.csv",sep=""), row.names = F)

###PCA
library(ggplot2)
library(ggbiplot)

res <- prcomp(t(x), center = TRUE, scale = FALSE)
names(res)
screeplot(res,type="lines")
temp<-predict(res) 
plot(temp[,1:2])
screeplot(temp,type="lines")

ggbiplot(res, obs.scale = 1, var.scale = 1,
         groups = group, ellipse = TRUE,var.axes = F,labels = names(alldata))

###PCA   PLOT

x_t <- t(alldata)

x <- as.data.frame(x_t)
x$group <- group
pca <- prcomp(t(rnaidata), scale=F)
percentVar <- pca$sdev^2 / sum( pca$sdev^2)
autoplot(pca, x, colour="group") + xlab(paste0("PC1 (", round(percentVar[1]*100), "% variance)")) + ylab(paste0("PC2 (", round(percentVar[2]*100), "% variance)")) + theme_bw() + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + theme(legend.position="right")
