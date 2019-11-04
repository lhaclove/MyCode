setwd("D:\\RNA3\\report\\expdata\\test1\\A7LCK3edgeR\\")
library(edgeR)
funaa<-read.csv("d:\\RNA3\\v4fun.csv" ,header = T)
names(funaa)<-c("Gene","desc")
#read count data

alldata <- read.csv("count-a7-a7l.csv" ,row.names=1)
#alldata <- read.csv("transcript_count_matrix.csv",row.names=1)
#alldata <- read.delim("allcountB2.txt",row.names=1)

filenamebase<-list("A7L-SRDX")
control<-list("CK1","CK2","CK3")


vp<-alldata[,c(1:4)]
srdx<-alldata[,c(5:9)]
rnai1<-alldata[,c(10,11)]
rnai2<-alldata[,c(12,13)]
a7loe<-alldata[,c(29:31)]
a7lrnai<-alldata[,c(32:34)]
oe<-alldata[,c(15,17,18)]
ck1<-alldata[,c(23,27,28)]
ck2<-alldata[,c(19,22,26)]
ck3<-alldata[,c(20,21,25)]
rnaia7<-alldata[,c(10,12:13)]
a7lsrdxck<-alldata[,c(38:40)]
a7lsrdx<-alldata[,c(35:37)]

sample<-list(a7lsrdx )
ck<-list(ck1,ck2,ck3)

for (i in 1:length(filenamebase)){
  for (m in 1:length(control))
{
  filename<-paste(filenamebase[[i]],control[[m]],sep="vs")

  
  x<-data.frame(sample[i],ck[m])
  
  
  group<-c(rep(filenamebase[[i]],length(sample[[i]])),
           rep(control[[m]],length(ck[[m]])))
  
  finddeg(x,group,filename)
 }
}


#define the data and group
#targets<-data.frame(TRAN=x[1:3],NT=x[4:6])


#defane the grorup and design
group<-factor(c(1,1,1,1,1,2,2,2,2))
finddeg<-function(x,group,filename){
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
de<-topTags(et,n=150000)



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
deg<-merge(alldegs,subset(funaa , funaa$Gene%in%alldegs[,1]),by="Gene",all=T)


write.csv(deg, file = paste(filename,"-DEG.csv",sep=""), row.names = F)
write.table(deg[,1], file=paste(filename,"-DEG-ID.txt",sep=""),quote=F,row.names = F)


}
###PCA
library(ggplot2)
library(ggbiplot)
plotpca<-function(x){

res <- prcomp(t(x), center = TRUE, scale = FALSE)
names(res)
#screeplot(res,type="lines")
temp<-predict(res) 
#plot(temp[,1:2])
#screeplot(temp,type="lines")

ggbiplot(res, obs.scale = 1, var.scale = 1,
         groups = group, ellipse = TRUE,var.axes = F,labels = names(x))
}
###PCA   PLOT

x_t <- t(alldata)

x <- as.data.frame(x_t)
x$group <- group
pca <- prcomp(t(rnaidata), scale=F)
percentVar <- pca$sdev^2 / sum( pca$sdev^2)
autoplot(pca, x, colour="group") + xlab(paste0("PC1 (", round(percentVar[1]*100), "% variance)")) + ylab(paste0("PC2 (", round(percentVar[2]*100), "% variance)")) + theme_bw() + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + theme(legend.position="right")

