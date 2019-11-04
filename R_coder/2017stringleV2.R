setwd("D:\\RNA3\\report\\exp\\origanal_group")
library(edgeR)
funaa<-read.csv("d:\\RNA3\\v4fun.csv" ,header = T)
names(funaa)<-c("Gene","desc")
#read count data

alldata <- read.csv("A7Lgene_count_matrix.csv" ,row.names=1)
#alldata <- read.csv("transcript_count_matrix.csv",row.names=1)
#alldata <- read.delim("allcountB2.txt",row.names=1)


##orignal group
{filenamebase<-list("VP","SRDX","RNAi","OE")
control<-list("CK")

vp<-alldata[,c(1:4)]
a7srdx<-alldata[,c(5:9)]
a7rnai<-alldata[,c(10:13)]
a7oe<-alldata[,c(14:18)]
nt<-alldata[,c(19:22)]



sample<-list(vp, a7srdx,a7rnai,a7oe )
ck<-list(nt)
}

##modify group
{
  
  filenamebase<-list("A7-SRDX")
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
}


##A7L
{filenamebase<-list("OE","RNAi","SRDX")
control<-list("CK")
loeck<-alldata[,c(13:15)]
loe<-alldata[,c(16:18)]
lrnai<-alldata[,c(19:21)]
lrnaick<-alldata[,c(22:24)]
lsrdx<-alldata[,c(25:27)]
lsrdxck<-alldata[,c(28:30)]
sample<-list(loe,lrnai,lsrdx)
ck<-list(loeck,lrnaick,lsrdxck)}
for (i in 1:length(filenamebase)){
    filename<-paste(filenamebase[[i]],control[[1]],sep="vs")
    x<-data.frame(sample[i],ck[i])
    group<-c(rep(filenamebase[[i]],length(sample[[i]])),
    rep(control[[1]],length(ck[[i]])))
    finddeg(x,group,filename)
  }




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
keep <- rowSums(degs>1) >= 3 & ((degpval[,1]>=1 ) | (degpval[,1]<=-1) & degpval[,3]<0.05)
degs<-degs[keep,]
degpval<-degpval[keep,]
dim(degs)
if(length(degs)>0){
  

Gene<-rownames(degs)
da<-data.frame(Gene,degs)
db<-data.frame(Gene,degpval)
alldegs<-merge(da,db,by="Gene")
deg<-merge(alldegs,subset(funaa , funaa$Gene%in%alldegs[,1]),by="Gene",all=T)


write.csv(deg, file = paste(filename,"-DEG.csv",sep=""), row.names = F)
write.table(deg[,1], file=paste(filename,"-DEG-ID.txt",sep=""),quote=F,row.names = F)
}else{
  write.table("", file=paste(filename,"-has-no-deg.txt",sep=""),quote=F,row.names = F)
}

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
pca <- prcomp(t(x), scale=F)
percentVar <- pca$sdev^2 / sum( pca$sdev^2)
plotPCA(alldata)
autoplot(pca, x, colour="group") + xlab(paste0("PC1 (", round(percentVar[1]*100), "% variance)")) + ylab(paste0("PC2 (", round(percentVar[2]*100), "% variance)")) + theme_bw() + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + theme(legend.position="right")

