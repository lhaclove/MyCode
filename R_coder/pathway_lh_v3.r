load("e:/ZXH/database/Zma_KEGG_LH.Rdata")
pathwayGeneFilter=5
library(HTSanalyzeR)


  {
    tasellopt<-("E:/ZXH/文章数据准备/GWAS result/17_Bareness_ratio.xls")
effect<-("E:/ZXH/文章数据准备/GWAS result/17_Bareness_ratio_effect.txt")
ldblockall<-read.delim("E:\\ZXH\\hmp\\linked.snp",stringsAsFactors = F)
single<-read.delim("E:\\ZXH\\hmp\\single.snp",stringsAsFactors = F,header = F)
##read tassel optfile for p value
rec <- read.table(tasellopt,header = T,sep = "\t",stringsAsFactors = F)[-1,]

##read effect file

eff<-read.table(effect,head = T,stringsAsFactors = F)
eff<-eff[which(eff$Effect != 0),]

##merge p and effect
assemble<-merge(rec,eff,by = "Marker", all=T,sort= F)
ass<-assemble
#[,c(1,2,3,4,5,6,7,14,15,16,17,18,22,23,24)]

##A or B
#ldblock<-ldblockall[which(ldblockall$CHR==5),]
ldblock<-ldblockall
##
for (i in 1:10){
  assign(paste("chr",i,sep=""),ass[which(ass$Chr==i),])
}

##process 2 linked snp (n=2) 


mergeBlock<-c("CHR","BP1","BP2","KB","N","marker","effect")
a=1  
for (i in 1:length(ldblock[,1])){
  if (i %/% 1000 > a){
    print(a)
    a =a+1
  }
  
  tmp<-unlist(strsplit(ldblock[i,]$SNPS,split = "[|]"))
  
  
  tmpp<-get(paste("chr",ldblock[i,]$CHR,sep=""))[which( get(paste("chr",ldblock[i,]$CHR,sep=""))$Marker %in% tmp),]
  
  
  
  #ass<-ass[which( ! ass$Marker %in% tmp),]
  if(ldblock[i,]$NSNPS ==2){    
    if(sum(tmpp$Effect >  0) == 2 ){
      mEffect<-tmpp[which.max(tmpp$Effect),]$Marker
    }
    
    if (sum(tmpp$Effect <  0) == 2 ){
      mEffect<-tmpp[which.min(tmpp$Effect),]$Marker
    }
    
    if (sum(tmpp$Effect <  0) ==1 ){
      mEffect<-tmpp[which.min(tmpp$p),]$Marker
    }
  }
  if(ldblock[i,]$NSNPS > 2){
    if(sum(tmpp$Effect >  0) > sum(tmpp$Effect <  0) ){
      
      mEffect<-tmpp[which.max(tmpp$Effect),]$Marker
    }else if(sum(tmpp$Effect >  0) < sum(tmpp$Effect <  0) ){
      
      mEffect<-tmpp[which.min(tmpp$Effect),]$Marker
      
    }else if(sum(tmpp$Effect >  0) == sum(tmpp$Effect <  0) ){
      
      mEffect<-tmpp[which.min(tmpp$P),]$Marker
    }
    
  }
  tmplist<-c(ldblock[i,]$CHR,ldblock[i,]$BP1,ldblock[i,]$BP2,ldblock[i,]$KB,ldblock[i,]$NSNPS,mEffect,tmpp[which(tmpp$Marker==mEffect),]$Effect)
  mergeBlock<-rbind(mergeBlock,tmplist)
  
}


colnames(mergeBlock)<-mergeBlock[1,]
mergeBlock<-mergeBlock[-1,]

#write.table(mergeBlock,file = "merBlock.np",row.names = F,sep = "\t",quote = F)  
#write.table(mergeBlock[,c(1:3,6,7)],file = "mer1Block.np",row.names = F,sep = "\t",quote = F)  

single.snp<-ass[ which( ass$Marker %in% single$V1),]
single.snp<-data.frame(single.snp$Chr,single.snp$Pos,single.snp$Marker,single.snp$Effect)
#write.table(single.snp[,c(1,2,2,3,4)],file = "singleSNP.np",row.names = F,sep = "\t",quote = F)



ssnp<-single.snp[,c(1,2,2,3,4)]

lsnp<-mergeBlock[,c(1:3,6,7)]

row.names(lsnp)<-NULL


lsnp<-as.data.frame(lsnp,stringsAsFactors = F)
names(ssnp)<-names(lsnp)

allsnp<-rbind(lsnp,ssnp, stringsAsFactors = F)
allsnp<-allsnp[order(as.numeric(allsnp[,1]),as.numeric(allsnp[,2]),decreasing=FALSE),]
write.table(allsnp,file = paste(tasellopt,".allsnp.np",sep = ""),row.names = F,sep = "\t",quote = F)
}
##
##need annotated peak
##

## need SNP_assign_effect.R

##the key problem is generation of genesyml2effect file
b<-read.delim("./test1.np")
names(b)<-c("v3","sxore")
c<-merge(idconvrt,b,by ="v3")
d<-c[,3:4]
d<-na.omit(d)
e<-d$sxore
names(e)<-d$ncbi
#length(a)







##KEGG and GO
hit<-names(e)
gscList <- list(GO_MF=GO_GL,PW_KEGG = PW_KEGG)
gscae <- new("GSCA", listOfGeneSetCollections=gscList,geneList = e, hits = hit)
gscae <- analyze(gscae, para = list(pValueCutoff = 0.05, pAdjustMethod = "BH", nPermutations = 1000, minGeneSetSize = pathwayGeneFilter,exponent = 1), doGSOA=TRUE, doGSEA=TRUE)
summarize(gscae)

gscae@result$GSEA.results$GO_MF

topGS_KEGG <- getTopGeneSets(gscae, "GSEA.results", c("PW_KEGG","GO_MF"), allSig=TRUE)
write.table(gscae@result$GSEA.results$GO_MF,file="GO",quote = F, sep="\t")
write.table(gscae@result$GSEA.results$PW_KEGG,file="KEGG",quote = F, sep="\t")

viewGSEA(gscae, "GO_MF", topGS_KEGG[["GO_MF"]][1])
viewGSEA(gscae, "PW_KEGG", topGS_KEGG[["PW_KEGG"]][1])

##MaizeCYC
f<-b$sxore
names(f)<-b$v3
hitf<-names(f)

gscList <- list(ZM_CYC=MC)
gscam <- new("GSCA", listOfGeneSetCollections=gscList,geneList = f, hits = hitf)
gscam <- analyze(gscam, para = list(pValueCutoff = 0.05, pAdjustMethod = "BH", nPermutations = 1000, minGeneSetSize = pathwayGeneFilter,exponent = 1), doGSOA=TRUE, doGSEA=TRUE)
topGS_MC<- getTopGeneSets(gscam, "GSEA.results", c("ZM_CYC"), allSig=TRUE)
#write.table(topGS_MC,file="topGS_MC",quote = F, s)
summarize(gscam)
viewGSEA(gscam, "ZM_CYC", topGS_MC[["ZM_CYC"]][1])
write.table(gscam@result$GSEA.results,file="MaizeCyc",quote = F,  sep="\t")


##debug
save(GO_GL,PW_KEGG,MC,idconvrt,file="e:/ZXH/database/Zma_KEGG_LH.Rdata")
a<-runif(length(d$ncbi), min=-5, max=5)
names(a)<-d$ncbi
hit<-names(a)
gsca <- new("GSCA", listOfGeneSetCollections=gscList,geneList = a, hits = hit)
gsca <- analyze(gsca, para = list(pValueCutoff = 0.05, pAdjustMethod = "BH", nPermutations = 1000, minGeneSetSize = pathwayGeneFilter,exponent = 1), doGSOA=TRUE, doGSEA=TRUE)
summarize(gsca)
topGS_GO_MF <- getTopGeneSets(gsca, "GSEA.results", c("PW_KEGG"), allSig=TRUE)
viewGSEA(gsca, "PW_KEGG", topGS_GO_MF[["PW_KEGG"]][1])
