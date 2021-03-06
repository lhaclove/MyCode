
library(profvis)
##read ld file
#ld<-read.table("E:\\TASSEL5\\R_0.8",header = T)
#blocktmp<-c("")


tasellopt<-("E:/LPF/2017/KGL-ratio-1.txt")
effect<-("E:/LPF/2017/KGL-ratio-1 effect.txt")
ldblockall<-read.delim("E:\\ZXH\\hmp\\linked_snps_1kb.txt",stringsAsFactors = F)
single<-read.delim("E:\\ZXH\\hmp\\single_snps_1kb.txt",stringsAsFactors = F,header = T)
##read tassel optfile for p value
rec <- read.table(tasellopt,header = T,sep = "\t",stringsAsFactors = F)[-1,]

##read effect file

eff<-read.table(effect,head = T,stringsAsFactors = F)
eff<-eff[which(eff$Effect != 0),]

##merge p and effect
assemble<-merge(rec,eff,by = "Marker", all=T,sort= F)
ass<-assemble[,c(1,2,3,4,5,6,7,14,15,16,17,18,22,23,24)]

##A or B
ldblock<-ldblockall[which(ldblockall$CHR==5),]
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
    
  
  tmpp<-ass[which( get(paste("chr",ldblock[i,]$CHR,sep=""))$Marker %in% tmp),]
    
    
    
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

single.snp<-ass[ which(! ass$Marker%in% single$SNP),]
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
