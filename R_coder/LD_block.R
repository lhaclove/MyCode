ldblockall<-read.delim("E:\\ZXH\\hmp\\linked_snps_50kb.txt",stringsAsFactors = F)
ldblock<-ldblockall[which(ldblockall$CHR==5),]

ass<-read.delim("E:\\ZXH\\GWAS-zhangxuhuan-0531\\2017-empty\ rate\\P_with_Eff.tab",stringsAsFactors = F)





##process 2 linked snp (n=2) 
mergedBlock<-c("")
for (i in 1:length(ldblock[,1])){
  if(ldblock[i,]$NSNPS ==2){
    tmp<-unlist(strsplit(ldblock[i,]$SNPS,split = "[|]"))
    tmpp<-ass[which( ass$Marker %in% tmp),]
    if(sum(tmpp$Effect >  0) == 2 ){
    mEffect<-tmpp[which.max(tmpp$Effect),]$Marker
    }
    
    if (sum(tmpp$Effect <  0) == 2 ){
    mEffect<-tmpp[which.min(tmpp$Effect),]$Marker
    }
    
    if (sum(tmpp$Effect <  0) ==1){
      mEffect<-tmpp[which.min(tmpp$p),]$Marker
    }
    tmplist<-c(ldblock[i,]$CHR,ldblock[i,]$BP1,ldblock[i,]$BP2,ldblock[i,]$KB,ldblock[i,]$NSNPS,mEffect,tmpp[which(tmpp$Marker==mEffect),]$Effect)
    mergedBlock<-rbind(mergedBlock,tmplist)  
  }
}

colnames(mergedBlock)<-c("CHR","BP1","BP2","KB","N","marker","effect")
write.table(mergedBlock,file = "mergedBlock",col.names = F, sep = "\t",quote = F)

##process more (n>2)

mergedBlock1<-c("")
for (i in 1:length(ldblock[,1])){
  if(ldblock[i,]$NSNPS > 2){
    tmp<-unlist(strsplit(ldblock[i,]$SNPS,split = "[|]"))
    tmpp<-ass[which( ass$Marker %in% tmp),] 
    if(sum(tmpp$Effect >  0) > sum(tmpp$Effect <  0) ){

      mEffect<-tmpp[which.max(tmpp$Effect),]$Marker
      }else if(sum(tmpp$Effect >  0) < sum(tmpp$Effect <  0) ){

        mEffect<-tmpp[which.min(tmpp$Effect),]$Marker
     
      }else if(sum(tmpp$Effect >  0) == sum(tmpp$Effect <  0) ){
      
        mEffect<-tmpp[which.min(tmpp$P),]$Marker
    }
    tmplist<-c(ldblock[i,]$CHR,ldblock[i,]$BP1,ldblock[i,]$BP2,ldblock[i,]$KB,ldblock[i,]$NSNPS,mEffect,tmpp[which(tmpp$Marker==mEffect),]$Effect)
    mergedBlock1<-rbind(mergedBlock1,tmplist)  
      }
}
colnames(mergedBlock1)<-c("CHR","BP1","BP2","KB","N","marker","effect")
  write.table(mergedBlock1,file = "mergedBlock1",col.names = F,sep = "\t",quote = F)

#merge two type together
merBlock<-rbind(mergedBlock,mergedBlock1)  
write.table(merBlock,file = "merBlock",col.names = F,sep = "\t",quote = F)  
write.table(merBlock[,1:3],file = "merBlock",row.names =               F,sep = "\t",quote = F)  

