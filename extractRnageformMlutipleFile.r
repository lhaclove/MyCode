setwd("e:/lpf/SNP CALLING/")
{
  fl<-list.files("nnn")    
needrange<-read.delim2("E:/LPF/SNP CALLING/2.txt",header = F)
colnames(needrange)<-c("chr","start","end")

  nnnnn<- ""
for (n in 1:length(fl)){
  print(fl[n])
  pValuedata <- read.delim2(paste("e:/lpf/SNP CALLING/nnn/",fl[n],sep=""), stringsAsFactors = F)[,c(1:4,7)]

for (i in 1:10){
  assign(paste("hmpChr",i,sep=""),pValuedata[which(pValuedata$Chr == i),])
}

  all<-""
  for (r in 1:length(needrange[,1])){
    tmpp<-get(paste("hmpChr",needrange[r,]$chr,sep=""))[which( get(paste("hmpChr",needrange[r,]$chr,sep=""))$Pos > needrange[r,]$start & get(paste("hmpChr",needrange[r,]$chr,sep=""))$Pos < needrange[r,]$end ),]
    #print(paste("extracting ",list(needrange[r,]),sep=""))
    all<-rbind(all,tmpp)
  }
all<-all[ order(all$Marker),]
print(length(all[,1]))
nnnnn<- cbind(nnnnn,all$p)
print(length(nnnnn[,1]))
print(length(all$p))
}
nnnnn<- cbind(all[,2:4],nnnnn[,-1])

colnames(nnnnn)<-c(names(all[,2:4]),fl)
write.table(nnnnn[-1,],"LPF-2.np",row.names = F,quote = F,sep="\t")
}
