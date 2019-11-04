suppressMessages(library('qqman'))
args<-commandArgs(T)
args<-c("E:\\temp\\KGL-4K.txt",5)
print("reading files")
gwas<-read.delim(args[1])
#####
gwasfilt<-gwas[-1,c(2:4,7)]
####


gwasfilt<-na.omit(gwasfilt)
print("omit")
names(gwasfilt)<-c("SNP","CHR","BP","P")
cbPalette <- c("#999999", "#E69F00", "#56B4E9", "#F0E442", "#0072B2", "#D55E00", "#CC79A7","#000000")

if(is.na(args[2])){
 print("Plotting ALL")
  
  pdf(paste(args[1],"-manhattan.PDF",sep=""))

manhattan(gwasfilt,col="gray",suggestiveline=FALSE,genomewideline=FALSE,highlight=snpOfIntersrNew)
dev.off() 
tiff(paste(args[1],"-manhattan.TIFF",sep=""),width = 1980,height = 1080,compression = "lzw",pointsize = 22)
manhattan(gwasfilt,col=cbPalette,suggestiveline=FALSE,genomewideline=FALSE,highlight=snpOfIntersrNew)
dev.off()

}else if(args[2]!=""){
  print("Plotting Chr")
  pdf(paste(args[1],"-CHR",args[2],"-manhattan.PDF",sep=""))
  
  manhattan(subset(gwasfilt,CHR==args[2]),col="#D55E00",suggestiveline=FALSE,genomewideline=FALSE)
  dev.off()
  
  tiff(paste(args[1],"-CHR",args[2],"-manhattan.TIFF",sep=""),width = 1980,height = 1080,compression = "lzw",pointsize = 22)
  manhattan(subset(gwasfilt,CHR==args[2]),col="#D55E00",suggestiveline=FALSE,genomewideline=FALSE)
  dev.off()
  
  
}


snp<-c("chr3.S_3400090","chr5.S_5064847","S2_237487658","chr6.S_131211595")

snpOfIntersr<-as.factor(subset(gwasfilt, BP> 3390100 & BP<3410100)[,-c(2:4)])

snpOfIntersr<-subset(gwasfilt,CHR==3& BP> 3390100 & BP<3410100,select =SNP)
class(snpOfIntersr)
snpsOfInterest<-as.vector(snpOfIntersr[,c(1:1)])
class(snpsOfInterest)
