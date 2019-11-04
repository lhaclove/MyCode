suppressMessages(library('qqman'))
cbPalette <- c("#999999", "#E69F00", "#56B4E9", "#F0E442", "#0072B2", "#D55E00", "#CC79A7","#000000")

args<-commandArgs(T)
args<-c("E:\\ZXH\\GWAS-zhangxuhuan-0531\\2017-empty\ rate\\2017_4k_emp_rate_wmlm.txt",5)
print("reading files")
plotcutoff=1
gwas<-read.delim(args[1])
gwasfilt<-gwas[-1,c(2:4,6)]
names(gwasfilt)<-c("SNP","CHR","BP","P")

#####
gwasfilt<-gwas[-1,c(2:4,6)]
eff<-abs(gwas[-1,6])
gwasfilt<-data.frame(gwasfilt,eff)
####

eff.file<-which(abs(gwasfilt["P"])>0 )
marker.sig <- which(gwasfilt["P"] < cutoff & gwasfilt["P"] != 0 )
gwasfilt.sig <- gwasfilt[eff.file,]
rec.sig<-rec[marker.sig,]
eff<-abs(gwasfilt.sig[,4])
gwasfilt<-data.frame(gwasfilt.sig[1:3],eff)
###

gwasfiltna<-na.omit(gwasfilt)
print("omit")
marker <- which(-log10(gwasfiltna["P"] )> plotcutoff )
gwasfilt <- gwasfiltna[marker,]



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


snp<-c("chr5.S_5064847","chr5.S_5064847","chr9.S_146918747")

snpOfIntersr<-as.factor(subset(gwasfilt, BP> 3390100 & BP<3410100)[,-c(2:4)])

snpOfIntersr<-subset(gwasfilt,select =SNP)
class(snpOfIntersr)
snpsOfInterest<-as.vector(snpOfIntersr[,c(1:1)])
class(snpsOfInterest)
