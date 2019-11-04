library('qqman')
args<-commandArgs(T)
gwas<-read.delim(args[1])
gwasfilt<-gwas[-1,c(2:4,7)]
names(gwasfilt)<-c("SNP","CHR","BP","P")
cbPalette <- c("#999999", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")
pdf(paste(args[1],"-manhattan.PDF",sep=""))

manhattan(gwasfilt,col=cbPalette)
dev.off()