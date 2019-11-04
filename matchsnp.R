setwd(".//")
snpneed<-read.table("f://新建文件夹//snp.txt",header=T)
snp<-read.delim2("F:\\zhangxuhuan\\287_inbreds_genotypic_data\\5DAP-sample_SNP-merge.new.sorted.MAF0.05.hmp",header=T)


e<-a[1,]
for(i in 1:length(b[[1]]))
  {c<-subset(snp , snp$chrom == snpneed[i,]$chr & snp$pos == snpneed[i,]$V3)
e<-rbind(e,c)} 


write.csv(e,"lengfengfei-snp.csv")

subset(snp , snp$chrom == "9" & snp$pos == "147417996")
