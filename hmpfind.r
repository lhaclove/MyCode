a<-read.delim("F:\\ZhangXuHuan\\hmp\\5DAP-sample_SNP-merge.new.sorted.MAF0.05.hmp")
b<-read.table("F:\\ZhangXuHuan\\range.txt",header=TRUE)
e<-a[1,]
for(i in 1:length(b[[1]]))
{c<-subset(a , a$chrom ==b[i,1:3]$CHR & a$pos > b[i,1:3]$Start & a$pos <b[i,1:3]$End)
e<-rbind(e,c)} 

write.csv(e,"zhangxuhuan-snp.csv")
