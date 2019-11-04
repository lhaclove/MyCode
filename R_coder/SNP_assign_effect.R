library(tidyverse)
chr=5

##read ld file
ld<-read.table("E:\\TASSEL5\\R_0.8",header = T)
blocktmp<-c("")

##read tassel optfile for p value
tasellopt<-("E:\\ZXH\\tmp/18_Ftime.p.txt")
rec <- read.table(tasellopt,header = T,sep = "\t")[-1,]


##read effect file
effect<-("E:\\ZXH\\tmp\\18_Ftime.effect.txt")
eff<-read.table(effect,head = T)
eff<-eff[which(eff$Effect != 0),]

##merge p and effect
assemble<-merge(rec,eff,by = "Marker", all=T,sort= F)
ass<-assemble[,c(1,2,3,4,5,6,7,14,15,16,17,18,22,23,24)]
write.table(ass,file="P_with_Eff.tab",quote = F,sep= "\t" )



## A or B
asschr1<-ass[which(ass$Chr==chr),]
asschr1<-ass


###

##linked snp means snp which has a R^2 larger than 0.8 to another snp
linked.snp<-unique(rbind(asschr1[which(asschr1$Chr %in% ld$Locus1 & asschr1$Pos %in% ld$Position1),],asschr1[which(asschr1$Chr %in% ld$Locus2 & asschr1$Pos %in% ld$Position2),]))

##single snp means not linked (type0), and these snps can assign the effect and p value directly
single.snp<-asschr1[which(!( asschr1$Marker%in% linked.snp$Marker)),]

chr1linkedLd<-ld[which(ld$Locus1 %in% linked.snp$chr & ld$Position1 %in% linked.snp$Pos),]
write.table(chr1linkedLd,file="chr5linkedLd",quote = F,sep= "\t" ,row.names = F)

##filter type1 linked snp
a<-chr1linkedLd[,2]
type1snp<-chr1linkedLd[which(a %in% unique(setdiff(a,a[duplicated(a)]))),]

##filter type1 linked snp
type2snp<-chr1linkedLd[which(!(a %in% unique(setdiff(a,a[duplicated(a)])))),]




plot(NULL,type="n",axes=F,main="",xlim=c(min(type2snp[,2])/100000,max(type2snp[,2])/100000),ylim=c(1,4),xlab="",ylab="")

rect_h=4
rect(min(type2snp[,2])/100000,rect_h,max(type2snp[,2])/100000,rect_h-1,lwd=0.8)
len_sig_m<-length(type2snp[,2])
for (n in 1:len_sig_m){
  segments(type2snp[,2][n]/100000,rect_h,(min(type2snp[,2])+((n-1)*(max(type2snp[,2])-min(type2snp[,2]))/len_sig_m))/100000,rect_h-1,lwd=0.8)
}




