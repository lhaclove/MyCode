setwd("D:\\RNA3\\tmp\\test1\\")
funaa<-read.csv("d:\\RNA3\\A7\\v4fun.csv" ,header = T)
names(funaa)<-c("Gene","desc")




a<-read.csv("VP-CK1-DEG.csv",header = T)
b<-read.csv("SRDX-CK1-DEG.csv",header = T)
c<-read.csv("OE-CK1-DEG.csv",header = T)
d<-read.csv("RNAi-CK1-DEG.csv",header = T)
e<-intersect(intersect(a[,1],b[,1]),intersect(c[,1],d[,1]))

suba<-subset(a,a$Gene %in% e)
subb<-subset(b,b$Gene %in% e)
subc<-subset(c,c$Gene %in% e)
subd<-subset(d,d$Gene %in% e)

write.csv(merge(merge(suba,subb, by="Gene") ,merge(subc,subd,by="Gene"),by="Gene"),"overlap1111.csv")

vpsrdx<-intersect(a[,1],b[,1])
subvp<-subset(a,a$Gene %in% vpsrdx)
subsrdx<-subset(b,b$Gene %in% vpsrdx)
write.csv(merge(subvp,subsrdx, by="Gene",all=T),"VP-SRDX.csv")

oernai<-intersect(c[,1],d[,1])
OE<-subset(c,c$Gene %in% oernai)
RNAi<-subset(d,d$Gene %in% oernai)
ov<-merge(OE,RNAi, by="Gene",all=T)

ovanno<-subset(funaa , funaa$Gene%in%ov$Gene)
final<-merge(ov,ovanno,by="Gene")
  
write.csv(final,"OE-RNAi.csv")


chippeak<-read.csv("515R-only_peaks.narrowPeak.anno.xls.csv",header=T)
ovanno<-subset(funaa , funaa$Gene%in%chippeak$Gene)
final<-merge(chippeak,ovanno,by="Gene",all = T)
write.csv(final,"515R-only_peaks.narrowPeak.anno.csv")
