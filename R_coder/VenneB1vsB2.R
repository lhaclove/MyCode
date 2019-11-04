setwd("D:\\lhac\\analysis\\Rtmp")
library(VennDiagram)
f1="RNAi"
f2="RNAi"
a<- read.csv(paste(f1,"B2normal-DEG.csv",sep=""))
b<- read.csv(paste(f2,"B1normal-DEG.csv",sep=""))

a1<-a[,1]
b1<-b[,1]
set<-list(a=a1,b=b1)

names(set)<-c(f1,f2) 

venn.diagram(set,fill=c("red","blue"),paste(f1,"B2vsB1",f2,"out.tiff",sep=""))
write.csv(intersect(a1,b1),file=paste(f1,"B2vsB1",f2,"out.csv",sep=""))
