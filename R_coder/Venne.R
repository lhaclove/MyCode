setwd("D:\\lhac\\analysis\\a7l\\B2")
library(VennDiagram)
f1="12OE"
f2="12RNAi"
{a<- read.csv(paste(f1,"-DEG.csv",sep=""))
b<- read.csv(paste(f2,"-DEG.csv",sep=""))

a1<-a[,1]
b1<-b[,1]
set<-list(a=a1,b=b1)

names(set)<-c(f1,f2) 

venn.diagram(set,fill=c("red","blue"),paste(f1,"vs",f2,"out.tiff",sep=""))
write.csv(intersect(a1,b1),file=paste(f1,"vs",f2,"out.csv",sep=""))
}
