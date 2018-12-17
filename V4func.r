funaa<-read.csv("d:\\RNA3\\A7\\v4fun.csv" ,header = T)
namea<-read.csv("d:\\RNA3\\A7\\v4Genename.csv" ,header = T)
names(funaa)<-c("Gene","desc")
names(namea)<-c("Gene","Name")

list1<-read.table("d:\\new 2.txt",header = T)

write.csv(subset(funaa , funaa$Gene%in%a$Gene),"anno11111.csv")


anno<-subset(funaa , funaa$Gene%in%list1$chip5k)
a<-merge(list1,anno,by="Gene",all=T)
write.csv(merge(funaa,namea
      ,by="Gene",all = T),"v4genedecname.csv")
