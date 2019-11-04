geneid<-read.delim2("D:\\BioTools\\Zea_mays.best.qPDB.download.txt" ,header = T)
id<-read.delim2("D:\\BioTools\\ID.txt",header = F)
names(id)<-"V3"
a<-subset(geneid,geneid$GeneID %in% id$V3 )


write.csv(a,file="qPrimerDB1.csv")
