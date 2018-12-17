list2<-read.delim2("d:\\RNA3\\tmp\\kegg.txt" ,header = T)

list1<-read.delim2("d:\\RNA3\\tmp\\a2.txt",header = T)

list3<-read.delim2("d:\\RNA3\\tmp\\v4keggid.txt",header = T)
list4<-read.delim2("d:\\RNA3\\tmp\\kegg.txt",header = T)
names(list1)<-c("ko","Aliases")
newlist<-merge(list2,list1
               ,by="kegg", all=T)

newlist1<-merge(newlist,list3
               ,by="Gene", all=T)

newlist2<-merge(newlist1,list4
                ,by="KEGG.Gene.ID",all = T)

write.csv(newlist,"a3.csv")
