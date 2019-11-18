all<-na.omit(read.delim2("e:/ZXH/文章数据准备/20190912/2017-B-all.txt",header = T))
names(all)<-c("marker","CHR","Pos","LD","HD","Ratio" )
snp<-unique(c("S9_144062061","chr5.S_5065457","chr5.S_5065465","chr5.S_5064847","chr5.S_5064848","chr5.S_5064850","chr5.S_5064851","chr2.S_2364890","S4_27786022","S4_27786044","S4_153528161","S4_197891503","S7_6545709","S7_6605771","chr8.S_157747971","chr9.S_146673482","chr9.S_146679111","chr9.S_146679247","chr9.S_146679248","chr9.S_146916031","chr9.S_146916678","chr9.S_146918747","chr9.S_146918803","chr9.S_146954380","chr9.S_146954523","S10_147968916","chr5.S_5064847","chr5.S_5064848","chr5.S_5064850","chr5.S_5064851"))
source("")

CM(all, plot.type="q", multracks=TRUE, ylim= c(4,10),threshold=c(7e-8,1e-5),threshold.lty=c(1,2),highlight = snp,
   threshold.lwd=c(1,1), threshold.col=c("black","grey"), amplify=TRUE,bin.size=1e6,
   chr.den.col=NULL, signal.col=c("red","green"),signal.cex=c(1,1),file="tiff",memo="",dpi=600)

CMplot(
  all, plot.type="m", multracks=TRUE, threshold=c(1e-7,1e-5),threshold.lty=c(1,2)
  threshold.lwd=c(1,1), ylim = c(4,10), threshold.col=c("black","grey"), amplify=TRUE,bin.size=1e6,bin.range=c(1,1024),
  chr.den.col=c("darkgreen", "yellow", "red"), signal.col=c("red","green"),signal.cex=c(1,1),file="jpg",memo="",dpi=600)

write.table(all,"2017-B-all.txt",row.names = FALSE, sep = "\t",quote = F)

CMplot(
  all, plot.type="m", multracks=TRUE, threshold=c(7e-8,1e-5),threshold.lty=c(1,2),
  threshold.lwd=c(1,1), ylim = c(4,10), threshold.col=c("black","grey"), amplify=TRUE,file="jpg",memo="",dpi=600)


CMplot(
  all, plot.type="m",multracks=TRUE, LOG10=TRUE, ylim= c(4,10), threshold=c(7e-8,1e-5),threshold.lty=c(1,2),highlight = snp, highlight.col = c("red"),
  threshold.lwd=c(1,1), threshold.col=c("black","grey"), amplify=TRUE,chr.den.col=NULL ,bin.size=1e6,
  signal.col=c("red","green"),signal.cex=c(1,1),
  signal.pch=c(19,19),file="jpg",memo="",dpi=300,file.output=TRUE,verbose=TRUE
)  


CMplot(  all,  plot.type="m", multracks=TRUE,
         ylim= c(4,10), threshold=c(7e-8,1e-5),threshold.lty=c(1,2),highlight = snp, highlight.col = c("red"),signal.col =c("red","yellow","green"),
         threshold.lwd=c(1,1), threshold.col=c("black","grey"), amplify=TRUE,chr.den.col=NULL ,bin.size=1e6,
        signal.cex=c(1,1),
         signal.pch=c(19,19),file="jpg",memo="",dpi=300,file.output=TRUE,verbose=TRUE
)  




###
CMplot(all,plot.type="q",col=c("dodgerblue1", "olivedrab3", "darkgoldenrod1"),threshold=1e-6,
       signal.pch=19,signal.cex=1.5,signal.col="red",conf.int=TRUE,box=FALSE,multracks=
         TRUE,file="jpg",memo="",dpi=300,file.output=TRUE,verbose=TRUE,ylim=c(0,8),width=5,height=5)
###


###
CMplot(all, plot.type="m",pch=1:3,multracks=TRUE,threshold=c(7e-8,1e-5),threshold.lty=c(1,2), 
       threshold.lwd=c(1,1), threshold.col=c("black","grey"), amplify=TRUE,bin.size=1e6,
       chr.den.col=NULL, signal.col=c("red","green"),signal.cex=c(1,1),
       file="jpg",memo="",dpi=300,file.output=TRUE,verbose=TRUE)
###