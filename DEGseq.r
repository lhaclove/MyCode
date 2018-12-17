library('DEGseq')
library('edgeR')

args<-commandArgs(TRUE)

#get extra paramters
if(length(args)<4){
    print("program <readcount file> <samples, ck,wll> <groups, c,w> <outdir>")
    q(save = "no", status = 1, runLast = FALSE)
}
infile <- args[1]
samples <- unlist(strsplit(args[2],','))
groups<-unlist(strsplit(args[3],','))
outdir<-args[4]
outdir<-"d://"
infile <- "d:\\transcript_count_matrix.csv"
samples <- unlist(strsplit(args[2],','))
groups<-unlist(strsplit(args[3],','))
outdir<-args[4]

samples<-names(rawcount)

if(!file.exists(infile)){
    print("Can not find input file")
    q(save = "no", status = 2, runLast = FALSE)
}
alldata <- read.csv("transcript_count_matrix.csv",row.names=1)
#readcount normalization
rawcount<-read.table(infile,sep="\t", row.names=1, header=TRUE, check.names=FALSE, comment.char = "")
rawcount = rawcount[c(19:22,1:4)]
groups <-factor(c(1,1,1,1,2,2,2,2))
d<-DGEList(counts = rawcount,groups )

normData<-calcNormFactors(d)
normCount<-rawcount
for(i in 1:length(groups)){
    normCount[,i]<- alldata[,i]*1000000/((normData$samples[i,2])*(normData$samples[i,3]))
}
outdir="d://"
normCountFile<-paste(outdir,'/normCount',sep='')
#colnames(normCount) <- c(groups[1], groups[2])
write.table(data.frame(geneID=rownames(normCount),normCount,check.names=FALSE),file=normCountFile, sep='\t',quote=FALSE, row.name=FALSE, col.name=TRUE)

title<-readLines(file(normCountFile,'r'),n=1)
titleArr<-unlist(strsplit(title, '\t'))
index1<-which(titleArr==samples[1])
index2<-which(titleArr==samples[2])
geneExpFile<-c(normCountFile)
geneExpMatrix1<-readGeneExp(file=geneExpFile, geneCol=1, valCol=c(index1))
geneExpMatrix2<-readGeneExp(file=geneExpFile, geneCol=1, valCol=c(index2))
DEGexp(geneExpMatrix1 = geneExpMatrix1, geneCol1 = 1, expCol1 = 2, groupLabel1 = groups[1], geneExpMatrix2 = geneExpMatrix2, geneCol2 = 1, expCol2 = 2, groupLabel2 = groups[2], method = "MARS", thresholdKind=5, outputDir=outdir,  normalMethod="none")

outScore<-read.table(paste(outdir,'/output_score.txt',sep=''), header=TRUE, check.names=FALSE, row.names = 1)
output <- merge(x = rawcount[,samples], y=outScore[,c('log2(Fold_change) normalized','p-value','q-value(Storey et al. 2003)')], by='row.names')
colnames(output) <- c('gene_id', groups, 'log2FoldChange', 'pval', 'padj')
write.table(output, file=paste(outdir,"/",groups[1],'vs',groups[2],"_raw.GeneCompare.xls", sep=''),quote=FALSE,sep='\t',row.name=FALSE)

