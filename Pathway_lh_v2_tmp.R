##NOTE
##need Zma_KEGG_LH.Rdata as KEGG data base.
##use NCBI id as genelist, and effect as value
load("e:/ZXH/database/Zma_KEGG_LH.Rdata")

## need SNP_assign_effect.R

##the key problem is generation of genesyml2effect file
b<-read.delim("E:\\ZXH\\pathway/\\CHR5_smaple_ID.txt")
names(b)<-c("v3","sxore")
c<-merge(idconvrt,b,by ="v3")
d<-c[,3:4]
d<-na.omit(d)
e<-d$sxore
names(e)<-d$ncbi
#length(a)







##KEGG and GO
hit<-names(e)
gscList <- list(GO_MF=GO_GL,PW_KEGG = PW_KEGG)
gscae <- new("GSCA", listOfGeneSetCollections=gscList,geneList = e, hits = hit)
gscae <- analyze(gscae, para = list(pValueCutoff = 0.05, pAdjustMethod = "BH", nPermutations = 1000, minGeneSetSize = 2,exponent = 1), doGSOA=TRUE, doGSEA=TRUE)
summarize(gscae)

gscae@result$GSEA.results$GO_MF

topGS_KEGG <- getTopGeneSets(gscae, "GSEA.results", c("PW_KEGG","GO_MF"), allSig=TRUE)
write.table(gscae@result$GSEA.results$GO_MF,file="GO",quote = F, sep="\t")
write.table(gscae@result$GSEA.results$PW_KEGG,file="KEGG",quote = F, sep="\t")

viewGSEA(gscae, "GO_MF", topGS_KEGG[["GO_MF"]][1])
viewGSEA(gscae, "PW_KEGG", topGS_KEGG[["PW_KEGG"]][1])

##MaizeCYC
f<-b$sxore
names(f)<-b$v3
hitf<-names(f)

gscList <- list(ZM_CYC=MC)
gscam <- new("GSCA", listOfGeneSetCollections=gscList,geneList = f, hits = hitf)
gscam <- analyze(gscam, para = list(pValueCutoff = 0.05, pAdjustMethod = "BH", nPermutations = 1000, minGeneSetSize = 2,exponent = 1), doGSOA=TRUE, doGSEA=TRUE)
topGS_MC<- getTopGeneSets(gscam, "GSEA.results", c("ZM_CYC"), allSig=TRUE)
#write.table(topGS_MC,file="topGS_MC",quote = F, s)
viewGSEA(gscam, "ZM_CYC", topGS_MC[["ZM_CYC"]][1])
write.table(gscam@result$GSEA.results,file="MaizeCyc",quote = F,  sep="\t")


##debug
save(GO_GL,PW_KEGG,MC,idconvrt,file="e:/ZXH/database/Zma_KEGG_LH.Rdata")
a<-runif(length(d$ncbi), min=-5, max=5)
names(a)<-d$ncbi
hit<-names(a)
gsca <- new("GSCA", listOfGeneSetCollections=gscList,geneList = a, hits = hit)
gsca <- analyze(gsca, para = list(pValueCutoff = 0.05, pAdjustMethod = "BH", nPermutations = 1000, minGeneSetSize = 10,exponent = 1), doGSOA=TRUE, doGSEA=TRUE)
summarize(gsca)
topGS_GO_MF <- getTopGeneSets(gsca, "GSEA.results", c("PW_KEGG"), allSig=TRUE)
viewGSEA(gsca, "PW_KEGG", topGS_GO_MF[["PW_KEGG"]][1])
