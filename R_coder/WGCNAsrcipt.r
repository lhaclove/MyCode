setwd("/public/home/liuhao/WGCNA/")
library(WGCNA)
enableWGCNAThreads()
allowWGCNAThreads()
zmData = read.csv("B2 RNASEQ ALL CPM.csv");
options(stringsAsFactors = FALSE);
zmData = read.csv("B2 RNASEQ ALL CPM.csv");

datExpr0 = as.data.frame(t(zmData[, -1]));
names(datExpr0) = zmData$Gene;
rownames(datExpr0) = names(zmData)[-1];

gsg = goodSamplesGenes(datExpr0, verbose = 3)
gsg$allOK
if (!gsg$allOK)
{
  if (sum(!gsg$goodGenes)>0)
    printFlush(paste("Removing genes:", paste(names(datExpr0)[!gsg$goodGenes], collapse = ", ")))
  if (sum(!gsg$goodSamples)>0)
    printFlush(paste("Removing samples:", paste(rownames(datExpr0)[!gsg$goodSamples], collapse = ", ")))
  datExpr0 = datExpr0[gsg$goodSamples, gsg$goodGenes]
}


keep <- apply(datExpr0, 2,  function(t){max(t)}) >4
datExpr0 <- datExpr0[keep]
keep <- apply(datExpr0, 2,  function(t){min(t)}) >4
datExpr0 <- datExpr0[keep]
dim(datExpr0)

powers = c(c(1:10), seq(from = 12, to=20, by=2))
sft = pickSoftThreshold(datExpr0, powerVector = powers, verbose = 3)

pdf("power.pdf")


plot(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
     xlab="Soft Threshold (power)",ylab="Scale Free Topology Model Fit,signed R^2",type="n",
     main = paste("Scale independence"));
text(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
     labels=powers,cex=cex1,col="red");
# this line corresponds to using an R^2 cut-off of h
abline(h=0.85,col="red")

dev.off()

pdf("sftpower.pdf")
plot(sft$fitIndices[,1], sft$fitIndices[,5],
     xlab="Soft Threshold (power)",ylab="Mean Connectivity", type="n",
     main = paste("Mean connectivity"))
text(sft$fitIndices[,1], sft$fitIndices[,5], labels=powers, cex=cex1,col="red")
dev.off()


net = blockwiseModules(datExpr0, power = 7,
                       TOMType = "unsigned", minModuleSize = 30,
                       reassignThreshold = 0, mergeCutHeight = 0.25,
                       numericLabels = TRUE, pamRespectsDendro = FALSE,
                       saveTOMs = TRUE,
                       saveTOMFileBase = "ZM",
                       verbose = 3,
                       maxBlockSize=20000)


table(net$colors)

mergedColors = labels2colors(net$colors)
pdf("modle.pdf")
plotDendroAndColors(net$dendrograms[[1]], mergedColors[net$blockGenes[[1]]],
                    "Module colors",
                    dendroLabels = FALSE, hang = 0.03,
                    addGuide = TRUE, guideHang = 0.05)
dev.off()
moduleLabels = net$colors
moduleColors = labels2colors(net$colors)
MEs = net$MEs;
geneTree = net$dendrograms[[1]];
save(MEs, moduleLabels, moduleColors, geneTree,
     file = "TEMP.RData")

TOM = TOMsimilarityFromExpr(datExpr0, power = 7);
nGenes = ncol(datExpr0)
nSamples = nrow(datExpr0)
dissTOM = 1-TOMsimilarityFromExpr(datExpr0, power = 7);
plotTOM = dissTOM^7
diag(plotTOM) = NA;
pdf("heatmap.pdf")
TOMplot(plotTOM, geneTree, moduleColors, main = "Network heatmap plot, all genes")
dev.off()

modules=c("black","blue","brown","cyan","darkgreen","darkgrey","darkorange","darkred","darkturquoise","green","greenyellow","grey","grey60","lightcyan","lightgreen","lightyellow","magenta","midnightblue","orange","pink","purple","red","royalblue","salmon","tan","turquoise","white","yellow")
probes = names(datExpr0) 

inModule = is.finite(match(moduleColors, modules));
modProbes = probes[inModule];
modTOM = TOM[inModule, inModule]; 
dimnames(modTOM) = list(modProbes, modProbes)
write.table(net$colors,file="model",quote=F)
write.table(probes,file="gene",quote=F)
cyt = exportNetworkToCytoscape(modTOM,  edgeFile = paste("CytoscapeInput-edges-all-070.txt", sep=""),  nodeFile = paste("CytoscapeInput-nodes-all-070.txt", sep=""), weighted = TRUE, threshold = 0.70, nodeNames = modProbes,  nodeAttr = moduleColors[inModule]);
save(TOM,file="TOM.RData") 
save.image(file="workspace.RData")
