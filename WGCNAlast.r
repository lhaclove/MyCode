setwd("~/WGCNA")
library(WGCNA)
enableWGCNAThreads()
allowWGCNAThreads()
options(stringsAsFactors = FALSE);
zmData = read.csv("B2B3gene_count_matrix.csv",header = TRUE);
datExpr0 = as.data.frame(t(zmData[, -c(1,2)]));
datExpr0=datExpr0*1.0
names(datExpr0) = zmData$gene_id;
rownames(datExpr0) = names(zmData)[-c(1)];


#####
#Detection if all the sample contain the goodgenes
#and discard the sample with too many miss expression value
#
#
#####


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

#####
#Prievew the cluster of all the sample
#
#
#####



dimnames(datExpr0)[[1]]=names(data.frame(zmData[,-1]))
sampleTree = hclust(dist(datExpr0), method = "average");

sizeGrWindow(12,9)

par(cex = 0.6);
par(mar = c(0,4,2,0))
plot(sampleTree, main = "Sample clustering to detect outliers", sub="", xlab="", cex.lab = 1.5,
     cex.axis = 1.5, cex.main = 2)

abline(h = 1500000, col = "red");

clust = cutreeStatic(sampleTree, cutHeight = 1500000, minSize = 10)
table(clust) 
keepSamples = (clust==1)  
datExpr0 = datExpr0[keepSamples, ] 
nGenes = ncol(datExpr) 
nSamples = nrow(datExpr)



#####
#Load trait file
#trait file are create with trait and value
#this step is useless
#####
traitData = read.csv("3.csv");
dim(traitData)
names(traitData)

femaleSamples = rownames(datExpr0);
traitRows = match(femaleSamples, traitData$source);

datTraits = traitData[traitRows, -1];
rownames(datTraits) = traitData[traitRows, 1];
collectGarbage();

sampleTree2 = hclust(dist(datExpr), method = "average")

traitColors = numbers2colors(datTraits, signed = FALSE);

plotDendroAndColors(sampleTree2, traitColors,
                    groupLabels = names(datTraits),
                    main = "Sample dendrogram and trait heatmap")


#####
#Select the Power used for block deteced
#
#####
powers = c(c(1:10), seq(from = 12, to=20, by=2))
sft = pickSoftThreshold(datExpr0, powerVector = powers, verbose = 3)


#####
#Select the Power used for block deteced
#
#####

sizeGrWindow(9, 5)
par(mfrow = c(1,2));
cex1 = 0.9;

plot(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
     xlab="Soft Threshold (power)",ylab="Scale Free Topology Model Fit,signed R^2",type="n",
     main = paste("Scale independence"));
text(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
     labels=powers,cex=cex1,col="red");
# this line corresponds to using an R^2 cut-off of h
abline(h=0.85,col="red")
# Mean connectivity as a function of the soft-thresholding power
plot(sft$fitIndices[,1], sft$fitIndices[,5],
     xlab="Soft Threshold (power)",ylab="Mean Connectivity", type="n",
     main = paste("Mean connectivity"))
text(sft$fitIndices[,1], sft$fitIndices[,5], labels=powers, cex=cex1,col="red")


net = blockwiseModules(datExpr0, power = sft$powerEstimate,
                       TOMType = "unsigned", minModuleSize = 30,
                       reassignThreshold = 0, mergeCutHeight = 0.25,
                       numericLabels = TRUE, pamRespectsDendro = FALSE,
                       saveTOMs = TRUE,
                       saveTOMFileBase = "zmTOM",
                       verbose = 3,
                       maxBlockSize=15000)
table(net$colors)
sizeGrWindow(12, 9)
mergedColors = labels2colors(net$colors)
plotDendroAndColors(net$dendrograms[[1]], mergedColors[net$blockGenes[[1]]],
                    "Module colors",
                    dendroLabels = FALSE, hang = 0.03,
                    addGuide = TRUE, guideHang = 0.05)
moduleLabels = net$colors
moduleColors = labels2colors(net$colors)
MEs = net$MEs;
geneTree = net$dendrograms[[1]];
save(MEs, moduleLabels, moduleColors, geneTree,
     file = "test-9362geme-model.RData")



TOM = TOMsimilarityFromExpr(datExpr0, power = 7);

 

modules=c("black","blue","brown","cyan","darkgreen","darkgrey","darkorange","darkred","darkturquoise","green","greenyellow","grey","grey60","lightcyan","lightgreen","lightyellow","magenta","midnightblue","orange","pink","purple","red","royalblue","salmon","tan","turquoise","white","yellow")
probes = names(datExpr0) 

inModule = is.finite(match(moduleColors, modules));
modProbes = probes[inModule];
modTOM = TOM[inModule, inModule]; 
dimnames(modTOM) = list(modProbes, modProbes)

write.table(probes,net$colors,file="0918",quote=F)
cyt = exportNetworkToCytoscape(modTOM,  edgeFile = paste("CytoscapeInput-edges-all-050.txt", sep=""),  nodeFile = paste("CytoscapeInput-nodes-all-050.txt", sep=""), weighted = TRUE, threshold = 0.50, nodeNames = modProbes,  nodeAttr = moduleColors[inModule]);
save(TOM,file="lastdata.RData") 