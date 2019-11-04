library(gmodels)


pca=fast.prcomp(t(datExpr0))
pca
summary(pca)$importance
sizeGrWindow(12,9)
biplot(pca, cex=c(1.3, 1.2));

datExpr0.pr <- princomp(datExpr0,cor=TRUE)
screeplot(datExpr0.pr,type="lines")
