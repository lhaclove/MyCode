setwd("/project/home/lhaclove/zxh/test3/")
args<-commandArgs(T)
source("./PAST_lh_mod.r")


load("./LD.RData")
gff<-"./Zea_mays.AGPv3.27.gene.gff"
MycPT<-"./Zm_Myc4PAST.txt"
tasellopt<-"./17_Bareness_4k.xls"
effect<-"./17_Bareness_4k.xls"

##core for script
nc<-4

gwas_data_lh <-load_GWAS_data(tasellopt, effect,
                              association_columns =c("Trait", "Marker", "Chr", "Pos", "p","MarkerR2"),
                              effects_columns = c("Trait", "Marker", "Locus", "Site","Effect"))


save(gwas_data, file=paste("17_Bareness_4k.RData",sep=""))

genes <-assign_SNPs_to_genes(gwas_data,LD,gff,1000,0.8,nc)
save(genes, file=paste("17_Bareness_4k_Final_genes.RData",sep=""))
write.csv(genes,paste("17_Bareness_4k_Final_genes.csv",sep=""),row.names = FALSE)



rugplots_data <- find_pathway_significance(genes,MycPT,5,"increasing",1000,nc)

plot_pathways(rugplots_data,"pvalue",0.02,"increasing","../pastout1")


print("Done")