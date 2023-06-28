#$ -S /LAB-DATA/BiRD/users/schevolleau/anaconda3/envs/wgcna/bin/Rscript
#$ -e /LAB-DATA/BiRD/users/schevolleau/rlog/
#$ -o /LAB-DATA/BiRD/users/schevolleau/rlog/

# WD
workingDir = "/SCRATCH-BIRD/users/schevolleau/mapped_datasets/correction_seed_123/430/WGCNA"
setwd(workingDir)

# Packages
source("https://gitlab.univ-nantes.fr/E198672Y/embryo-functions/-/raw/master/functions.R?inline=false")
source("https://gitlab.univ-nantes.fr/E114424Z/veneR/raw/master/loadFun.R?inline=false")
library(WGCNA)
library(ComplexHeatmap)
library(RJSONIO)
library(ggplot2)
library(tidyr)

# Data
load("rsave/all.RData")
TFdata <- readTable("/LAB-DATA/BiRD/users/schevolleau/git/Human_transcription_factor.tsv", rnames=FALSE)
TFdata <- TFdata[,1]

MEs <- net$MEs

print(paste0("There are ", length(table(genesMod$Module)), " modules"))
sdME<-apply(MEs,2,sd)
sortedMEs<-apply(MEs,2,sort,decreasing=TRUE)
diffMEs<-sortedMEs[1,]-sortedMEs[2,]

pdf("wgcna_cleaning.pdf")
barplot(sort(diffMEs),las=3,cex.names = .8)

goodModules<-names(diffMEs)[diffMEs<0.8]
goodModules <-  substr(goodModules, 3, nchar(goodModules))
genePerModule<-list()
tfPerModule <- list()
tfDatPerPerModule<-list()

for(module in goodModules){
	genePerModule[[module]]<-as.character(genesMod$Gene)[genesMod$Module==module]
	genePerModule[[module]]<-genePerModule[[module]][order(geneModuleMembership[genePerModule[[module]],module],decreasing = TRUE)]
	tfDatPerPerModule[[module]]<-intersect(genePerModule[[module]], TFdata)
	tfPerModule[[module]] <- tfDatPerPerModule[[module]]
	tfDatPerPerModule[[module]]$Membership<-geneModuleMembership[tfDatPerPerModule[[module]],module]
}

exportJSON <- toJSON(tfDatPerPerModule)
exportJSON <- toJSON(tfPerModule)
write(exportJSON, "tfPerModule.json")
write(exportJSON, "tfDatPerPerModule.json")
dev.off()

for (iter in 2:20){
	pdf(paste0("MEs_", iter, ".pdf"))
	rsave <- paste0("/SCRATCH-BIRD/users/schevolleau/mapped_datasets/correction_seed_123/430/WGCNA/rsave/all_", iter, ".RData")
	load(rsave)

	MEs <- net$MEs
	MEs_gathered <- gather(MEs, "Module_eigen", "value")

	ggplot(MEs_gathered, aes(x=value, fill=Module_eigen)) + geom_density(alpha=.5)

	for (ME in colnames(MEs)){
		plot(density(MEs[,ME]), main = ME)
	}
	dev.off()
}
getwd()