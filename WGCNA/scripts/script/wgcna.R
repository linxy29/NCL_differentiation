#$ -S /CONDAS/users/schevolleau/embryo-wgcna/bin/Rscript
#$ -e /LAB-DATA/BiRD/users/schevolleau/rlog/
#$ -o /LAB-DATA/BiRD/users/schevolleau/rlog/
#$ -pe smp 5
#$ -M simon.chevolleau@gmail.com
#$ -m beas

# PACKAGES
library(WGCNA)
library(jsonlite)
library(tidyr)
library(ggplot2)
library(gridExtra)

# FUNCTIONS
source("https://gitlab.univ-nantes.fr/E114424Z/veneR/raw/master/loadFun.R?inline=false")
source("https://gitlab.univ-nantes.fr/E198672Y/embryo-functions/-/raw/master/functions.R?inline=false")

# WD
setwd(snakemake@input[["WD"]])

# DIRECTORIES
dir.create("geneRanking", recursive=T)
dir.create("results", recursive=T)
dir.create("rsave", recursive=T)

# SETTINGS
options(stringsAsFactors = FALSE)
softPower <- as.numeric(snakemake@params[["SOFTPOWER"]])
datExpr <- t(readTable(snakemake@input[["EXPR"]]))
sample_annot <- readTable(snakemake@input[["ANNOTATION"]])

# WGCNA
# Major part of net parameters values come from D.Meistermann implementation
net = blockwiseModules(
	# Input data
	datExpr, 
	# Data checking options
	checkMissingData = TRUE,
	# Options for splitting data into blocks
	blocks = NULL,
	maxBlockSize = as.numeric(snakemake@params[["MAX_BLOCKSIZE"]]),
	blockSizePenaltyPower = 5,
	nPreclusteringCenters = as.integer(min(ncol(datExpr)/20, 100*ncol(datExpr)/as.numeric(snakemake@params[["MAX_BLOCKSIZE"]]))),
	randomSeed = 12345,
	# load TOM from previously saved file?
	loadTOM = FALSE,
	# Network construction arguments: correlation options
	corType = "pearson",
	maxPOutliers = 1, 
	quickCor = 0,
	pearsonFallback = "individual",
	cosineCorrelation = FALSE,
	# Adjacency function options
	power = softPower,
	networkType = "signed",
	replaceMissingAdjacencies = FALSE,
	# Topological overlap options
	TOMType = "signed",
	TOMDenom = "min",
	# Saving or returning TOM
	getTOMs = NULL,
	saveTOMs = FALSE, 
	saveTOMFileBase = "rsave/blockwiseTOM.RData",
	# Basic tree cut options
	deepSplit = 2,
	detectCutHeight = 0.995, 
	minModuleSize = min(20, ncol(datExpr)/2 ),
	# Advanced tree cut options
	maxCoreScatter = NULL, minGap = NULL,
	maxAbsCoreScatter = NULL, minAbsGap = NULL,
	minSplitHeight = NULL, minAbsSplitHeight = NULL,
	useBranchEigennodeDissim = FALSE,
	minBranchEigennodeDissim = 0.15,
	stabilityLabels = NULL,
	minStabilityDissim = NULL,
	pamStage = TRUE, pamRespectsDendro = TRUE,
	# Gene reassignment, module trimming, and module "significance" criteria
	reassignThreshold = 1e-6,
	minCoreKME = 0.5, 
	minCoreKMESize = min(20, ncol(datExpr)/2 )/3,
	minKMEtoStay = 0.3,
	# Module merging options
	mergeCutHeight = 0.15, 
	impute = TRUE, 
	trapErrors = FALSE, 
	# Output options
	numericLabels = FALSE,
	# Options controlling behaviour
	nThreads = 1,
	verbose = 1, indent = 0
)

save(net,file=paste0("rsave/net_", softPower,".RData"))

geneNames <- colnames(datExpr)
modNames<-substr(cn(net$MEs),3,60)
genesMod<-data.frame(Gene=colnames(datExpr),Module=net$colors)
write.table(genesMod , snakemake@output[["output_file"]],sep="\t",row.names=FALSE,col.names=TRUE,quote=FALSE)

# MEMBERSHIP OF GENES
geneModuleMembership <- as.data.frame(cor(datExpr, net$MEs, use = "p"))
colnames(geneModuleMembership) <- modNames

# MEMBERSHIP SAVE
ecrire(geneModuleMembership,paste0("results/geneModuleMembership_", softPower, ".tsv"))

# CONNECTIVITY
Alldegrees=intramodularConnectivity.fromExpr(datExpr, net$colors,networkType = "signed",power=softPower)
rownames(Alldegrees)<-geneNames
ecrire(Alldegrees,paste0("results/Alldegrees_",softPower, ".tsv"))
for(module in modNames ){
  genes<-genesMod$Gene[which(genesMod$Module==module)]
  d<-data.frame(row.names = genes, IntraConnectivity=Alldegrees[genes,"kWithin"],
                InterConnectivity=Alldegrees[genes,"kOut"],Membership=geneModuleMembership[genes,module])
  d<-d[order(d$Membership,decreasing = TRUE),]
  ecrire(d,paste0("geneRanking/",module,"_geneRank_", softPower,".tsv"))
}

## PLOT MODULE EIGEN GENES VALUES
pdf(paste0("MEs_", softPower, ".pdf"))
ModEig <- net$MEs
MEs_gathered <- gather(ModEig, "Module_eigen", "value")

plot(ggplot(MEs_gathered, aes(x=value, fill=Module_eigen)) + geom_density(alpha=.5))
plots <- list()
for (cluster in sample_annot$finalClusters){
	samples_cluster <- rownames(sample_annot)[sample_annot$finalClusters == cluster]
	mod_eig_cluster <-  ModEig[samples_cluster,]
	MEs_gathered_clusters <- gather(mod_eig_cluster, "Module_eigen", "value")
	plots[[cluster]] <- ggplot(MEs_gathered_clusters, aes(x=value, color=Module_eigen)) + geom_density() + labs(title=cluster)
}
marrangeGrob(plots, ncol=2, nrow=2)

save.image(paste0("rsave/all_", softPower,".RData"))
