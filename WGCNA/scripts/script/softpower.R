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
datExpr <- t(readTable(snakemake@input[["EXPR"]]))

# SOFTPOWER ANALYSIS
powers <- c(as.numeric(snakemake@params[["MIN_SOFTPOWER"]]):as.numeric(snakemake@params[["MAX_SOFTPOWER"]]))
sft <- pickSoftThreshold(datExpr, powerVector = powers, verbose = 5,networkType = "signed",blockSize = as.numeric(snakemake@params[["MAX_BLOCKSIZE"]]))
sftAbstract<- -sign(sft$fitIndices[,3])*sft$fitIndices[,2]
names(sftAbstract)<-sft$fitIndices[,1]

pdf(paste0(snakemake@output))
	# Plot the results:
	cex1 = 0.9;
	# Scale-free topology fit index as a function of the soft-thresholding power
	plot(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2], xlab="Soft Threshold (power)",ylab="Scale Free Topology Model Fit,signed R^2", type="n",main = paste("Scale independence"));
	text(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2], labels=powers,cex=cex1,col="red"); 
	# this line corresponds to using an R^2 cut-off of h
	abline(h=0.90,col="red")
	# Mean connectivity as a function of the soft-thresholding power
	plot(sft$fitIndices[,1], sft$fitIndices[,5], xlab="Soft Threshold (power)",ylab="Mean Connectivity", type="n", main = paste("Mean connectivity"))
	text(sft$fitIndices[,1], sft$fitIndices[,5], labels=powers, cex=cex1,col="red")
dev.off()
