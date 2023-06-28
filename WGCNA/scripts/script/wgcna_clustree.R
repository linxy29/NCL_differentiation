#$ -S /CONDAS/users/schevolleau/embryo-batch-effect_correction/bin/Rscript
#$ -e /LAB-DATA/BiRD/users/schevolleau/rlog/
#$ -o /LAB-DATA/BiRD/users/schevolleau/rlog/

setwd("/SCRATCH-BIRD/users/schevolleau/mapped_datasets/correction_seed_123/430/WGCNA/results")
print(paste0("WD : ", getwd()))

source("https://gitlab.univ-nantes.fr/E114424Z/veneR/raw/master/loadFun.R?inline=false")
source("https://gitlab.univ-nantes.fr/E198672Y/embryo-functions/-/raw/master/functions.R?inline=false")

# Packages
library(monocle3)
library(igraph)
library(clustree) # Not installed
library(RcppParallel)
library(reticulate)
library(RColorBrewer)
library(stringr)
library(tidyr)
print("Packages have been loaded")

# Data
wgcna_list <- list()
for (SF in 2:19){
    wgcna_list[[SF]] <- paste0("genesModulesCorrespondence_TF-FALSE_SF-", SF, ".tsv")
}
wgcna_genes <- readTable(wgcna_list[[2]], rnames=FALSE)
colnames(wgcna_genes)[2] <- "Module2"
for (i in 3:19){
    wgcna_genes[, paste0("Module",i)] <- readTable(paste0("genesModulesCorrespondence_TF-FALSE_SF-", i, ".tsv"), rnames=FALSE)$Module
}
writeTable(wgcna_genes, "leiden_cluster.tsv")
print("Datas have been loaded")

pdf("cluster_tree_wgcna.pdf")
clustree(wgcna_genes,prefix = "Module", node_text_size=0)
print("Clustree has been plotted")
dev.off()