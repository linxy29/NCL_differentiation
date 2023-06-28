#$ -S /CONDAS/users/schevolleau/embryo-wgcna/bin/Rscript
#$ -e /LAB-DATA/BiRD/users/schevolleau/rlog/
#$ -o /LAB-DATA/BiRD/users/schevolleau/rlog/

source("https://gitlab.univ-nantes.fr/E114424Z/veneR/raw/master/loadFun.R?inline=false")
source("https://gitlab.univ-nantes.fr/E198672Y/embryo-functions/-/raw/master/functions.R?inline=false")

# Load the WGCNA package
library(WGCNA)
library(jsonlite)
library(tidyr)
library(ggplot2)
library(gridExtra)

# VARIABLES
setwd("/SCRATCH-BIRD/users/schevolleau/mapped_datasets/correction_seed_123/430"); qPrint("WD : ", getwd())

## PATH
expr_Mat_path <- "/SCRATCH-BIRD/users/schevolleau/mapped_datasets/correction_seed_123/430/exprDat.normalized.log.corrected.tsv"
sample_annot_path <- "/SCRATCH-BIRD/users/schevolleau/mapped_datasets/correction_seed_123/430/sampleAnnotation.normalized.log.corrected.tsv"
genesMod_Mb_path <- "WGCNA/results/geneModuleMembership_14.tsv"
genesMod_Cr_path <- "WGCNA/results/genesModulesCorrespondence_TF-FALSE_SF-14.tsv"
genes_rank_path <- "WGCNA/geneRanking/iter_14/bisque4_geneRank_14.tsv"
rsave <- "WGCNA/rsave/all_14.RData"

# DATA
expr_Mat <- readTable(expr_Mat_path)
sample_annot <- readTable(sample_annot_path)
genesMod_Mb <- readTable(genesMod_Mb_path)
genesMod_Cr <- readTable(genesMod_Cr_path)
genes_rank <- readTable(genes_rank_path)

load(rsave)
colnames(MEs) <- substring(colnames(MEs), 3)

dir_plot <- "WGCNA/geneRanking/iter_14/plot"
dir.create(dir_plot)

genes_rank_path <- "WGCNA/geneRanking/iter_14"
list_mods <- setdiff(list.files(genes_rank_path), list.dirs(genes_rank_path, recursive = FALSE, full.names = FALSE))
mod_name <- list_mods[1]
for (mod_name in list_mods){
  name <- strsplit(mod_name, "_")[[1]][1]
  mod <- readTable(paste0(genes_rank_path, "/", mod_name))
  mod <- mod[order(mod$Membership, decreasing=TRUE),]
  mod$rn <- factor(rownames(mod), levels=rownames(mod))
 
  pdf(paste0(dir_plot, "/", name, ".pdf"))
  plot(ggplot(mod, aes(x= rn, y=Membership)) +
    geom_point(size=2, shape=23) + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank(), axis.line = element_line(colour = "black"), axis.title.x=element_blank(), axis.text.x=element_blank(), axis.ticks.x=element_blank()) + labs(x = "Genes", y = name))
  dev.off()
}
