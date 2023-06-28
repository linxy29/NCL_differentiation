#$ -S /LAB-DATA/BiRD/users/schevolleau/anaconda3/envs/batch-effect_correction/bin/Rscript
#$ -e /LAB-DATA/BiRD/users/schevolleau/rlog
#$ -o /LAB-DATA/BiRD/users/schevolleau/rlog

# PACKAGES
library(WGCNA)

# FUNCTIONS
source("https://raw.githubusercontent.com/DimitriMeistermann/veneR/main/loadFun.R?inline=false") 
source("https://gitlab.univ-nantes.fr/E198672Y/embryo-functions/-/raw/master/functions.R?inline=false")

setwd("/SCRATCH-BIRD/users/schevolleau/mapped_datasets/correction_seed_123/430/WGCNA/WGCNA/rsave_all_SF")
load("all_14.RData")

exprMat <- datExpr
sample_annotation <- lire("/SCRATCH-BIRD/users/schevolleau/mapped_datasets/correction_seed_123/430/sampleAnnotation.normalized.log.corrected.tsv")
res <- bicor(exprMat)
ecrire(res, "/SCRATCH-BIRD/users/schevolleau/mapped_datasets/correction_seed_123/430/WGCNA/WGCNA/14_correlation_from_bicor.tsv")

klf17_cor <- sort(res[, "KLF17"], decreasing = T)

klf17_cor_res <- data.frame(genes = names(klf17_cor), values = klf17_cor)
ecrire(klf17_cor_res, "/SCRATCH-BIRD/users/schevolleau/mapped_datasets/correction_seed_123/430/WGCNA/WGCNA/klf17_correlation.tsv")

sample_annotation_subset_EPI <- sample_annotation[sample_annotation$res.0.75 == "Epiblast",]
exprMat_susbset_EPI <- exprMat[rownames(sample_annotation_subset_EPI) ,]

res_subset_EPI <- bicor(exprMat_susbset_EPI)
ecrire(res_subset_EPI, "/SCRATCH-BIRD/users/schevolleau/mapped_datasets/correction_seed_123/430/WGCNA/WGCNA/correlation_from_bicor_susbset_EPI.tsv")

klf17_cor_subset_EPI <- sort(res_subset_EPI[, "KLF17"], decreasing = T)

klf17_cor_res_subset_EPI <- data.frame(genes = names(klf17_cor_subset_EPI), values = klf17_cor_subset_EPI)
ecrire(klf17_cor_res, "/SCRATCH-BIRD/users/schevolleau/mapped_datasets/correction_seed_123/430/WGCNA/WGCNA/klf17_correlation_subset_EPI.tsv")

dir.create("/SCRATCH-BIRD/users/schevolleau/mapped_datasets/correction_seed_123/430/WGCNA/correlation", recursive = T)
setwd("/SCRATCH-BIRD/users/schevolleau/mapped_datasets/correction_seed_123/430/WGCNA/correlation")

exprMat <- datExpr
sample_annotation <- lire("/SCRATCH-BIRD/users/schevolleau/mapped_datasets/correction_seed_123/430/sampleAnnotation.normalized.log.corrected.tsv")
mods <- readTable("/SCRATCH-BIRD/users/schevolleau/mapped_datasets/correction_seed_123/430/WGCNA/WGCNA/wgcna_modules_names_genes_14.tsv")

exprMat_row_mods <- exprMat[, rownames(mods)]
exprMat_row_mods <- cbind(mods["name"], t(exprMat_row_mods))
exprMat_row_mods[1:5,1:5]

# https://stats.stackexchange.com/questions/9918/how-to-compute-correlation-between-within-groups-of-variables
