# PACKAGES
# https://github.com/jokergoo/ComplexHeatmap/issues/32
library(rjson)
library(ComplexHeatmap)
library(dendextend)

setwd("/SCRATCH-BIRD/users/schevolleau/mapped_datasets/correction_seed_123/430/WGCNA/SF-14")

load("heatmap_WGCNA_TF-FALSE_SF-14env.Rdata")


# Retrieve trophoblast
pdf("garbage")
ht_1 <- htList[[1]]
temp <- column_dend(ht_1)
dev.off()

dendextend::cutree(temp$TE, k = 3)
