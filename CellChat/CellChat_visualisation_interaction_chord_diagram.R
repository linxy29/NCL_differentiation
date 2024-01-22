## ---------------------------
##
## Script name: CellChat_visualiation_interaction_chord_diagram
##
## Purpose of script: 
## Visualisation of interactions inside signalling pathways : 
## COLLAGEN, HH, SEMA3 and SPP1, under NT, NOTO and NOTO SB conditions.
##
## Author: Bluwen Guidoux d'Halluin
##
## Date Created: 2023-
##
## ---------------------------
##
## Notes:
## R version 4.2.2
## Seurat 4.2.0
## tidyverse 1.3.2
## CellChat 1.6.0
## patchwork 1.1.2.9000
## NMF 0.25
## ggalluvial 0.12.3
##
## ---------------------------
##
## Package

library(Seurat)
library(CellChat)
library(patchwork)
options(stringsAsFactors = FALSE)
library(NMF)
library(ggalluvial)

## ---------------------------

## Importing Dataset

#CellChat object

cellchat_NT_D7 <- readRDS(file = "data/cellchat_NT_iPS_NLC_trimean.rds")
cellchat_NOTO_D7 <- readRDS(file = "data/cellchat_NOTO_iPS_NLC_trimean.rds")
cellchat_NOTO_SB_D7 <- readRDS(file = "data/cellchat_NOTO_SB_iPS_NLC_trimean.rds")

## ---------------------------

col <- c('steelblue2', "tomato1", "orange2", "grey40", "ivory4", "indianred", "green4", "darkkhaki", "dodgerblue4")

## Draw and save the chord diagram

# NT condition, there is no Notochord cluster and no SPP1 pathway

pdf(file = "results/COLLAGEN/Chord_COLLAGEN_NT_iPS_NLC.pdf", height = 15, width = 15)
netVisual_aggregate(cellchat_NT_D7, signaling = c("COLLAGEN"), layout = "chord", color.use = col)
dev.off()

pdf(file = "results/HH/Chord_HH_NT_iPS_NLC.pdf", height = 15, width = 15)
netVisual_aggregate(cellchat_NT_D7, signaling = c("HH"), layout = "chord", color.use = col)
dev.off()

pdf(file = "results/SEMA3/Chord_SEMA3_NT_iPS_NLC.pdf", height = 15, width = 15)
netVisual_aggregate(cellchat_NT_D7, signaling = c("SEMA3"), layout = "chord", color.use = col)
dev.off()

# NOTO condition, clusters who don't receive or send ligand to the Notochord cluster are in ligthgrey
# you can change that by just use color.use = col

pdf(file = "results/COLLAGEN/Chord_COLLAGEN_NOTO_iPS_NLC.pdf", height = 15, width = 15)
netVisual_aggregate(cellchat_NOTO_D7, signaling = c("COLLAGEN"), layout = "chord", color.use = col)
dev.off()

pdf(file = "results/HH/Chord_HH_NOTO_iPS_NLC.pdf", height = 15, width = 15)
netVisual_aggregate(cellchat_NOTO_D7, signaling = c("HH"), layout = "chord", color.use = c('steelblue2', "lightgrey", "orange2", "grey40", "ivory4", "lightgrey", "green4", "darkkhaki", "lightgrey"))
dev.off()

pdf(file = "results/SEMA3/Chord_SEMA3_NOTO_iPS_NLC.pdf", height = 15, width = 15)
netVisual_aggregate(cellchat_NOTO_D7, signaling = c("SEMA3"), layout = "chord", color.use = c('lightgrey', "tomato1", "orange2", "grey40", "lightgrey", "indianred", "green4", "darkkhaki", "lightgrey"))
dev.off()

pdf(file = "results/SPP1/Chord_SPP1_NOTO_iPS_NLC.pdf", height = 15, width = 15)
netVisual_aggregate(cellchat_NOTO_D7, signaling = c("SPP1"), layout = "chord", color.use = c('lightgrey', "tomato1", "orange2", "lightgrey", "lightgrey", "indianred", "green4", "darkkhaki", "dodgerblue4"))
dev.off()

# NOTO SB condition, clusters who don't receive or send ligand to the Notochord cluster are in ligthgrey
# you can change that by just use color.use = col

pdf(file = "results/COLLAGEN/Chord_COLLAGEN_NOTO_SB_iPS_NLC.pdf", height = 15, width = 15)
netVisual_aggregate(cellchat_NOTO_SB_D7, signaling = c("COLLAGEN"), layout = "chord", color.use = c('steelblue2', "tomato1", "orange2", "lightgrey", "ivory4", "lightgrey", "green4", "darkkhaki", "dodgerblue4"))
dev.off()

pdf(file = "results/HH/Chord_HH_NOTO_SB_iPS_NLC.pdf", height = 15, width = 15)
netVisual_aggregate(cellchat_NOTO_SB_D7, signaling = c("HH"), layout = "chord", color.use = c('steelblue2', "lightgrey", "orange2", "lightgrey", "lightgrey", "lightgrey", "green4", "darkkhaki", "dodgerblue4"))
dev.off()

pdf(file = "results/SEMA3/Chord_SEMA3_NOTO_SB_iPS_NLC.pdf", height = 15, width = 15)
netVisual_aggregate(cellchat_NOTO_SB_D7, signaling = c("SEMA3"), layout = "chord", color.use = c('lightgrey', "tomato1", "lightgrey", "lightgrey", "lightgrey", "lightgrey", "green4", "lightgrey", "lightgrey"))
dev.off()

pdf(file = "results/SPP1/Chord_SPP1_NOTO_SB_iPS_NLC.pdf", height = 15, width = 15)
netVisual_aggregate(cellchat_NOTO_SB_D7, signaling = c("SPP1"), layout = "chord", color.use = c('lightgrey', "tomato1", "lightgrey", "lightgrey", "lightgrey", "lightgrey", "green4", "lightgrey", "lightgrey"))
dev.off()