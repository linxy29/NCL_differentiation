## ---------------------------
##
## Script name: CellChat_visualiation_LR_interaction_chord_diagram
##
## Purpose of script: 
## Visualisation of ligands emitted by the notochord towards other clusters in the signalling pathways : 
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

pdf(file = "results/COLLAGEN/Chord_COLLAGEN_LR_NT_iPS_NLC.pdf", height = 15, width = 15)
netVisual_chord_gene(cellchat_NT_D7, signaling = c("COLLAGEN"),legend.pos.x = 8, color.use = col)
dev.off()

pdf(file = "results/HH/Chord_HH_LR_NT_iPS_NLC.pdf", height = 15, width = 15)
netVisual_chord_gene(cellchat_NT_D7, signaling = c("HH"),legend.pos.x = 8, color.use = col)
dev.off()

pdf(file = "results/SEMA3/Chord_SEMA3_LR_NT_iPS_NLC.pdf", height = 15, width = 15)
netVisual_chord_gene(cellchat_NT_D7, signaling = c("SEMA3"),legend.pos.x = 8, color.use = col)
dev.off()

# NOTO condition

pdf(file = "results/COLLAGEN/Chord_COLLAGEN_LR_NOTO_iPS_NLC.pdf", height = 15, width = 15)
netVisual_chord_gene(cellchat_NOTO_D7, sources.use = c("Notochord"), signaling = c("COLLAGEN"),legend.pos.x = 8, color.use = col)
dev.off()

pdf(file = "results/HH/Chord_HH_LR_NOTO_iPS_NLC.pdf", height = 15, width = 15)
netVisual_chord_gene(cellchat_NOTO_D7, sources.use = c("Notochord"), signaling = c("HH"),legend.pos.x = 8, color.use = col)
dev.off()

pdf(file = "results/SEMA3/Chord_SEMA3_LR_NOTO_iPS_NLC.pdf", height = 15, width = 15)
netVisual_chord_gene(cellchat_NOTO_D7, sources.use = c("Notochord"), signaling = c("SEMA3"),legend.pos.x = 8, color.use = col)
dev.off()

pdf(file = "results/SPP1/Chord_SPP1_LR_NOTO_iPS_NLC.pdf", height = 15, width = 15)
netVisual_chord_gene(cellchat_NOTO_D7, sources.use = c("Notochord"), signaling = c("SPP1"),legend.pos.x = 8, color.use = col)
dev.off()

# NOTO SB condition, before running the netVisual_chord_gene2 for SEMA3 be sure to run the script : CellChat_netVisual_chord_bug_correction

pdf(file = "results/COLLAGEN/Chord_COLLAGEN_LR_NOTO_SB_iPS_NLC.pdf", height = 15, width = 15)
netVisual_chord_gene(cellchat_NOTO_SB_D7, sources.use = c("Notochord"), signaling = c("COLLAGEN"),legend.pos.x = 8, color.use = col)
dev.off()

pdf(file = "results/HH/Chord_HH_LR_NOTO_SB_iPS_NLC.pdf", height = 15, width = 15)
netVisual_chord_gene(cellchat_NOTO_SB_D7, sources.use = c("Notochord"), signaling = c("HH"),legend.pos.x = 8, color.use = col)
dev.off()

pdf(file = "results/SEMA3/Chord_SEMA3_LR_NOTO_SB_iPS_NLC.pdf", height = 15, width = 15)
netVisual_chord_gene2(cellchat_NOTO_SB_D7, sources.use = c("Notochord"), signaling = c("SEMA3"),legend.pos.x = 8, color.use = col)
dev.off()

pdf(file = "results/SPP1/Chord_SPP1_LR_NOTO_SB_iPS_NLC.pdf", height = 15, width = 15)
netVisual_chord_gene(cellchat_NOTO_SB_D7, sources.use = c("Notochord"), signaling = c("SPP1"),legend.pos.x = 8, color.use = col)
dev.off()
