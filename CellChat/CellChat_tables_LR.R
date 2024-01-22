## ---------------------------
##
## Script name: CellChat_tables_LR
##
## Purpose of script: 
## Obtain a table of all the receptor ligands of the signalling pathways: COLLAGEN, HH, SEMA3 and SPP1.
## The tables come from the entire CellChat database, NT, NOTO and NOTO SB conditions.
##
## Author: Bluwen Guidoux d'Halluin
##
## Date Created: 2023-
##
## ---------------------------
##
## Notes:
## R version 4.2.2
## CellChat 1.6.0
##
## ---------------------------
##
## Package

library(CellChat)

## ---------------------------

## Importing Dataset

#CellChat object

cellchat_NT_D7 <- readRDS(file = "data/cellchat_NT_iPS_NLC_trimean.rds")
cellchat_NOTO_D7 <- readRDS(file = "data/cellchat_NOTO_iPS_NLC_trimean.rds")
cellchat_NOTO_SB_D7 <- readRDS(file = "data/cellchat_NOTO_SB_iPS_NLC_trimean.rds")
CellChatDB <- CellChatDB.human

## ---------------------------

## Process the data
interaction_NT_D7 <- cellchat_NT_D7@LR$LRsig
interaction_NOTO_D7 <- cellchat_NOTO_D7@LR$LRsig
interaction_NOTO_SB_D7 <- cellchat_NOTO_SB_D7@LR$LRsig
interaction_all <- CellChatDB$interaction

COL9_NT_D7 <- interaction_NT_D7[grep("COL9", interaction_NT_D7$ligand),]
HH_NT_D7 <- interaction_NT_D7[grep("HH", interaction_NT_D7$ligand),]
SEMA3C_NT_D7 <- interaction_NT_D7[grep("SEMA3C", interaction_NT_D7$ligand),]
SPP1_NT_D7 <- interaction_NT_D7[grep("SPP1", interaction_NT_D7$ligand),]

COL9_NOTO_D7 <- interaction_NOTO_D7[grep("COL9", interaction_NOTO_D7$ligand),]
HH_NOTO_D7 <- interaction_NOTO_D7[grep("HH", interaction_NOTO_D7$ligand),]
SEMA3C_NOTO_D7 <- interaction_NOTO_D7[grep("SEMA3C", interaction_NOTO_D7$ligand),]
SPP1_NOTO_D7 <- interaction_NOTO_D7[grep("SPP1", interaction_NOTO_D7$ligand),]

COL9_NOTO_SB_D7 <- interaction_NOTO_SB_D7[grep("COL9", interaction_NOTO_SB_D7$ligand),]
HH_NOTO_SB_D7 <- interaction_NOTO_SB_D7[grep("HH", interaction_NOTO_SB_D7$ligand),]
SEMA3C_NOTO_SB_D7 <- interaction_NOTO_SB_D7[grep("SEMA3C", interaction_NOTO_SB_D7$ligand),]
SPP1_NOTO_SB_D7 <- interaction_NOTO_SB_D7[grep("SPP1", interaction_NOTO_SB_D7$ligand),]

COL9_all <- interaction_all[grep("COL9", interaction_all$ligand),]
HH_all <- interaction_all[grep("HH", interaction_all$ligand),]
SEMA3C_all <- interaction_all[grep("SEMA3C", interaction_all$ligand),]
SPP1_all <- interaction_all[grep("SPP1", interaction_all$ligand),]

## ---------------------------

## Save the data

write.csv(COL9_NT_D7, "results/Dataset_LR/NT_iPSC_NLC_COL9_interaction_CellChat.csv")
write.csv(HH_NT_D7, "results/Dataset_LR/NT_iPSC_NLC_HH_interaction_CellChat.csv")
write.csv(SEMA3C_NT_D7, "results/Dataset_LR/NT_iPSC_NLC_SEMA3C_interaction_CellChat.csv")
write.csv(SPP1_NT_D7, "results/Dataset_LR/NT_iPSC_NLC_SPP1_interaction_CellChat.csv")

write.csv(COL9_NOTO_D7, "results/Dataset_LR/NOTO_iPSC_NLC_COL9_interaction_CellChat.csv")
write.csv(HH_NOTO_D7, "results/Dataset_LR/NOTO_iPSC_NLC_HH_interaction_CellChat.csv")
write.csv(SEMA3C_NOTO_D7, "results/Dataset_LR/NOTO_iPSC_NLC_SEMA3C_interaction_CellChat.csv")
write.csv(SPP1_NOTO_D7, "results/Dataset_LR/NOTO_iPSC_NLC_SPP1_interaction_CellChat.csv")

write.csv(COL9_NOTO_SB_D7, "results/Dataset_LR/NOTO_SB_iPSC_NLC_COL9_interaction_CellChat.csv")
write.csv(HH_NOTO_SB_D7, "results/Dataset_LR/NOTO_SB_iPSC_NLC_HH_interaction_CellChat.csv")
write.csv(SEMA3C_NOTO_SB_D7, "results/Dataset_LR/NOTO_SB_iPSC_NLC_SEMA3C_interaction_CellChat.csv")
write.csv(SPP1_NOTO_SB_D7, "results/Dataset_LR/NOTO_SB_iPSC_NLC_SPP1_interaction_CellChat.csv")

write.csv(COL9_all, "results/Dataset_LR/All_COL9_interaction_CellChat.csv")
write.csv(HH_all, "results/Dataset_LR/All_HH_interaction_CellChat.csv")
write.csv(SEMA3C_all, "results/Dataset_LR/All_SEMA3C_interaction_CellChat.csv")
write.csv(SPP1_all, "results/Dataset_LR/All_SPP1_interaction_CellChat.csv")