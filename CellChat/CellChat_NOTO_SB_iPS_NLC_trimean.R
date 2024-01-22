## ---------------------------
##
## Script name: CellChat_NOTO_SB_iPSC_NLC_trimean
##
## Purpose of script:
## Creation and analyse of a cellchat object for the NOTO SB condition with trimean method
##
## Author: Bluwen Guidoux d'Halluin
##
## Date Created: 2023-22-03
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

#Seurat object
D7 <- readRDS("data/ipsc_nc_diff_multi_annotated_sctrans_updated_22_03_23.RDS", refhook = NULL)

## ---------------------------

##Extract the CellChat input files from a Seurat V4 object

#Retrieve the NOTO_SB_D7 population
NOTO_SB_D7.subset <- subset(x = D7, subset = treatment %in% c("NOTOSB"))

#Data extraction
data_NOTO_SB_D7.input <- GetAssayData(NOTO_SB_D7.subset, assay = "RNA", slot = "data") # normalized data matrix
labels_NOTO_SB_D7 <- Idents(NOTO_SB_D7.subset)
group_NOTO_SB_D7  = NOTO_SB_D7.subset@meta.data[["annotation"]] # recover clusters with a resolution of res.0.4
meta_NOTO_SB_D7  <- data.frame(group = group_NOTO_SB_D7 , row.names = names(labels_NOTO_SB_D7)) # create a dataframe of the cell labels

#Create a CellChat object
cellchat_NOTO_SB_D7  <- createCellChat(object = data_NOTO_SB_D7.input, meta = meta_NOTO_SB_D7, group.by = "group")

## ---------------------------

##Set the ligand-receptor interaction database

CellChatDB <- CellChatDB.human

#Set the used database in the object (here we use all CellChatDB)
cellchat_NOTO_SB_D7@DB <- CellChatDB

## ---------------------------

##Preprocessing the expression data for cell-cell communication analysis

#Subset the expression data of signaling genes for saving computation cost
cellchat_NOTO_SB_D7 <- subsetData(cellchat_NOTO_SB_D7) # This step is necessary even if using the whole database
cellchat_NOTO_SB_D7 <- identifyOverExpressedGenes(cellchat_NOTO_SB_D7)
cellchat_NOTO_SB_D7 <- identifyOverExpressedInteractions(cellchat_NOTO_SB_D7)

## ---------------------------

##Compute the communication probability and infer cellular communication network

cellchat_NOTO_SB_D7 <- computeCommunProb(cellchat_NOTO_SB_D7)

#Filter out the cell-cell communication if there are only few number of cells in certain cell groups (here Primitive Streak-3 is filtered)
cellchat_NOTO_SB_D7 <- filterCommunication(cellchat_NOTO_SB_D7, min.cells = 10)

## ---------------------------

##Extract the inferred cellular communication network as a data frame

#Returns a data frame consisting of all the inferred cell-cell communications at the level of ligands/receptors.
df.net_NOTO_SB_D7 <- subsetCommunication(cellchat_NOTO_SB_D7)

#Returns a data frame consisting of all the inferred cell-cell communications at the level of signaling pathways.
df.net_NOTO_SB_D7_pathway <- subsetCommunication(cellchat_NOTO_SB_D7, slot.name = "netP")

## ---------------------------

##Infer the cell-cell communication at a signaling pathway level

cellchat_NOTO_SB_D7 <- computeCommunProbPathway(cellchat_NOTO_SB_D7)

## ---------------------------

##Calculate the aggregated cell-cell communication network

cellchat_NOTO_SB_D7 <- aggregateNet(cellchat_NOTO_SB_D7) 
groupSize_NOTO_SB_D7  <- as.numeric(table(cellchat_NOTO_SB_D7@idents)) # number of cells in each cell group

col <- c('steelblue2', "tomato1", "orange2", "grey40", "ivory4", "indianred", "green4", "darkkhaki", "dodgerblue4")

#Visualize
par(mfrow = c(1,2), xpd=TRUE)
netVisual_circle(cellchat_NOTO_SB_D7@net$count, vertex.weight = groupSize_NOTO_SB_D7, weight.scale = T, label.edge= F, title.name = "Number of interactions", color.use = col)
netVisual_circle(cellchat_NOTO_SB_D7@net$weight, vertex.weight = groupSize_NOTO_SB_D7, weight.scale = T, label.edge= F, title.name = "Interaction weights/strength", color.use = col)

#Visualize for each cluster
mat_NOTO_SB_D7 <- cellchat_NOTO_SB_D7@net$weight
par(mfrow = c(3,4), xpd=TRUE)
for (i in 1:nrow(mat_NOTO_SB_D7)) {
  mat2 <- matrix(0, nrow = nrow(mat_NOTO_SB_D7), ncol = ncol(mat_NOTO_SB_D7), dimnames = dimnames(mat_NOTO_SB_D7))
  mat2[i, ] <- mat_NOTO_SB_D7[i, ]
  netVisual_circle(mat2, vertex.weight = groupSize_NOTO_SB_D7, weight.scale = T, edge.weight.max = max(mat_NOTO_SB_D7), title.name = rownames(mat_NOTO_SB_D7)[i], color.use = col)
}

## ---------------------------

##Visualize each signaling pathway using Circle plot

#Access all the signaling pathways showing significant communications
pathways_NOTO_SB_D7.show.all <- cellchat_NOTO_SB_D7@netP$pathways

#Check the order of cell identity to set suitable vertex.receiver
levels(cellchat_NOTO_SB_D7@idents)

for (i in 1:length(pathways_NOTO_SB_D7.show.all)) {
  #Visualize communication network associated with both signaling pathway and individual L-R pairs
  netVisual(cellchat_NOTO_SB_D7, signaling = pathways_NOTO_SB_D7.show.all[i], layout = "circle",height=7, out.format = "pdf", color.use = col)
  
  #Compute and visualize the contribution of each ligand-receptor pair to the overall signaling pathway
  gg <- netAnalysis_contribution(cellchat_NOTO_SB_D7, signaling = pathways_NOTO_SB_D7.show.all[i])
  ggsave(filename=paste0("results/Signaling_pathways/NOTO_SB_D7_trimean/",pathways_NOTO_SB_D7.show.all[i], "_L-R_contribution.pdf"), plot=gg, width = 6, height = 4, units = 'in', dpi = 300)
}

## ---------------------------

##Identify signaling roles (e.g., dominant senders, receivers) of cell groups as well as the major contributing signaling

#Compute the network centrality scores
cellchat_NOTO_SB_D7 <- netAnalysis_computeCentrality(cellchat_NOTO_SB_D7, slot.name = "netP") # the slot 'netP' means the inferred intercellular communication network of signaling pathways

#Visualize the computed centrality scores using heatmap, allowing ready identification of major signaling roles of cell groups
for (i in 1:length(pathways_NOTO_SB_D7.show.all)) {
  pdf(file = paste0("results/Major_signaling_roles/NOTO_SB_D7_trimean/all_pathways/",pathways_NOTO_SB_D7.show.all[i], "_major_signaling_roles.pdf"))
  netAnalysis_signalingRole_network(cellchat_NOTO_SB_D7, signaling = pathways_NOTO_SB_D7.show.all[i], width = 8, height = 2.5, font.size = 10, color.use = col)
  dev.off()
}

## ---------------------------

##Visualize the dominant senders (sources) and receivers (targets) in a 2D space

#Signaling role analysis on the aggregated cell-cell communication network from all signaling pathways
gg1 <- netAnalysis_signalingRole_scatter(cellchat_NOTO_SB_D7, color.use = col)
ggsave(filename=paste0("results/Major_signaling_roles/NOTO_SB_D7_trimean/dominant_senders_and_receivers.pdf"), plot=gg1, width = 6, height = 4, units = 'in', dpi = 300)

##Identify signals contributing most to outgoing or incoming signaling of certain cell groups

# Signaling role analysis on the aggregated cell-cell communication network from all signaling pathways
pdf(file = "results/Major_signaling_roles/NOTO_SB_D7_trimean/Heatmap_outgoing_signaling.pdf",height = 10)
netAnalysis_signalingRole_heatmap(cellchat_NOTO_SB_D7, pattern = "outgoing",height = 15, color.use = col)
dev.off()

pdf(file = "results/Major_signaling_roles/NOTO_SB_D7_trimean/Heatmap_incoming_signaling.pdf",height = 10)
netAnalysis_signalingRole_heatmap(cellchat_NOTO_SB_D7, pattern = "incoming",height = 15, color.use = col)
dev.off()

## ---------------------------

##Identify global communication patterns to explore how multiple cell types and signaling pathways coordinate together

pdf(file="results/Communication_patterns/NOTO_SB_D7_trimean/SelectK_outgoing.pdf", width=15)
selectK(cellchat_NOTO_SB_D7, pattern = "outgoing")
dev.off()

nPatterns = 3 # Both Cophenetic and Silhouette values begin to drop suddenly when the number of outgoing patterns is 3
cellchat_NOTO_SB_D7 <- identifyCommunicationPatterns(cellchat_NOTO_SB_D7, pattern = "outgoing", k = nPatterns)

pdf(file = "results/Communication_patterns/NOTO_SB_D7_trimean/outgoing_patern.pdf", width=15, height = 15)
identifyCommunicationPatterns(cellchat_NOTO_SB_D7, pattern = "outgoing", k = nPatterns,height = 15, color.use = c('steelblue2', "tomato1", "orange2", "grey40", "green4", "darkkhaki", "dodgerblue4"))
dev.off()

#River plot
pdf("results/Communication_patterns/NOTO_SB_D7_trimean/rivier_plot_outgoing_patern.pdf", width=15)
netAnalysis_river(cellchat_NOTO_SB_D7, pattern = "outgoing", color.use = c('steelblue2', "tomato1", "orange2", "grey40", "green4", "darkkhaki", "dodgerblue4"))
dev.off()

#Dot plot
pdf("results/Communication_patterns/NOTO_SB_D7_trimean/dot_plot_outgoing_patern.pdf", width=15)
netAnalysis_dot(cellchat_NOTO_SB_D7, pattern = "outgoing", color.use = c('steelblue2', "tomato1", "orange2", "grey40", "green4", "darkkhaki", "dodgerblue4"))
dev.off()

## ---------------------------

##Identify and visualize incoming communication pattern of target cells

pdf(file="results/Communication_patterns/NOTO_SB_D7_trimean/SelectK_incoming.pdf", width=15)
selectK(cellchat_NOTO_SB_D7, pattern = "incoming")
dev.off()

nPatterns = 3 # Both Cophenetic and Silhouette values begin to drop suddenly when the number of outgoing patterns is 3
cellchat_NOTO_SB_D7 <- identifyCommunicationPatterns(cellchat_NOTO_SB_D7, pattern = "incoming", k = nPatterns)

pdf(file = "results/Communication_patterns/NOTO_SB_D7_trimean/incoming_patern.pdf", width=15, height = 15)
identifyCommunicationPatterns(cellchat_NOTO_SB_D7, pattern = "incoming", k = nPatterns, height = 15, color.use = c('steelblue2', "tomato1", "orange2", "grey40", "green4", "darkkhaki", "dodgerblue4"))
dev.off()

#River plot
pdf("results/Communication_patterns/NOTO_SB_D7_trimean/rivier_plot_incoming_patern.pdf", width=15)
netAnalysis_river(cellchat_NOTO_SB_D7, pattern = "incoming", color.use = c('steelblue2', "tomato1", "orange2", "grey40", "green4", "darkkhaki", "dodgerblue4"))
dev.off()

#Dot plot
pdf("results/Communication_patterns/NOTO_SB_D7_trimean/dot_plot_incoming_patern.pdf", width=15)
netAnalysis_dot(cellchat_NOTO_SB_D7, pattern = "incoming", color.use = c('steelblue2', "tomato1", "orange2", "grey40", "green4", "darkkhaki", "dodgerblue4"))
dev.off()

## ---------------------------

##Identify signaling groups based on their functional similarity
cellchat_NOTO_SB_D7 <- computeNetSimilarity(cellchat_NOTO_SB_D7, type = "functional")
cellchat_NOTO_SB_D7 <- netEmbedding(cellchat_NOTO_SB_D7, type = "functional")
cellchat_NOTO_SB_D7 <- netClustering(cellchat_NOTO_SB_D7, type = "functional", do.parallel = FALSE)

# Visualization in 2D-space
pdf(file = "results/Functional_structural/NOTO_SB_D7_trimean/functional_similarity.pdf")
netVisual_embedding(cellchat_NOTO_SB_D7, type = "functional", label.size = 3.5)
dev.off()

##Identify signaling groups based on structure similarity
cellchat_NOTO_SB_D7 <- computeNetSimilarity(cellchat_NOTO_SB_D7, type = "structural")
cellchat_NOTO_SB_D7 <- netEmbedding(cellchat_NOTO_SB_D7, type = "structural")
cellchat_NOTO_SB_D7 <- netClustering(cellchat_NOTO_SB_D7, type = "structural", do.parallel = FALSE)

# Visualization in 2D-space
pdf(file = "results/Functional_structural/NOTO_SB_D7_trimean/structural_similarity.pdf")
netVisual_embedding(cellchat_NOTO_SB_D7, type = "structural", label.size = 3.5)
dev.off()

## ---------------------------

##Save the CellChat object
saveRDS(cellchat_NOTO_SB_D7, file = "data/cellchat_NOTO_SB_iPS_NLC_trimean.rds")
