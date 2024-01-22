## ---------------------------
##
## Script name: CellChat_iPS_NLC_trimean
##
## Purpose of script:
## Creation and analyse of a cellchat object for iPS-NLC dataset with trimean method
##
## Author: Bluwen Guidoux d'Halluin
##
## Date Created: 2023-19-01
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

#Data extraction
data.input <- GetAssayData(D7, assay = "RNA", slot = "data") # normalized data matrix
labels <- Idents(D7)
group = D7@meta.data[["annotation"]] # recover clusters with a resolution of res.0.4
meta <- data.frame(group = group, row.names = names(labels)) # create a dataframe of the cell labels

#Create a CellChat object
cellchat_d7 <- createCellChat(object = data.input, meta = meta, group.by = "group")

## ---------------------------

##Set the ligand-receptor interaction database

CellChatDB <- CellChatDB.human

#Set the used database in the object (here we use all CellChatDB)
cellchat_d7@DB <- CellChatDB

## ---------------------------

##Preprocessing the expression data for cell-cell communication analysis

#Subset the expression data of signaling genes for saving computation cost
cellchat_d7 <- subsetData(cellchat_d7) # This step is necessary even if using the whole database
cellchat_d7 <- identifyOverExpressedGenes(cellchat_d7)
cellchat_d7 <- identifyOverExpressedInteractions(cellchat_d7)

## ---------------------------

##Compute the communication probability and infer cellular communication network
cellchat_d7 <- computeCommunProb(cellchat_d7)

#Filter out the cell-cell communication if there are only few number of cells in certain cell groups (here nothing is filtered)
cellchat_d7 <- filterCommunication(cellchat_d7, min.cells = 10)

## ---------------------------

##Extract the inferred cellular communication network as a data frame

#Returns a data frame consisting of all the inferred cell-cell communications at the level of ligands/receptors.
df.net_pilot <- subsetCommunication(cellchat_d7)

#Returns a data frame consisting of all the inferred cell-cell communications at the level of signaling pathways.
df.net_pilot_pathway <- subsetCommunication(cellchat_d7, slot.name = "netP")

## ---------------------------

##Infer the cell-cell communication at a signaling pathway level
cellchat_d7 <- computeCommunProbPathway(cellchat_d7)

## ---------------------------

##Calculate the aggregated cell-cell communication network

cellchat_d7 <- aggregateNet(cellchat_d7) 
groupSize_pilot <- as.numeric(table(cellchat_d7@idents)) # number of cells in each cell group

col <- c('steelblue2', "tomato1", "orange2", "grey40", "ivory4", "indianred", "green4", "darkkhaki", "dodgerblue4")

#Visualize
par(mfrow = c(1,2), xpd=TRUE)
netVisual_circle(cellchat_d7@net$count, vertex.weight = groupSize_pilot, weight.scale = T, label.edge= F, title.name = "Number of interactions", color.use = col)
netVisual_circle(cellchat_d7@net$weight, vertex.weight = groupSize_pilot, weight.scale = T, label.edge= F, title.name = "Interaction weights/strength", color.use = col)

#Visualize for each cluster
mat <- cellchat_d7@net$weight
par(mfrow = c(3,4), xpd=TRUE)
for (i in 1:nrow(mat)) {
  mat2 <- matrix(0, nrow = nrow(mat), ncol = ncol(mat), dimnames = dimnames(mat))
  mat2[i, ] <- mat[i, ]
  netVisual_circle(mat2, vertex.weight = groupSize_pilot, weight.scale = T, edge.weight.max = max(mat), title.name = rownames(mat)[i], color.use = col)
}

## ---------------------------

##Visualize each signaling pathway using Circle plot

#Access all the signaling pathways showing significant communications
pathways.show.all <- cellchat_d7@netP$pathways

#Check the order of cell identity to set suitable vertex.receiver
levels(cellchat_d7@idents)

for (i in 1:length(pathways.show.all)) {
  #Visualize communication network associated with both signaling pathway and individual L-R pairs
  netVisual(cellchat_d7, signaling = pathways.show.all[i], layout = "circle",height=7, out.format = "pdf", color.use = col)
  
  #Compute and visualize the contribution of each ligand-receptor pair to the overall signaling pathway
  gg <- netAnalysis_contribution(cellchat_d7, signaling = pathways.show.all[i])
  ggsave(filename=paste0("results/Signaling_pathways/iPS_NLC_trimean/",pathways.show.all[i], "_L-R_contribution.pdf"), plot=gg, width = 6, height = 4, units = 'in', dpi = 300)
}

## ---------------------------

##Identify signaling roles (e.g., dominant senders, receivers) of cell groups as well as the major contributing signaling

#Compute the network centrality scores
cellchat_d7 <- netAnalysis_computeCentrality(cellchat_d7, slot.name = "netP") # the slot 'netP' means the inferred intercellular communication network of signaling pathways

#Visualize the computed centrality scores using heatmap, allowing ready identification of major signaling roles of cell groups
for (i in 1:length(pathways.show.all)) {
  pdf(file = paste0("results/Major_signaling_roles/iPS_NLC_trimean/all_pathways/",pathways.show.all[i],"_major_signaling_roles.pdf"))
  netAnalysis_signalingRole_network(cellchat_d7, signaling = pathways.show.all[i], width = 8, height = 2.5, font.size = 10, color.use = col)
  dev.off()
}

## ---------------------------

##Visualize the dominant senders (sources) and receivers (targets) in a 2D space

#Signaling role analysis on the aggregated cell-cell communication network from all signaling pathways
gg1 <- netAnalysis_signalingRole_scatter(cellchat_d7, color.use = col)
ggsave(filename=paste0("results/Major_signaling_roles/iPS_NLC_trimean/dominant_senders_and_receivers.pdf"), plot=gg1, width = 6, height = 4, units = 'in', dpi = 300)

##Identify signals contributing most to outgoing or incoming signaling of certain cell groups

# Signaling role analysis on the aggregated cell-cell communication network from all signaling pathways
pdf(file = "results/Major_signaling_roles/iPS_NLC_trimean/Heatmap_outgoing_signaling.pdf",height = 10)
netAnalysis_signalingRole_heatmap(cellchat_d7, pattern = "outgoing",height = 15, color.use = col)
dev.off()

pdf(file = "results/Major_signaling_roles/iPS_NLC_trimean/Heatmap_incoming_signaling.pdf",height = 10)
netAnalysis_signalingRole_heatmap(cellchat_d7, pattern = "incoming",height = 15, color.use = col)
dev.off()

## ---------------------------

##Identify global communication patterns to explore how multiple cell types and signaling pathways coordinate together

pdf(file="results/Communication_patterns/iPS_NLC_trimean/SelectK_outgoing.pdf", width=15)
selectK(cellchat_d7, pattern = "outgoing")
dev.off()

nPatterns = 3 # Both Cophenetic and Silhouette values begin to drop suddenly when the number of outgoing patterns is 3

cellchat_d7 <- identifyCommunicationPatterns(cellchat_d7, pattern = "outgoing", k = nPatterns,height = 15)

pdf(file = "results/Communication_patterns/iPS_NLC_trimean/outgoing_patern.pdf", width=15, height = 15)
identifyCommunicationPatterns(cellchat_d7, pattern = "outgoing", k = nPatterns,height = 15, color.use = col)
dev.off()

#River plot
pdf("results/Communication_patterns/iPS_NLC_trimean/rivier_plot_outgoing_patern.pdf", width=15)
netAnalysis_river(cellchat_d7, pattern = "outgoing", color.use = col)
dev.off()

#Dot plot
pdf("results/Communication_patterns/iPS_NLC_trimean/dot_plot_outgoing_patern.pdf", width=15)
netAnalysis_dot(cellchat_d7, pattern = "outgoing", color.use = col) 
dev.off()

## ---------------------------

##Identify and visualize incoming communication pattern of target cells

pdf(file="results/Communication_patterns/iPS_NLC_trimean/SelectK_incoming.pdf", width=15)
selectK(cellchat_d7, pattern = "incoming")
dev.off()

nPatterns = 3 # Both Cophenetic and Silhouette values begin to drop suddenly when the number of incoming patterns is 3

cellchat_d7 <- identifyCommunicationPatterns(cellchat_d7, pattern = "incoming", k = nPatterns,height = 15)

pdf(file = "results/Communication_patterns/iPS_NLC_trimean/incoming_patern.pdf", width=15, height = 15)
identifyCommunicationPatterns(cellchat_d7, pattern = "incoming", k = nPatterns, height = 15, color.use = col)
dev.off()

#River plot
pdf("results/Communication_patterns/iPS_NLC_trimean/rivier_plot_incoming_patern.pdf", width=15)
netAnalysis_river(cellchat_d7, pattern = "incoming", color.use = col)
dev.off()

#Dot plot
pdf("results/Communication_patterns/iPS_NLC_trimean/dot_plot_incoming_patern.pdf", width=15)
netAnalysis_dot(cellchat_d7, pattern = "incoming", color.use = col)
dev.off()

## ---------------------------

##Identify signaling groups based on their functional similarity
cellchat_d7 <- computeNetSimilarity(cellchat_d7, type = "functional")
cellchat_d7 <- netEmbedding(cellchat_d7, type = "functional")
cellchat_d7 <- netClustering(cellchat_d7, type = "functional", do.parallel = FALSE)

# Visualization in 2D-space
pdf(file = "results/Functional_structural/iPS_NLC_trimean/functional_similarity.pdf")
netVisual_embedding(cellchat_d7, type = "functional", label.size = 3.5)
dev.off()

##Identify signaling groups based on structure similarity
cellchat_d7 <- computeNetSimilarity(cellchat_d7, type = "structural")
cellchat_d7 <- netEmbedding(cellchat_d7, type = "structural")
cellchat_d7 <- netClustering(cellchat_d7, type = "structural", do.parallel = FALSE)

# Visualization in 2D-space
pdf(file = "results/Functional_structural/iPS_NLC_trimean/structural_similarity.pdf")
netVisual_embedding(cellchat_d7, type = "structural", label.size = 3.5)
dev.off()

## ---------------------------

##Save the CellChat object
saveRDS(cellchat_d7, file = "data/cellchat_iPS_NLC_trimean.rds")
