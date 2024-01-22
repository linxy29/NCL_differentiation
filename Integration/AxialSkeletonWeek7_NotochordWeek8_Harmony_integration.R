# TBXT is in wk3 and wk4, not wk2 (no TBXT at all), not variable gene in sample 1 as well (4 cells only)
#filtered means filtering with this:  subset = nFeature_RNA > 200 & nFeature_RNA < 6000 & percent.mt < 10

unfiltered_7wk2 <- readRDS("/usersdata/share/linxy/2023-04-01_unfiltered_sevenweek2.RDS")


# normalize 
normalized_7wk2 <- NormalizeData(unfiltered_7wk2, normalization.method = "LogNormalize")
normalized_7wk2 <- FindVariableFeatures(normalized_7wk2, selection.method = "vst", nfeatures = 2000)

# scale
all.genes.wk2 <- rownames(normalized_7wk2)
normalized_7wk2 <- ScaleData(normalized_7wk2, features = all.genes.wk2)



> sum(GetAssayData(object = unfiltered_7wk4, slot = "data", assay="RNA")["TBXT",]>0)
[1] 10
> sum(GetAssayData(object = unfiltered_7wk4, slot = "data", assay="RNA")["SHH",]>0)
[1] 28
> sum(GetAssayData(object = unfiltered_7wk4, slot = "data", assay="RNA")["SOX10",]>0)
[1] 555

after SCT
> sum(GetAssayData(object = wk4, slot = "data", assay="SCT")["TBXT",]>0)
[1] 9
> sum(GetAssayData(object = wk4, slot = "data", assay="SCT")["SHH",]>0)
[1] 26
> sum(GetAssayData(object = wk4, slot = "data", assay="SCT")["SOX10",]>0)
[1] 511


after RNA normalization/
  > sum(GetAssayData(object = normalized_7wk4, slot = "scale.data", assay="RNA")["TBXT",]>0)
[1] 10
> sum(GetAssayData(object = normalized_7wk4, slot = "scale.data", assay="RNA")["SHH",]>0)
[1] 28
> sum(GetAssayData(object = normalized_7wk4, slot = "scale.data", assay="RNA")["SOX10",]>0)





# 

#raw_seurat_list


# ncc sample using seurat norm and scaling
ncc_annotated

DefaultAssay(ncc_annotated) <- "RNA"

ncc_annotated <- NormalizeData(ncc_annotated)
ncc_annotated <- FindVariableFeatures(ncc_annotated, selection.method = "vst")
ncc_annotated <- ScaleData(ncc_annotated)

s.genes <- cc.genes$s.genes
g2m.genes <- cc.genes$g2m.genes
ncc_annotated <- CellCycleScoring(ncc_annotated, s.features = s.genes, g2m.features = g2m.genes, set.ident = TRUE)
ncc_annotated <- ScaleData(ncc_annotated, vars.to.regress = c("S.Score", "G2M.Score"), features = rownames(ncc_annotated))
ncc_annotated <- RunPCA(ncc_annotated, features = VariableFeatures(ncc_annotated))
ncc_annotated  <- FindNeighbors(ncc_annotated , dims = 1:15)
ncc_annotated  <- FindClusters(ncc_annotated , resolution = 0.3)
ncc_annotated <- RunUMAP(ncc_annotated, dims = 1:15)
DimPlot(ncc_annotated, label=TRUE)


levels(ncc_annotated@meta.data$SCT.CC.reg_snn_res.0.2) <- new.cluster.ids
# we will choose 0.3 in RNA-assay with cell-cycle regressed
# clean assays of SCT and SCT reg before proceed
library(dplyr)
library(harmony)

obj.list <- c(filtered_ncc, filtered_7wk1, filtered_7wk2, filtered_7wk3, filtered_7wk4)
merged_seurat <- merge(x = obj.list[[1]],
                       y = obj.list[2:length(obj.list)], 
                       add.cell.ids = c("NCC", "7week1", "7week2", "7week3", "7week4"), 
                       merge.data = TRUE)



merged_seurat <- merged_seurat %>%
  NormalizeData() %>%
  FindVariableFeatures(selection.method = "vst", nfeatures = 2000) %>% 
  ScaleData() 




merged_seurat <- RunPCA(merged_seurat, verbose = FALSE, npcs = 50)


###
pct <- merged_seurat[["pca"]]@stdev / sum(merged_seurat[["pca"]]@stdev) * 100

# Calculate cumulative percents for each PC
cumu <- cumsum(pct)

# Determine which PC exhibits cumulative percent greater than 90% and % variation associated with the PC as less than 5
co1 <- which(cumu > 90 & pct < 5)[1]

co1

# Determine the difference between variation of PC and subsequent PC
co2 <- sort(which((pct[1:length(pct) - 1] - pct[2:length(pct)]) > 0.1), decreasing = T)[1] + 1

# last point where change of % of variation is more than 0.1%.
co2

pcs <- min(co1, co2)

pcs

plot_df <- data.frame(pct = pct, 
                      cumu = cumu, 
                      rank = 1:length(pct))

# Elbow plot to visualize 
ggplot(plot_df, aes(cumu, pct, label = rank, color = rank > pcs)) + 
  geom_text() + 
  geom_vline(xintercept = 90, color = "grey") + 
  geom_hline(yintercept = min(pct[pct > 5]), color = "grey") +
  theme_bw()


###




# select first 14 PCs i.e. 1:! this is 04-saved object
harmonized_seurat <- RunHarmony(merged_seurat, group.by.vars = "orig.ident")
harmonized_seurat <- RunUMAP(harmonized_seurat, reduction = "harmony",dims = 1:15)



harmonized_seurat <- FindNeighbors(object = harmonized_seurat, reduction = "harmony", dims = 1:15)
harmonized_seurat <- FindClusters(harmonized_seurat, resolution = seq(0.1, 0.5, 0.1))

FeaturePlot(harmonized_seurat, c("TBXT", "SHH", "SOX10", "CTSK", "NPY", "KRT15", "PDGFRA", "APOA1"))

saveRDS(merged_seurat, file="2023-4-9-Merged_Seurat_Object_NCC.7week1-4.PCA.RDS")
saveRDS(harmonized_seurat, file="2023-4-9-HarmonyIntegrated_Seurat_Object_NCC.7week1-4.clustered.RDS")

DimPlot(harmonized_seurat, reduction = "umap", group.by = c("RNA_snn_res.0.2", "RNA_snn_res.0.4", "RNA_snn_res.0.8", "RNA_snn_res.0.9"), label=TRUE)


c4 <- subset(x = harmonized_seurat, subset = RNA_snn_res.0.2 == "4")

rownames(c4_4@meta.data[GetAssayData(object = c4_4, slot = "data", assay="RNA")["TBXT",]>0,])


tbxt_cluster_24_c4 <- RenameIdents(tbxt_cluster_24_c4,3 = 2)


refiltered_7wk1_test <- subset(filtered_7wk1, subset = nFeature_RNA > 500 & nFeature_RNA < 4000 & percent.mt < 15)

# 4-11
# 27 based on harmony on 7week samples
harmonized_seurat_7wk <- FindNeighbors(object = harmonized_seurat_7wk, reduction = "harmony", dims = 1:27)

# ~27 is a good number by ElbowPlot

hNC_cluster_11 <- subset(harmonized_seurat_7wk, RNA_snn_res.0.2 == "11")

hNC_cluster_11 <- FindNeighbors(object = hNC_cluster_11, reduction = "pca", dims = 1:27)
hNC_cluster_11 <- FindClusters(hNC_cluster_11, resolution = seq(0.1, 0.5, 0.1))

FindAllMarkers(pbmc,  min.pct = 0.25, logfc.threshold = 0.75)



# 04-17
# this is your object
filtered_7wk3 <- NormalizeData(filtered_7wk3, normalization.method = "LogNormalize")
filtered_7wk3 <- FindVariableFeatures(filtered_7wk3, selection.method = "vst", nfeatures = 2000)
all.genes.7wk3 <- rownames(filtered_7wk3)
filtered_7wk3 <- ScaleData(filtered_7wk3, features = all.genes.7wk3)
filtered_7wk3 <- RunPCA(filtered_7wk3, features = VariableFeatures(filtered_7wk3))

# determine how many PC to use for downstream analysis
pct <- filtered_7wk3[["pca"]]@stdev / sum(filtered_7wk3[["pca"]]@stdev) * 100
# Calculate cumulative percents for each PC
cumu <- cumsum(pct)
# Determine which PC exhibits cumulative percent greater than 90% and % variation associated with the PC as less than 5
co1 <- which(cumu > 90 & pct < 5)[1]
co1
co2 <- sort(which((pct[1:length(pct) - 1] - pct[2:length(pct)]) > 0.1), decreasing = T)[1] + 1
# last point where change of % of variation is more than 0.1%.
co2
pcs <- min(co1, co2)
pcs
# pcs gives 18 so I use 18 PC
# https://hbctraining.github.io/scRNA-seq_online/lessons/elbow_plot_metric.html
filtered_7wk3 <- RunUMAP(filtered_7wk3, dims = 1:18)
filtered_7wk3  <- FindNeighbors(filtered_7wk3 , dims = 1:18)
filtered_7wk3  <- FindClusters(filtered_7wk3 , resolution = seq(0.1,0.5,0.1))
DimPlot(filtered_7wk3, label=TRUE)


notochord_feature <- list(c("KRT18","KRT8","KRT19","SHH"))
filtered_7wk3 <- AddModuleScore(object = filtered_7wk3,features = notochord_feature,name = 'notochord_feature'
)
muscle_feature <- list(c("MYOD1","MYF5","PAX7"))
filtered_7wk3 <- AddModuleScore(object = filtered_7wk3,features = muscle_feature,name = 'muscle_feature'
)
chondrocyte_feature <- list(c("ACAN","SOX9","COL2A1","MIA","MATN4","HAPLN1"))
filtered_7wk3 <- AddModuleScore(object = filtered_7wk3,features = chondrocyte_feature,name = 'chondrocyte_feature'
)
conntissue_feature <- list(c("CXCL12","PRRX1","GJA1","FOXC2","TWIST1","EGFL6"))
filtered_7wk3 <- AddModuleScore(object = filtered_7wk3,features = conntissue_feature,name = 'conntissue_feature'
)
cellcycle_feature <- list(c("TOP2A","HIST1H4C","HIST1H1A"))
filtered_7wk3 <- AddModuleScore(object = filtered_7wk3,features = cellcycle_feature,name = 'cellcycle_feature'
)
neuron_feature <- list(c("TTYH1","HES6","STMN2"))
filtered_7wk3 <- AddModuleScore(object = filtered_7wk3,features = neuron_feature,name = 'neuron_feature'
)
immune_feature <- list(c("CD74","FCER1G","CSF1R"))
filtered_7wk3 <- AddModuleScore(object = filtered_7wk3,features = immune_feature,name = 'immune_feature'
)
endothelium_feature <- list(c("CDH5","EMCN","VWF"))
filtered_7wk3 <- AddModuleScore(object = filtered_7wk3,features = endothelium_feature,name = 'endothelium_feature'
)
pericyte_feature <- list(c("ACTA2","MYH11"))
filtered_7wk3 <- AddModuleScore(object = filtered_7wk3,features = pericyte_feature,name = 'pericyte_feature'
)
FeaturePlot(object = filtered_7wk3, features = c('muscle_feature1',"chondrocyte_feature1","conntissue_feature1","cellcycle_feature1","neuron_feature1","immune_feature1","endothelium_feature1"))

DimPlot(filtered_7wk3, label=TRUE, group.by = c("RNA_snn_res.0.1","RNA_snn_res.0.2","RNA_snn_res.0.3","RNA_snn_res.0.4"))
filtered_7wek3.0.3.markers <- FindAllMarkers(filtered_7wk3, min.pct = 0.25, logfc.threshold = 0.75)
filtered_7week3.0.3.markers_top50 <- filtered_7wek3.0.3.markers %>%  group_by(cluster) %>% slice_max(n = 50, order_by = avg_log2FC)
filtered_7week3.0.3.markers <- FindAllMarkers(filtered_7wk3, min.pct = 0.25, logfc.threshold = 0.75)        

write.csv(filtered_7week3.0.3.markers_top50, file="2023-04-17-filtered_7week3_resolution.0.3.markers.0.75_top50.csv")
saveRDS(filtered_7wk3, file = "2023-4-17-7week3_clustered.RDS")


# 04-18
# this is your object
filtered_7wk4 <- NormalizeData(filtered_7wk4, normalization.method = "LogNormalize")
filtered_7wk4 <- FindVariableFeatures(filtered_7wk4, selection.method = "vst", nfeatures = 2000)
all.genes.7wk4 <- rownames(filtered_7wk4)
filtered_7wk4 <- ScaleData(filtered_7wk4, features = all.genes.7wk3)
filtered_7wk4 <- RunPCA(filtered_7wk4, features = VariableFeatures(filtered_7wk4))

# determine how many PC to use for downstream analysis
pct <- filtered_7wk4[["pca"]]@stdev / sum(filtered_7wk4[["pca"]]@stdev) * 100
# Calculate cumulative percents for each PC
cumu <- cumsum(pct)
# Determine which PC exhibits cumulative percent greater than 90% and % variation associated with the PC as less than 5
co1 <- which(cumu > 90 & pct < 5)[1]
co1
co2 <- sort(which((pct[1:length(pct) - 1] - pct[2:length(pct)]) > 0.1), decreasing = T)[1] + 1
# last point where change of % of variation is more than 0.1%.
co2
pcs <- min(co1, co2)
pcs
# pcs gives 22 so I use 22 PC
# https://hbctraining.github.io/scRNA-seq_online/lessons/elbow_plot_metric.html
filtered_7wk4 <- RunUMAP(filtered_7wk4, dims = 1:22)
filtered_7wk4  <- FindNeighbors(filtered_7wk4 , dims = 1:22)
filtered_7wk4  <- FindClusters(filtered_7wk4 , resolution = seq(0.1,0.5,0.1))
DimPlot(filtered_7wk4, label=TRUE)


notochord_feature <- list(c("KRT18","KRT8","KRT19","SHH"))
filtered_7wk4 <- AddModuleScore(object = filtered_7wk4,features = notochord_feature,name = 'notochord_feature'
)
muscle_feature <- list(c("MYOD1","MYF5","PAX7"))
filtered_7wk4 <- AddModuleScore(object = filtered_7wk4,features = muscle_feature,name = 'muscle_feature'
)
chondrocyte_feature <- list(c("ACAN","SOX9","COL2A1","MIA","MATN4","HAPLN1"))
filtered_7wk4 <- AddModuleScore(object = filtered_7wk4,features = chondrocyte_feature,name = 'chondrocyte_feature'
)
conntissue_feature <- list(c("CXCL12","PRRX1","GJA1","FOXC2","TWIST1","EGFL6"))
filtered_7wk4 <- AddModuleScore(object = filtered_7wk4,features = conntissue_feature,name = 'conntissue_feature'
)
cellcycle_feature <- list(c("TOP2A","HIST1H4C","HIST1H1A"))
filtered_7wk4 <- AddModuleScore(object = filtered_7wk4,features = cellcycle_feature,name = 'cellcycle_feature'
)
neuron_feature <- list(c("TTYH1","HES6","STMN2"))
filtered_7wk4 <- AddModuleScore(object = filtered_7wk4,features = neuron_feature,name = 'neuron_feature'
)
immune_feature <- list(c("CD74","FCER1G","CSF1R"))
filtered_7wk4 <- AddModuleScore(object = filtered_7wk4,features = immune_feature,name = 'immune_feature'
)
endothelium_feature <- list(c("CDH5","EMCN","VWF"))
filtered_7wk4 <- AddModuleScore(object = filtered_7wk4,features = endothelium_feature,name = 'endothelium_feature'
)
pericyte_feature <- list(c("ACTA2","MYH11"))
filtered_7wk4 <- AddModuleScore(object = filtered_7wk4,features = pericyte_feature,name = 'pericyte_feature'
)
FeaturePlot(object = filtered_7wk4, features = c('muscle_feature1',"chondrocyte_feature1","conntissue_feature1","cellcycle_feature1","neuron_feature1","immune_feature1","endothelium_feature1","notochord_feature1"))

DimPlot(filtered_7wk4, label=TRUE, group.by = c("RNA_snn_res.0.1","RNA_snn_res.0.2","RNA_snn_res.0.3","RNA_snn_res.0.4"))
Idents(filtered_7wk4) <- "RNA_snn_res.0.3"
filtered_7week4.0.3.markers <- FindAllMarkers(filtered_7wk4, min.pct = 0.25, logfc.threshold = 0.75)
filtered_7week4.0.3.markers_top50 <- filtered_7week4.0.3.markers %>%  group_by(cluster) %>% slice_max(n = 50, order_by = avg_log2FC)
write.csv(filtered_7week4.0.3.markers_top50, file="2023-04-17-filtered_7week4_resolution.0.3.markers.0.75_top50.csv")
saveRDS(filtered_7wk4, file = "2023-4-17-7week4_clustered.RDS")


# week 1
