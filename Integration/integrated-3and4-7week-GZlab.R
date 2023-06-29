

obj.list <- c(filtered_7wk3, filtered_7wk4)
merged_seurat <- merge(x = obj.list[[1]],
                       y = obj.list[2:length(obj.list)], 
                       add.cell.ids = c("7week3", "7week4"), 
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

harmonized_seurat <- RunHarmony(merged_seurat, group.by.vars = "orig.ident")

harmonized_seurat <- RunUMAP(harmonized_seurat, reduction = "harmony", dims = 1:20)

harmonized_seurat <- FindNeighbors(object = harmonized_seurat, reduction = "harmony", dims = 1:20)

###
notochord_feature <- list(c("KRT18","KRT8","KRT19","SHH"))
harmonized_seurat <- AddModuleScore(object = harmonized_seurat,features = notochord_feature,name = 'notochord_feature'
)
muscle_feature <- list(c("MYOD1","MYF5","PAX7"))
harmonized_seurat <- AddModuleScore(object = harmonized_seurat,features = muscle_feature,name = 'muscle_feature'
)
chondrocyte_feature <- list(c("ACAN","SOX9","COL2A1","MIA","MATN4","HAPLN1"))
harmonized_seurat <- AddModuleScore(object = harmonized_seurat,features = chondrocyte_feature,name = 'chondrocyte_feature'
)
conntissue_feature <- list(c("CXCL12","PRRX1","GJA1","FOXC2","TWIST1","EGFL6"))
harmonized_seurat <- AddModuleScore(object = harmonized_seurat,features = conntissue_feature,name = 'conntissue_feature'
)
cellcycle_feature <- list(c("TOP2A","HIST1H4C","HIST1H1A"))
harmonized_seurat <- AddModuleScore(object = harmonized_seurat,features = cellcycle_feature,name = 'cellcycle_feature'
)
neuron_feature <- list(c("TTYH1","HES6","STMN2"))
harmonized_seurat <- AddModuleScore(object = harmonized_seurat,features = neuron_feature,name = 'neuron_feature'
)
immune_feature <- list(c("CD74","FCER1G","CSF1R"))
harmonized_seurat <- AddModuleScore(object = harmonized_seurat,features = immune_feature,name = 'immune_feature'
)
endothelium_feature <- list(c("CDH5","EMCN","VWF"))
harmonized_seurat <- AddModuleScore(object = harmonized_seurat,features = endothelium_feature,name = 'endothelium_feature'
)
pericyte_feature <- list(c("ACTA2","MYH11"))
harmonized_seurat <- AddModuleScore(object = harmonized_seurat,features = pericyte_feature,name = 'pericyte_feature'
)
FeaturePlot(object = harmonized_seurat, features = c('muscle_feature1',"chondrocyte_feature1","conntissue_feature1","cellcycle_feature1","neuron_feature1","immune_feature1","endothelium_feature1", "notochord_feature1"))


harmonized_seurat <- RenameIdents(object = harmonized_seurat,  
                                  '3' = 'hChon',  '4' = 'hChon',
                                  '6' = 'hMuscle', 
                                  '2'='hCT','5'='hCT','0'='hCT',
                                  '10'='hCycle', '1'='hCycle',
                                  '8'='hNeuron','7'='hNeuron','9'='hNeuron','11'='hNeuron',
                                  '12'='hImmune','13'='hEC','14'='hNC')