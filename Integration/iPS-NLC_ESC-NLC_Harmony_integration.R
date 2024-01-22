
merged_objects <- readRDS("2023-03-27-Merged_CFS.EU.SCT.CC.only.unintegrated.RDS")
# split first
# Normal Integration followed from SCTransform
# https://satijalab.org/seurat/articles/sctransform_v2_vignette.html Perform integration using Pearson residuals

features <- SelectIntegrationFeatures(object.list = merged_objects, nfeatures = 3000)
integration.list <- PrepSCTIntegration(object.list = merged_objects, anchor.features = features)

anchors <- FindIntegrationAnchors(object.list = integration.list, normalization.method = "SCT",
                                         anchor.features = features)anchors <- FindIntegrationAnchors(object.list = integration.list, normalization.method = "SCT",
                                         anchor.features = features)

combined.sct <- IntegrateData(anchorset = anchors, normalization.method = "SCT")
# "integrated" assay is now by default
combined.sct <- RunPCA(combined.sct, verbose = FALSE)
combined.sct <- RunUMAP(combined.sct, reduction = "pca", dims = 1:15)

combined.sct <- FindNeighbors(combined.sct, reduction = "pca", dims = 1:15)
combined.sct <- FindClusters(combined.sct, resolution = c(seq(0.1,1,0.1)))

# pca
DimPlot(combined.sct, group.by = 'orig.ident', pt.size = 0.1, reduction="pca")
# umap
DimPlot(combined.sct, group.by = 'orig.ident', pt.size = 0.1, reduction="umap")

clustree(combined.sct@meta.data, prefix="integrated_snn_res.")

### Harmony 
# ! 2023-03-29 A bug in harmony need to fix median_umi in SCTModel - lost during integration
# https://hbctraining.github.io/scRNA-seq_online/lessons/06a_integration_harmony.html
# Calculate PCs using manually set variable features
merged_seurat <- RunPCA(merged_objects, assay = "SCT", npcs = 50, features = VariableFeatures(merged_objects))

# ElbowPlot suggest 30 is a good number
ElbowPlot(merged_seurat, ndims = 50)

# harmony integration
DimPlot(merged_seurat, group.by = 'orig.ident', pt.size = 0.1)
merged_seurat <- RunUMAP(merged_seurat, reduction = "harmony", assay = "SCT", dims = 1:15)
hm.integrated <- RunHarmony(merged_seurat, 
                                group.by.vars = "orig.ident", 
                                reduction = "pca", assay.use = "SCT", reduction.save = "harmony", dims= 1:15)

hm.integrated <- FindNeighbors(object = hm.integrated, reduction = "harmony", graph.name = "integrated")

hm.integrated <- FindClusters(graph.name = "integrated", 
  object = hm.integrated,
  resolution = c(seq(0.1,0.8,0.1))
) 

# maybe 0.2 or 0.3?
clustree(hm.integrated@meta.data, prefix = "integrated_res.")

# pca
DimPlot(hm.integrated, group.by = 'orig.ident', pt.size = 0.1, reduction="harmony")
# umap
DimPlot(hm.integrated, group.by = 'orig.ident', pt.size = 0.1, reduction="umap")

# use 0.3??
DimPlot(hm.integrated, group.by = 'SCT_snn_res.0.3', pt.size = 0.1, reduction="umap")

# 
clustree(hm.integrated@meta.data, prefix = "SCT_snn_res.")


# Pure merge without integration
srt.merged <- RunPCA(merged_objects, assay = "SCT", npcs = 50, features = VariableFeatures(merged_objects))

srt.merged <- RunUMAP(srt.merged, dims = 1:15, reduction = "pca", assay = "SCT")
srt.merged <- FindNeighbors(object = srt.merged, dims = 1:15, reduction="pca", graph.name = "merged")
srt.merged <- FindClusters(graph.name="merged",
  object = srt.merged,
  resolution = c(seq(0.1,0.8,0.1))
) 

DimPlot(srt.merged, group.by = 'orig.ident', pt.size = 0.1, reduction="pca")

DimPlot(srt.merged, group.by = 'orig.ident', pt.size = 0.1, reduction="umap")

DimPlot(srt.merged, group.by = 'SCT_snn_res.0.3', pt.size = 0.1, label=TRUE, split.by = "orig.ident")

clustree(srt.merged@meta.data, prefix = "merged_res.")
