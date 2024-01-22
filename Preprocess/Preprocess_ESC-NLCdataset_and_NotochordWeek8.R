## Processing of 8 week fetal notochord data

#Load data into Seurat.
NCC.raw.data <- Read10X_h5("~/R/data/2022-07-26_ncc_filtered_feature_bc_matrix.h5", use.names = TRUE, unique.features = TRUE)


# Make Seurat object
NCC.CR.6.1.2 <- CreateSeuratObject(counts = NCC.raw.data, project = "NCC", min.cells = 3, min.features = 200)

## Filter the cells

# Calculate the proportion of transcripts mapping to mitochondrial genes and add them to the metadata
NCC.CR.6.1.2[["percent.mt"]] <- PercentageFeatureSet(NCC.CR.6.1.2, pattern = "^MT-")

# Visualize QC metrics as a violin plot
VlnPlot(NCC.CR.6.1.2, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)

# Filter the cells that have unique feature counts over 200 and less than 5000, and also filter cells that have >10% mitochondrial counts

NCC.filtered.5000nF <- subset(NCC.CR.6.1.2, subset = nFeature_RNA > 200 & nFeature_RNA < 5000 & percent.mt < 10)

# Normalize data with SCTransform and perform cell cycle regression
NCC.SCT <- SCTransform(NCC.filtered.5000nF, verbose = FALSE, do.scale = TRUE, return.only.var.genes = FALSE)

NCC.SCT <- CellCycleScoring(NCC.SCT, s.features = s.genes, g2m.features = g2m.genes, assay = 'SCT', set.ident = TRUE)

NCC.SCT <- SCTransform(NCC.SCT, assay = 'RNA', new.assay.name = 'SCT.CC.reg' , vars.to.regress = c('S.Score', 'G2M.Score'), verbose = FALSE, do.scale = TRUE, return.only.var.genes = FALSE)

#Save the file
saveRDS(NCC.SCT, '2022-10-13_NCC.SCT.CC.RDS')


# run PCA
NCC.SCT.CC <- RunPCA(NCC.SCT.CC, features = VariableFeatures(object = NCC.SCT.CC))


# find neighbours
NCC.SCT.CC <- FindNeighbors(NCC.SCT.CC, dims = 1:10)

# find clusters
NCC.SCT.CC <- FindClusters(NCC.SCT.CC, resolution = c(0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0, 1.2, 1.4, 1.6))

# Run umap
NCC.SCT.CC <- RunUMAP(NCC.SCT.CC, dims = 1:10)

# Set resolution 0.2
Idents(object = NCC.SCT.CC) <- "SCT.CC.reg_snn_res.0.2"




## Processing of the H1-CFS differentiated cell data

# Load data into Seurat.
CFSEUdiff.data <- Read10X_h5("~/20221123CFSEU/outs/filtered_feature_bc_matrix.h5", use.names = TRUE, unique.features = TRUE)

# Make Seurat object
CFSEUdiff <- CreateSeuratObject(counts = CFSEUdiff.data, project = "CFSEU", min.cells = 3, min.features = 200)

## Filter the cells
# Calculate the proportion of transcripts mapping to mitochondrial genes and add them to the metadata
CFSEUdiff[["percent.mt"]] <- PercentageFeatureSet(CFSEUdiff, pattern = "^MT-")

# Visualize QC metrics as a violin plot
VlnPlot(CFSEUdiff, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)

# Filter the cells that have unique feature counts over 200 and less than 7500, and also filter cells that have >10% mitochondrial counts
CFSEUdiff.filtered <- subset(CFSEUdiff, subset = nFeature_RNA > 200 & nFeature_RNA < 7500 & percent.mt < 10)

# Normalize data with SCTransform and perform cell cycle regression
CFSEUdiff.SCT <- SCTransform(CFSEUdiff.filtered, verbose = FALSE, do.scale = TRUE, return.only.var.genes = FALSE)

CFSEUdiff.SCT <- CellCycleScoring(CFSEUdiff.SCT, s.features = s.genes, g2m.features = g2m.genes, assay = 'SCT', set.ident = TRUE)

CFSEUdiff.SCT <- SCTransform(CFSEUdiff.SCT, assay = 'RNA', new.assay.name = 'SCTCCreg' , vars.to.regress = c('S.Score', 'G2M.Score'), verbose = FALSE, do.scale = TRUE, return.only.var.genes = FALSE)

#Save the file
saveRDS(CFSEUdiff.SCT, '2022-11-24_CFSEUdiff.SCT.RDS')

# Run PCA
CFSEUdiff.SCT <- RunPCA(CFSEUdiff.SCT, verbose = FALSE)

# Run UMAP and find neighbours
CFSEUdiff.SCT <- RunUMAP(CFSEUdiff.SCT, dims = 1:30, verbose = FALSE)

CFSEUdiff.SCT <- FindNeighbors(CFSEUdiff.SCT, dims = 1:30, verbose = FALSE)

#Find clusters at different resolutions 
CFSEUdiff.SCT <- FindClusters(CFSEUdiff.SCT, resolution = c(0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0, 1.2))

#choose resolution 0.4
Idents(object = CFSEUdiff.SCT) <- "SCTCCreg_snn_res.0.4"
