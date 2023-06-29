## NCC and NLCs heatmap

# Load the object
NCC.SCT.CC <- readRDS("~/R/data/2022-10-13_NCC.SCT.CC.RDS")

# run PCA
NCC.SCT.CC <- RunPCA(NCC.SCT.CC, features = VariableFeatures(object = NCC.SCT.CC))

# Clustering the cells
# find neighbours
NCC.SCT.CC <- FindNeighbors(NCC.SCT.CC, dims = 1:10)

# find clusters
NCC.SCT.CC <- FindClusters(NCC.SCT.CC, resolution = c(0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0, 1.2, 1.4, 1.6))

# Run umap
NCC.SCT.CC <- RunUMAP(NCC.SCT.CC, dims = 1:10)

# Set resolution 0.2
Idents(object = NCC.SCT.CC) <- "SCT.CC.reg_snn_res.0.2"


# identify all the marker genes
NCC.0.2.allgenes <- FindAllMarkers(NCC.SCT.CC, only.pos = FALSE, min.pct = 0.25, logfc.threshold = 0.25)

# This list was  matched to the matrisome genes by Naba et al, and exported as an excel file "NCC.matrisome"

# Load the integrated object for iPS and CFS NLCs
CFS.EU.SCT.CC.clustered <- readRDS("~/R/data/2023-06-08_CFSxEU_clustered_annotated.RDS")

# set resolution at 0.3 for the integrated object
Idents(CFS.EU.SCT.CC.clustered) = CFS.EU.SCT.CC.clustered$integrated_snn_res.0.3

# find all markers
CFS.EU.SCT.0.3.all.markers <-FindAllMarkers(CFS.EU.SCT.CC.clustered, only.pos=FALSE, min.pct = 0.1, logfc.threshold = 0.25)

# This list was also matched to the matrisome genes by Naba et al, and exported as an excel file "NLC.matrisome"
 
# The common matrisome genes between NCC (notochord cluster) and NLC (notochord cluster) lists were selected and saved as common_matrisome_NCC_CFS_EU in excel
# Load common matrisome list

common_matrisome_NCC_CFS_EU <- read_excel("2023-04-25_common_matrisome_NCC_CFS.EU.xlsx")

# Subset cluster 4 (notochord cluster) from NCC object
cluster4NCC <- subset(x = NCC.SCT.CC, idents = "4")

# Plot heatmap of the matrisome genes for notochord cluster in NCC
DoHeatmap(cluster4NCC, features = common_matrisome_NCC_CFS_EU$Common_matrisome) + theme(axis.text.y = element_text(size = 7))

# subset the somitic mesoderm cluster in NCC
SomiticMesoderm <- subset(x = NCC.SCT.CC, idents = "1")

#  do heatmap for common genes in somitic mesoderm cluster in NCC
DoHeatmap(SomiticMesoderm, features = common_matrisome_NCC_CFS_EU$Common_matrisome) + theme(axis.text.y = element_text(size = 7))



# Subset cluster 4 (notochord cluster) from integrated obcject
cluster4CFS.EU <- subset(x = CFS.EU.SCT.CC.clustered, idents = "4")

# Plot heatmap of the matrisome genes for notochord cluster of integrated object
DoHeatmap(cluster4CFS.EU, features = common_matrisome_NCC_CFS_EU$Common_matrisome, group.by = "orig.ident") + theme(axis.text.y = element_text(size = 7))

# Subset cluster 7 (endoderm) from integrated object 
Cluster7CFS.EU <- subset(x = CFS.EU.SCT.CC.clustered, idents = "7")

# Plot heatmap for endoderm of integrated object
DoHeatmap(Cluster7CFS.EU, features = common_matrisome_NCC_CFS_EU$Common_matrisome, group.by = "orig.ident") + theme(axis.text.y = element_text(size = 7))


