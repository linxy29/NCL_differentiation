library(Seurat)
library(tidyverse)
library(RColorBrewer)
library(patchwork)
library(kableExtra)
source("functions01.R")
set.seed(2084) # for reproducibility set a seed

# place holder for prerocessing steps
# currently is is not used (eval=FALSE) so it contains the table of conditions and tags
D7NT_ThermoFiPSC,CMO301
D7NOTO_ThermoFiPSC,CMO302
D7NOTOSB_ThermoFiPSC,CMO303
D7NT_190PBMC4FiPSC,CMO304
D7NOTO_190PBMC4FiPSC,CMO305
D7NOTOSB_190PBMC4FiPSC,CMO306

## Seurat analysis {.tabset}

# set the path to reflect the location of your data
data_path <- "../data/"
seuratRaw <- readRDS(paste0(data_path,"ipsc_nc_diff_multi_annotated.RDS"), refhook = NULL)
#head(seuratRaw@meta.data)
seuratRaw$orig.ident <- c("ipsc_nc_diff")
seuratRaw[["percent.mt"]] <- PercentageFeatureSet(seuratRaw, pattern = "^MT-")

### Plot raw data
# define some filter values 
my_nFeature_RNA_min <- 1000
my_nFeature_RNA_max <- 8000
my_nCount_RNA_min <- 1000
my_nCount_RNA_max <- 50000
my_percent_MT <- 10

# plot the raw data
vp <- VlnPlot(seuratRaw, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), group.by = "MULTI_ID", ncol = 3, pt.size = 0.5, cols = color_labels, combine = FALSE)
vp <- lapply(X = vp, FUN = function(x) x + theme(
  plot.title = element_text(size = 8),
  axis.text=element_text(size=6),
  axis.title=element_text(size=8)))
vp[[1]] <- vp[[1]] + geom_hline(yintercept = c(my_nFeature_RNA_min, my_nFeature_RNA_max), colour = "#990000", linetype = "dashed", size=0.75)
vp[[2]] <- vp[[2]] + geom_hline(yintercept = c(my_nCount_RNA_min, my_nCount_RNA_max), colour = "#990000", linetype = "dashed", size=0.75)
vp[[3]] <- vp[[3]] + geom_hline(yintercept = c(my_percent_MT), colour = "#990000", linetype = "dashed", size=0.75)

wrap_plots(vp, ncol = 3) & NoLegend()

### scatterplot counts and features
FeatureScatter(seuratRaw, "nCount_RNA", "nFeature_RNA", group.by = "MULTI_ID", pt.size = 0.5, cols = color_labels)+ geom_vline(xintercept = c(my_nCount_RNA_min,my_nCount_RNA_max), colour = "#990000", linetype = "dashed") + geom_hline(yintercept = c(my_nFeature_RNA_min,my_nFeature_RNA_max), colour = "#990000", linetype = "dashed") + guides(colour = guide_legend(override.aes = list(size=5)))

### scatterplot counts and percent MT
FeatureScatter(seuratRaw, feature1 = "nCount_RNA", feature2 = "percent.mt", pt.size = 0.5, cols = color_labels)+ geom_vline(xintercept = c(my_nCount_RNA_min, my_nCount_RNA_max), colour = "#990000", linetype = "dashed") + geom_hline(yintercept = c(my_percent_MT), colour = "#990000", linetype = "dashed") + guides(colour = guide_legend(override.aes = list(size=5)))

### scatterplot features and percent MT
FeatureScatter(seuratRaw, feature1 = "nFeature_RNA", feature2 = "percent.mt", pt.size = 0.5, cols = color_labels)+ geom_vline(xintercept = c(my_nFeature_RNA_min,my_nFeature_RNA_max), colour = "#990000", linetype = "dashed") + geom_hline(yintercept = c(my_percent_MT), colour = "#990000", linetype = "dashed") + guides(colour = guide_legend(override.aes = list(size=5)))

## Filtering {.tabset}
# make a new seurat object with filtered data.
ips2nc_fltrd <- subset(seuratRaw, subset = nFeature_RNA > my_nFeature_RNA_min & nFeature_RNA < my_nFeature_RNA_max &
                          nCount_RNA > my_nCount_RNA_min & nCount_RNA < my_nCount_RNA_max &
                          percent.mt < my_percent_MT)

tbl_01 <- as.data.frame(table(ips2nc_fltrd$MULTI_ID))
colnames(tbl_01) <- c("Tag","Cell number")


# also remove the Doublets and Negative tagged cells
ips2nc_fltrd <- subset(ips2nc_fltrd, subset = MULTI_ID == "Negative" | MULTI_ID == "Doublet", invert = TRUE)

tbl_02 <- as.data.frame(table(ips2nc_fltrd$MULTI_ID))
colnames(tbl_02) <- c("Tag","Cell number")

# drop the unused levels of Doublet and Negative
ips2nc_fltrd@meta.data$MULTI_ID <- droplevels(ips2nc_fltrd@meta.data$MULTI_ID)

tbl_03 <- as.data.frame(table(ips2nc_fltrd$MULTI_ID))
colnames(tbl_03) <- c("Tag","Cell number")

tbl_04 <- as.data.frame(ips2nc_fltrd@meta.data %>% group_by(sample) %>% summarise(sum(nCount_RNA),sum(nFeature_RNA)))

Some filtering is needed. 

Filters used:

* Minimum number of features(genes) `r my_nFeature_RNA_min` and maximum number of features(genes) `r format(my_nFeature_RNA_max, scientific = FALSE)`. 
* Counts between `r my_nCount_RNA_min`  and `r format(my_nCount_RNA_max, scientific = FALSE)`. 
* The maximum percentage of mitochondrial counts was set to `r my_percent_MT` .


Data still contains all the tags from the demultiplexing step including the doublet and Negative cells:

`r print(knitr::kable(tbl_01, caption = "Number of cells per tag after filtering") %>% kable_styling())`

After removing these cells the makeup is

`r print(knitr::kable(tbl_03, caption = "Number of cells per tag after filtering and removal of Doublet and Negative tags") %>% kable_styling())`

After this the data looked like this:
`r print(ips2nc_fltrd)`

The total counts and features per sample:

`r print(knitr::kable(tbl_04, caption = "Total counts and features per sample after filtering and removal of Doublet and Negative tags") %>% kable_styling())`

### Plot filtered data
```{r plot_filt01, echo=FALSE, message=FALSE, warning=FALSE, fig.width=9, fig.cap="Violin plot of the filterd data that will be used for further analysis."}
vp <- VlnPlot(ips2nc_fltrd, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), group.by = "MULTI_ID", ncol = 3, combine = FALSE, cols = color_labels)
vp <- lapply(X = vp, FUN = function(x) x + theme(
  plot.title = element_text(size = 8),
  axis.text=element_text(size=6),
  axis.title=element_text(size=8)))

wrap_plots(vp, ncol = 3) & NoLegend()
```

### scatterplot counts and features
```{r plot_filterd_02, echo=FALSE, message=FALSE, warning=FALSE, fig.width=9, fig.cap="Relation between counts and features."}
FeatureScatter(ips2nc_fltrd, "nCount_RNA", "nFeature_RNA",  group.by = "MULTI_ID", pt.size = 0.5, cols = color_labels) + guides(colour = guide_legend(override.aes = list(size=5)))
```

## Clustering {.tabset}

Before further analysis the raw counts are normalized and variance stabilized using the SCTransform method.
This method uses a regularized binomial model to model the UMI counts and remove variation due to sequencing depth i.e. total n UMIs per cell. At the same time adjusting the variance based on pooling expressions across genes with similar expressions. The model returns residuals that are the normalized expression levels for the tested genes.  
The scTransform method regresses out the sequencing depth automatically. 

**References: **

1. Hafemeister, C., Satija, R. Normalization and variance stabilization of single-cell RNA-seq data using regularized negative binomial regression. Genome Biol 20, 296 (2019). <https://doi.org/10.1186/s13059-019-1874-1>

2. <https://satijalab.org/seurat/articles/sctransform_vignette.html>


After This normalization step the dimensionality of the data is reduced using principal component analysis (**PCA**) and the standard deviations of each of the principle components (PC) are plotted in an elbow plot. This indicates the most significant principal components that can be easily identified by the "elbow" bend in the graph. This number of principal components is subsequently used in the K-nearest neighbor (**KNN**) embedding to generate the cell clusters. Lastly, the clusters are visualized with t-distributed stochastic neighbor embedding (**t-SNE**).

```{r clustering01, echo=FALSE, message=FALSE, warning=FALSE, fig.width=9, fig.cap="Elbow plot plotting the standard deviations for each of the principal components. The red triangle indicates the number of (informative) PCs used in the downstream analysis."}
ips2nc_sctrans <- SCTransform(ips2nc_fltrd, verbose = FALSE, do.scale = TRUE, return.only.var.genes = FALSE)#, vars.to.regress = "percent.mt")# default settings uses no scaling so set do.scale = TRUE, default return only variable genes set return.only.var.genes = FALSE to get all genes into the scale.data)

# perform dimentionality reduction
ips2nc_sctrans <- RunPCA(ips2nc_sctrans, verbose = FALSE)
#DimPlot(ips2nc_sctrans, reduction = "pca", group.by = "tissue_source")

# Determine percent of variation associated with each PC
pct <- ips2nc_sctrans@reductions$pca@stdev / sum(ips2nc_sctrans@reductions$pca@stdev) * 100
# Calculate cumulative percents for each PC
cum <- cumsum(pct)
# Determine which PC exhibits cumulative percent greater than 90% and % variation associated with the PC as less than 5
co1 <- which(cum > 90 & pct < 5)[1]
# Determine the difference between variation of PC and subsequent PC
co2 <- sort(which((pct[1:length(pct)-1] - pct[2:length(pct)]) > 0.1),  decreasing = T)[1] + 1 # last point where change of % of variation is more than 0.1%.
#co2
# Minimum of the two calculation
pcs <- min(co1, co2) # change to any other number
#pcs

ElbowPlot(ips2nc_sctrans)+ geom_point(aes(pcs,ips2nc_sctrans@reductions$pca@stdev[pcs]), color="red", size=5, shape=6)

## clustering using the first pcs PCs (ori: resolutio = 0.4)
ips2nc_sctrans <- FindNeighbors(ips2nc_sctrans, dims = 1:pcs)
ips2nc_sctrans <- FindClusters(object = ips2nc_sctrans, dims = 1:pcs, verbose = TRUE, save.SNN = TRUE, resolution = 0.4)

# using the first 10 dims for umap and tsne
ips2nc_sctrans <- RunUMAP(ips2nc_sctrans, dims = 1:pcs, verbose = FALSE)
ips2nc_sctrans <- RunTSNE(ips2nc_sctrans, dims = 1:pcs, verbose = FALSE)

### Table of cells per cluster
tbl_05 <- as.data.frame(table(ips2nc_sctrans$seurat_clusters))
tbl_05 <- t(tbl_05)
tbl_05 <- as.data.frame(tbl_05)

rownames(tbl_05) <- c("Cluster","Cell number")
colnames(tbl_05) <- NULL
knitr::kable(tbl_05, caption = "Number of cells per cluster") %>% kable_styling()

### StackedBar
sample_cluster <- as.data.frame(table(ips2nc_sctrans$seurat_clusters, ips2nc_sctrans$sample))
colnames(sample_cluster) <- c("cluster", "group", "Freq")
color_labels_2 <- c("#80B1D3", "#999999", "#FDB462", "#BEBADA", "#FB8072", "#8DD3C7")


ggplot(sample_cluster, aes(x=cluster, y=Freq, fill=group)) + geom_bar(stat = 'identity')+ scale_fill_manual(values = color_labels_2)+ NoGrid() + WhiteBackground() + theme(axis.line = element_line(size = 0.5, linetype = "solid", colour = "black"))

## Cluster visualization, tSNE or umap plots. {.tabset}

### Basic UMAP
reduction_used <- "umap" # specify what reduction to use (tsne or umap) for all plots 
DimPlot(object = ips2nc_sctrans, pt.size = 0.5, reduction = reduction_used, label = TRUE)


### Basic tSNE
reduction_used <- "tsne" # specify what reduction to use (tsne or umap) for all plots
DimPlot(object = ips2nc_sctrans, pt.size = 0.5, reduction = reduction_used, label = TRUE)


### Tag
base_plot <- DimPlot(ips2nc_sctrans, reduction = reduction_used, label = TRUE, pt.size = 0.2)+NoLegend()

plot_tag <- DimPlot(object = ips2nc_sctrans, group.by = 'MULTI_ID', pt.size = 0.5, reduction = reduction_used, label = FALSE, cols = color_labels)

plot_tag + inset_element(base_plot, left = 0.8, bottom = 0.6, right = 1.2, top = 1.1)


### Sample
plot_sample <- DimPlot(object = ips2nc_sctrans, group.by = 'sample', pt.size = 0.5, reduction = reduction_used, label = FALSE, cols = color_labels_2)
plot_sample + inset_element(base_plot, left = 0.8, bottom = 0.6, right = 1.2, top = 1.1)


### Treat
base_plot <- DimPlot(ips2nc_sctrans, reduction = reduction_used, label = TRUE, pt.size = 0.2)+NoLegend()

plot_treatment <- DimPlot(object = ips2nc_sctrans, group.by = 'treatment', ncol = 1 , pt.size = 0.5, reduction = reduction_used, label = FALSE, cols = color_labels)

plot_treatment + inset_element(base_plot, left = 0.8, bottom = 0.6, right = 1.2, top = 1.1)


### Cellline
base_plot <- DimPlot(ips2nc_sctrans, reduction = reduction_used, label = TRUE, pt.size = 0.2)+NoLegend()
plot_cellline <- DimPlot(object = ips2nc_sctrans, group.by = 'cellline', ncol = 1 , pt.size = 0.5, reduction = reduction_used, label = FALSE, cols = color_labels)

plot_cellline + inset_element(base_plot, left = 0.8, bottom = 0.6, right = 1.2, top = 1.1)
