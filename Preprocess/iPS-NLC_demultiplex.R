# load the libraries needed.
library(Seurat)
library(dplyr)
library(ggplot2)
library(patchwork)
set.seed(2084) # for reproducability, set a seed value

# set tha data path, cahnge to suit your location.
data_path <- c("raw_feature_bc_matrix/")
list.files(data_path) # list to check you've got the proper path
data10x <- Read10X(data.dir = data_path) # load the 10x data

data_mux <- data10x$`Multiplexing Capture` # make a separate object of the multiplexing data
dim(data_mux)
# load the gene expression data into a Seurat object with some options:
# min.cells = 100 : a gene (feature) needs to be expressed in this many cells
# min.features = 100 : a cell must have this may genes (features) expressed
seurat_object <- CreateSeuratObject(counts = data10x$`Gene Expression`, min.cells = 100, min.features = 500)
head(seurat_object@meta.data)
seurat_object

cells_in_seurat <- colnames(seurat_object) # make a vector of all the cells in the seurat object

data_mux <- data_mux[,colnames(data_mux)%in%cells_in_seurat] # subset\filter the mux data for the cells that are in the seuratobject
dim(data_mux)
data_mux <- data_mux[1:6,] # keep only the cmo301 - 306 

seurat_object[['mux']] <- CreateAssayObject(counts = data_mux) # add the mux data to the seurat object
seurat_object

# Normalize MUX data, here we use centered log-ratio (CLR) transformation
seurat_object <- NormalizeData(seurat_object, assay = "mux", normalization.method = "CLR")

# Demultiplex cells based on CMO tag enrichment
seurat_object <- MULTIseqDemux(seurat_object, assay="mux", maxiter = 8, autoThresh = TRUE, verbose = TRUE) # do a quantile assignment of singlets, doublets and negative cells from multiplexing experiments. Annotate singlets by tags.
# have a look at the results.
table(seurat_object$MULTI_ID)
table(seurat_object$MULTI_classification)
head(seurat_object@meta.data)

# add additional metadata based on the CMO tags
seurat_object$sample <- ifelse(seurat_object$MULTI_ID=="CMO301", "D7NT_ThermoFiPSC", 
                               ifelse(seurat_object$MULTI_ID=="CMO302", "D7NOTO_ThermoFiPSC",
                                      ifelse(seurat_object$MULTI_ID=="CMO303", "D7NOTOSB_ThermoFiPSC",
                                             ifelse(seurat_object$MULTI_ID=="CMO304", "D7NT_190PBMC4FiPSC",
                                                    ifelse(seurat_object$MULTI_ID=="CMO305", "D7NOTO_190PBMC4FiPSC",
                                                           ifelse(seurat_object$MULTI_ID=="CMO306", "D7NOTOSB_190PBMC4FiPSC", seurat_object$MULTI_ID))))))

seurat_object$sample <- ifelse(seurat_object$sample==7, "Doublet", ifelse(seurat_object$sample==8, "Negative", seurat_object$sample))
head(seurat_object@meta.data)
table(seurat_object$sample)

seurat_object$treatment <- ifelse(seurat_object$MULTI_ID=="CMO301", "NT", 
                                  ifelse(seurat_object$MULTI_ID=="CMO302", "NOTO",
                                         ifelse(seurat_object$MULTI_ID=="CMO303", "NOTOSB",
                                                ifelse(seurat_object$MULTI_ID=="CMO304", "NT",
                                                       ifelse(seurat_object$MULTI_ID=="CMO305", "NOTO",
                                                              ifelse(seurat_object$MULTI_ID=="CMO306", "NOTOSB", "None"))))))

head(seurat_object@meta.data)
table(seurat_object$treatment)

seurat_object$cellline <- ifelse(seurat_object$MULTI_ID=="CMO301", "ThermoFiPSC", 
                                 ifelse(seurat_object$MULTI_ID=="CMO302", "ThermoFiPSC",
                                        ifelse(seurat_object$MULTI_ID=="CMO303", "ThermoFiPSC",
                                               ifelse(seurat_object$MULTI_ID=="CMO304", "PBMC4FiPSC",
                                                      ifelse(seurat_object$MULTI_ID=="CMO305", "PBMC4FiPSC",
                                                             ifelse(seurat_object$MULTI_ID=="CMO306", "PBMC4FiPSC", "None"))))))
head(seurat_object@meta.data)
table(seurat_object$cellline)

# Save the resulting seurat object as a RDS data file for further processing and analysis.
# the line has been commented to prevent overwriting it by accident. 
# Uncomment (remove #) and probably change the name of the file to be written if you want to save a new file
saveRDS(seurat_object, "ipsc_nc_diff_multi_annotated.RDS")

# the other way in Seurat to demultiplex the cells based on their CMO tag is by HTODemux function.
# below some code to do that and to compare it to the MULTIseDemux method used previously.
Idents(seurat_object) <- "orig.indent"
# do the demultiplexing with the HTODemux method
seurat_object <- HTODemux(seurat_object, assay = "mux", positive.quantile = 0.99)
head(seurat_object@meta.data)

table(seurat_object$MULTI_ID) # table of cells per tag MULTIseDemux
table(seurat_object$hash.ID) # table of cells per tag HTODemux
table(seurat_object$MULTI_ID, seurat_object$hash.ID) # table of cells per tag MULTIseDemux vs HTODemux
table(seurat_object$MULTI_ID==seurat_object$hash.ID) # how many cells did receive the same tag? (sum of the diagonal of previous table)

# set the identity of each cell to the tag fount with MULTIseDemux ("MULTI_ID")
Idents(seurat_object) <- "MULTI_ID"
# make a ridge plot for counts of each CMO tag
rp <- RidgePlot(seurat_object, assay = "mux", features = rownames(seurat_object[["mux"]])[1:6], combine = FALSE)
rp <- lapply(X = rp, FUN = function(x) x + theme(plot.title = element_text(size = 8), axis.text=element_text(size=8),axis.title=element_text(size=8)))
wrap_plots(rp)+plot_layout(guides = 'collect')+plot_annotation(title = "Idents=MULTI_ID")

# set the identity of each cell to the tag fount with HTODemux ("hash.ID")
Idents(seurat_object) <- "hash.ID"
# make a ridge plot for counts of each CMO tag
rp <- RidgePlot(seurat_object, assay = "mux", features = rownames(seurat_object[["mux"]])[1:6], combine = FALSE)
rp <- lapply(X = rp, FUN = function(x) x + theme(plot.title = element_text(size = 8), axis.text=element_text(size=8),axis.title=element_text(size=8)))
wrap_plots(rp)+plot_layout(guides = 'collect')+plot_annotation(title = "Idents=hash.ID")


# Compare number of UMIs for singlets, doublets and negative cells
Idents(seurat_object) <- "MULTI_ID"
VlnPlot(seurat_object, features = c("nCount_RNA", "nFeature_RNA", "nCount_mux"), pt.size = 0.1, log = TRUE)

Idents(seurat_object) <- "hash.ID"
VlnPlot(seurat_object, features = c("nCount_RNA", "nFeature_RNA", "nCount_mux"), pt.size = 0.1, log = TRUE)

