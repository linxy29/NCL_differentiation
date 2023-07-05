## Median scatterplot

# Subset the notochord cluster (cluster 4) from the integrated differentiated cell dataset

cluster4CFS.EU <- subset(x = CFS.EU.SCT.CC.clustered, idents = "4")

# Calculate the medians of all the genes for the iPSC and CFS cells in the integrated object
Cluster4.CFS.EU_median4 <- data.frame(row.names = row.names(Cluster4.CFS.EU@assays$integrated@data[ , Cluster4.CFS.EU$orig.ident == 'CFSEU']), median.CFS = rowMedians(as.matrix(Cluster4.CFS.EU@assays$integrated@data[ , Cluster4.CFS.EU$orig.ident == 'CFSEU'])), median.EU = rowMedians(as.matrix(Cluster4.CFS.EU@assays$integrated@data[ , Cluster4.CFS.EU$orig.ident == 'ipsc_nc_diff'])))

# data was transformed by the addition of a constant of 1.303329024 to each median value manually in excel.
# load the file
adjusted_all_genes_notochordcluster_CFSxEU <- read_excel("2023-06-06_adjusted_all_genes_notochordcluster_CFSxEU.xlsx")


# Create scatterplot of medians of all the genes for the notochord cluster
ggplot(data= adjusted_all_genes_notochordcluster_CFSxEU, mapping = aes(x = median.CFS, y = median.EU)) + 
  geom_point() + 
  stat_cor(method = "pearson", label.x = 0, label.y = 5) + 
  theme(panel.background = element_rect(fill = "white")) +
  theme(axis.line = element_line(size = 0.5, color = "black")) + 
  scale_y_continuous(limits = c(0, 5)) +
  scale_x_continuous(limits = c(0, 5))




## Matrisome genes scatterplot
# Select the genes from the integrated object for the two different cell types, and then select the commone genes with the matrisome genes (from the Naba matrisome)
integ.cfs3 <- Cluster4.CFS.EU@assays$integrated@data[, Cluster4.CFS.EU@meta.data$orig.ident=="CFSEU"]
integ.eu3 <- Cluster4.CFS.EU@assays$integrated@data[, Cluster4.CFS.EU@meta.data$orig.ident=="ipsc_nc_diff"]

matrisome_genes_cfs3 <- intersect(matrisome, rownames(Cluster4.CFS.EU@assays$integrated@data[, Cluster4.CFS.EU@meta.data$orig.ident=="CFSEU"]))
matrisome_genes_eu3 <- intersect(matrisome, rownames(Cluster4.CFS.EU@assays$integrated@data[, Cluster4.CFS.EU@meta.data$orig.ident=="ipsc_nc_diff"]))

# select the common matrisome cfs and eu matrisome genes by intersecting matrisome cfs and eu
Common_matrisome_genes_cfs_eu3 <- intersect(matrisome_genes_cfs3, matrisome_genes_eu3)

# Subset the integrated expression data to keep only the common genes for both groups
Integ_common_matrisome_genes_cfs3 <- Cluster4.CFS.EU@assays$integrated@data[Common_matrisome_genes_cfs_eu3, Cluster4.CFS.EU@meta.data$orig.ident=="CFSEU"]
Integ_common_matrisome_genes_eu3 <- Cluster4.CFS.EU@assays$integrated@data[Common_matrisome_genes_cfs_eu3, Cluster4.CFS.EU@meta.data$orig.ident=="ipsc_nc_diff"]

# Select the medians of matrisome genes
cfs.eu.integ.median_matrisome_only <- data.frame(row.names = row.names(Integ_common_matrisome_genes_cfs3), median.CFS = rowMedians(as.matrix(Integ_common_matrisome_genes_cfs3)), median.EU = rowMedians(as.matrix(Integ_common_matrisome_genes_eu3)))

# a constant of 0.849284686 was added to each median value manually in excel.
# load the file
adjusted_matrisome_notochordcluster_CFSxEU <- read_excel("2023-06-06_adjusted_matrisome_notochordcluster_CFSxEU.xlsx")

# plot the matriome gene scatterplot
ggplot(data= adjusted_matrisome_notochordcluster_CFSxEU, mapping = aes(x = median.CFS, y = median.EU)) + 
  geom_point() + 
  stat_cor(method = "pearson", label.x = 0, label.y = 5) + 
  theme(panel.background = element_rect(fill = "white")) +
  theme(axis.line = element_line(size = 0.5, color = "black"))
