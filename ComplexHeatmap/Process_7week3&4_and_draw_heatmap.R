## ---------------------------
##
## Script name: Process_7week3_and_draw_heatmap
##
## Purpose of script:
## Process datastets 7week3 to make a heatmap with WGCNA brown and cyan modules
##
## Author: Bluwen Guidoux d'Halluin
##
## Date Created: 2023-06-28
##
## ---------------------------
##
## Notes:
## R version 4.2.2
## ComplexHeatmap 2.15.4
## Seurat 4.2.0
## dplyr 1.1.1
## readxl 1.4.2
## circlize 0.4.15
## readr 2.1.4
##
## ---------------------------
##
## Package

library(ComplexHeatmap)
library(Seurat)
library(dplyr)
library(readxl)
library(circlize)
library(readr)

## ---------------------------

## Importing Dataset

S3 <-readRDS("data/s3_annotated.rds")

## ---------------------------

## Process data

S3 <- subset(S3, subset = cluster.origin %in% c("hNC_Axial_skeleton_Week7 ", "hCT_Axial_skeleton_Week7 "))

## Scale data

all_genes <- rownames(S3)
S3 <- ScaleData(S3, features = all_genes)

## ---------------------------

## Prepare annotations, order and colours for the heatmap

meta_data <- S3@meta.data

ordered_meta_data <- meta_data[order(meta_data$orig.ident), ]

ordered_meta_data <- ordered_meta_data[, c("data.origin",
                                           "type.origin",
                                           "cluster.origin")]

ordered_meta_data <- ordered_meta_data[order(ordered_meta_data$data.origin), ]

annotation_colors <- list("data.origin" = c("Zhou" = "#6670b5"),
                          "type.origin" = c("Axial_skeleton_Week7 " = "#a5a9c6"),
                          "cluster.origin" = c("hNC_Axial_skeleton_Week7 " = "green",
                                               "hCT_Axial_skeleton_Week7 " = "#a26d85"))

ha <- HeatmapAnnotation(df = ordered_meta_data,
                        show_annotation_name = TRUE, col = annotation_colors)

genes_to_use <- c("ALCAM",
                  "SHH",
                  "CHRM2",
                  "SPP1",
                  "RAB3B",
                  "TBXT",
                  "CA3",
                  "HOPX",
                  "JUNB",
                  "IER2",
                  "FOS",
                  "UCHL1",
                  "VAMP5",
                  "TSPAN7",
                  "IFI6",
                  "TAGLN2",
                  "TUBB6",
                  "FRZB",
                  "TNFRSF11B",
                  "CNTNAP2",
                  "NKD1",
                  "PEG10",
                  "NTM",
                  "EDNRB",
                  "NOG",
                  "SFRP5",
                  "SCRG1",
                  "MPZ",
                  "CADM2",
                  "MIA",
                  "COL25A1",
                  "COL4A1",
                  "DIAPH3",
                  "SERPINE2",
                  "MCAM",
                  "FSTL5",
                  "EGFL7",
                  "CYBA",
                  "HSP90B1",
                  "PLEKHB1",
                  "PERP",
                  "SEMA3C",
                  "S100A6",
                  "SOX9",
                  "SORBS2",
                  "COL9A3",
                  "ACAN",
                  "COL11A2",
                  "COL2A1",
                  "COL9A1",
                  "COL11A1",
                  "COL9A2",
                  "LRP2",
                  "BGN",
                  "MARCKS",
                  "FLRT2",
                  "CYP1B1",
                  "COLEC12",
                  "FBN2",
                  "SPOCK3",
                  "FOXC2",
                  "FOXC1",
                  "FOXD1",
                  "SIX1",
                  "NR2F1",
                  "TWIST1",
                  "VIM",
                  "TCF15")


S3 <- ScaleData(S3, genes.use = genes_to_use)

my_data <- S3@assays[["RNA"]]@scale.data

my_data <- my_data[, rownames(ordered_meta_data)]

my_data <- subset(my_data, rownames(my_data) %in% genes_to_use)

gene_list = rep(c("brown_gene", "cyan_gene"), times = c(53, 15))

row_gene <- data.frame(gene_list, row.names = genes_to_use)

row_gene <- subset(row_gene, rownames(row_gene) %in% rownames(my_data))

color_annotation <- list("gene_list" = c("brown_gene" = "brown", "cyan_gene" ="cyan"))

row_ha <- rowAnnotation(df = row_gene,
                        show_annotation_name = TRUE,
                        col = color_annotation)

my_data <- my_data[order(match(row.names(my_data), row.names(row_gene))), ]

row_gene <- arrange(row_gene, order(match(row.names(row_gene), row.names(my_data))))

my_data <- my_data[rownames(row_gene), ]

row_gene$gene_name <- c("ALCAM",
                        "SHH",
                        "CHRM2",
                        "SPP1",
                        "RAB3B",
                        "TBXT",
                        "CA3",
                        "HOPX",
                        "JUNB",
                        "IER2",
                        "FOS",
                        "UCHL1",
                        "VAMP5",
                        "TSPAN7",
                        "IFI6",
                        "TAGLN2",
                        "TUBB6",
                        "FRZB",
                        "TNFRSF11B",
                        "CNTNAP2",
                        "NKD1",
                        "PEG10",
                        "NTM",
                        "EDNRB",
                        "NOG",
                        "SFRP5",
                        "SCRG1",
                        "MPZ",
                        "CADM2",
                        "MIA",
                        "COL25A1",
                        "COL4A1",
                        "DIAPH3",
                        "SERPINE2",
                        "MCAM",
                        "FSTL5",
                        "EGFL7",
                        "CYBA",
                        "HSP90B1",
                        "PLEKHB1",
                        "PERP",
                        "SEMA3C",
                        "S100A6",
                        "SOX9",
                        "SORBS2",
                        "COL9A3",
                        "ACAN",
                        "COL11A2",
                        "COL2A1",
                        "COL9A1",
                        "COL11A1",
                        "COL9A2",
                        "LRP2",
                        "BGN",
                        "MARCKS",
                        "FLRT2",
                        "CYP1B1",
                        "COLEC12",
                        "FBN2",
                        "SPOCK3",
                        "FOXC2",
                        "FOXC1",
                        "FOXD1",
                        "SIX1",
                        "NR2F1",
                        "TWIST1",
                        "VIM",
                        "TCF15")
## ---------------------------

## Create the heatmap object

set.seed(123)

h_cluster <- Heatmap(my_data,
                     heatmap_legend_param = list(
                       title = "gene_expression", at = c(-4, -2, 0, 2, 4, 6, 8, 10)),
                     show_column_names = FALSE,
                     top_annotation = ha,
                     right_annotation = row_ha,
                     column_split = factor(ordered_meta_data$cluster.origin,
                                           levels = c("hNC_Axial_skeleton_Week7 ",
                                                      "hCT_Axial_skeleton_Week7 ")),
                     cluster_columns = TRUE,
                     cluster_column_slices = FALSE,
                     cluster_row_slices = FALSE,
                     row_split = factor(row_gene$gene_name,
                                        levels = c("ALCAM",
                                                   "SHH",
                                                   "CHRM2",
                                                   "SPP1",
                                                   "RAB3B",
                                                   "TBXT",
                                                   "CA3",
                                                   "HOPX",
                                                   "JUNB",
                                                   "IER2",
                                                   "FOS",
                                                   "UCHL1",
                                                   "VAMP5",
                                                   "TSPAN7",
                                                   "IFI6",
                                                   "TAGLN2",
                                                   "TUBB6",
                                                   "FRZB",
                                                   "TNFRSF11B",
                                                   "CNTNAP2",
                                                   "NKD1",
                                                   "PEG10",
                                                   "NTM",
                                                   "EDNRB",
                                                   "NOG",
                                                   "SFRP5",
                                                   "SCRG1",
                                                   "MPZ",
                                                   "CADM2",
                                                   "MIA",
                                                   "COL25A1",
                                                   "COL4A1",
                                                   "DIAPH3",
                                                   "SERPINE2",
                                                   "MCAM",
                                                   "FSTL5",
                                                   "EGFL7",
                                                   "CYBA",
                                                   "HSP90B1",
                                                   "PLEKHB1",
                                                   "PERP",
                                                   "SEMA3C",
                                                   "S100A6",
                                                   "SOX9",
                                                   "SORBS2",
                                                   "COL9A3",
                                                   "ACAN",
                                                   "COL11A2",
                                                   "COL2A1",
                                                   "COL9A1",
                                                   "COL11A1",
                                                   "COL9A2",
                                                   "LRP2",
                                                   "BGN",
                                                   "MARCKS",
                                                   "FLRT2",
                                                   "CYP1B1",
                                                   "COLEC12",
                                                   "FBN2",
                                                   "SPOCK3",
                                                   "FOXC2",
                                                   "FOXC1",
                                                   "FOXD1",
                                                   "SIX1",
                                                   "NR2F1",
                                                   "TWIST1",
                                                   "VIM",
                                                   "TCF15")),
                     raster_device = c("png"),
                     cluster_rows = TRUE,
                     row_names_gp = gpar(fontsize = 8),
                     row_title_rot = 0,
                     row_title_gp = gpar(fontsize = 6),
                     column_title_gp = gpar(fontsize = 6),
                     column_title_rot = 45,
                     column_gap = unit(5, "mm"), 
                     row_gap = unit(0, "mm"))

## ---------------------------

## Draw and save the heatmap

png(file = "results/heatmap_7week3.png",
    width = 10000,
    height = 10000, res = 400)
print(h_cluster)
dev.off()

pdf(file = "results/heatmap_7week3.pdf", width = 20,
    height = 30)
print(h_cluster)
dev.off()


