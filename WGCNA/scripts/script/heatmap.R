#$ -S /CONDAS/users/schevolleau/embryo-heatmap/bin/Rscript
#$ -e	/LAB-DATA/BiRD/users/schevolleau/rlog/
#$ -o /LAB-DATA/BiRD/users/schevolleau/rlog/

# PACKAGES
library(rjson)
library(ComplexHeatmap)

# FUNCTIONS
source("https://gitlab.univ-nantes.fr/E114424Z/veneR/raw/master/loadFun.R?inline=false")
source("https://gitlab.univ-nantes.fr/E198672Y/embryo-functions/-/raw/master/functions.R?inline=false")

# VARIABLES
## IF FILE IS LOAD WITH SNAKEMAKE
if(exists("snakemake")){
    settings <- fromJSON(file = snakemake@input[["heatmap_config_path"]])
    markers <- readTable(snakemake@input[["markers_path"]], rnames=FALSE)
    exprMat <- readTable(snakemake@input[["expr_Mat_path"]])
    try(annotation <- readTable(snakemake@input[["annotation_path"]]), silent = T)
    output_folder <- snakemake@params[["output_folder"]]
}else{
    settings <- fromJSON(file = "/LAB-DATA/BiRD/users/schevolleau/git/embryo-heatmap/settings_blastoid.json")
    markers <- readTable(settings$markers_path, rnames=FALSE)
    exprMat <- readTable(settings$expr_Mat_path)
    try(annotation <- readTable(settings$annotation_path), silent = T)
    output_folder <- settings$output_folder
}

## DATA
try(TF <- readTable(settings$TF_path), silent = T)
try(annot <- annotation[, settings$annotation], silent = T)
if (settings$group_cells_order == "None"){
    try(group <- annotation[, settings$group_cells], silent = T)
}else{
    try(group <- factor(annotation[, settings$group_cells], levels=settings$group_cells_order), silent = T)
}
try(markers <- markers[markers[, settings$marker_selection], ], silent = T)

## PREPROCESSING 
if (settings$log == "True")
    exprMat <- log2(exprMat + 1)
if (!(settings$TF_path == "None"))
  markers <- markers[markers[, settings$genes_col] %in% rownames(TF),]

if (settings$marker_split == "None"){
    markers <- markers[settings$genes_col]
    li_markers <- list()
    li_markers[[names(markers)]] <- markers[, settings$genes_col]
    markers <- li_markers
}else{
    markers <- split(markers[, settings$genes_col], markers[, settings$marker_split])
}

dir.create(output_folder, recursive = T)

# PLOTS
## HEATMAP OPTIONS
options_list <- list(showGrid = F,
                     use_raster = T,
                     raster_quality = 5,
                     returnHeatmap = TRUE,
                     scale=ifelse (settings$scale == "True", T, F),
                     center= ifelse (settings$center == "True", T, F),
                     cluster_column_slices=F,
                     cluster_columns = ifelse (settings$cluster_col == "True", T, F),
                     column_title_rot=90,
                     row_title_gp=gpar(fontsize=4),
                     column_title_gp=gpar(fontsize=5),
                     border = TRUE, name=ifelse (settings$log == "True", "Centred log2(x+1) Expression", "Centred Expression"),
                     row_title_rot = 0, show_row_names = ifelse(settings$printRownames=="True", TRUE, FALSE),
                     autoFontSizeRow = T)

if (!(settings$samplesMean == "None"))
    options_list[["bottom_annotation"]] <- columnAnnotation(sample=anno_lines(colSums(exprMat[genes,]), axis_param = list(direction = "reverse")), annotation_name_gp= gpar(fontsize = 0))
if (!(settings$genesMean == "None"))
    options_list[["right_annotation"]] <- rowAnnotation(genes=anno_lines(rowSums(exprMat[genes,])), annotation_name_gp= gpar(fontsize = 0))
if (settings$row_order == "True"){
    options_list[["row_order"]] <- genes
    options_list[["cluster_rows"]] <- F
}else{
    options_list[["cluster_rows"]] <- T
}
if (!(settings$annotation_path == "None"))
    if (!(settings$group_cells == "None"))
        options_list[["column_split"]] <- group
if (!(settings$color_annotation == "None")){
    color_annotation <- settings$color_annotation
    names(color_annotation) <- settings$color_order
}
if (!(settings$annotation == "None"))
    options_list[["top_annotation"]] <- HeatmapAnnotation(annot = factor(annot, levels=c(settings[["color_order"]])), col = list(annot = color_annotation), annotation_name_side = "right")
if (!(settings$colnames == "None")){
    options_list[["show_column_names"]] <- T
}else{
    options_list[["show_column_names"]] <- F
}
options_list[["autoFontSizeRow"]] <- T
options_list_separately <- options_list

## HEATMAP CREATION
htList <- list()
for(marker_name in names(markers)){
    uppercase_genes <- unlist(lapply(markers[[marker_name]], toupper))
    genes <- intersect(uppercase_genes, rownames(exprMat))
    options_list["row_title"] <- paste0(marker_name,"\n",length(genes)," genes")
    options_list_separately["row_title"] <- paste0(marker_name,"\n",length(genes)," genes")

    if (length(genes) == 0){
        next
    }else{
        # CREATE HEATMAP
        ## ALL MODULES IN 1 HEATMAP
        htList[[marker_name]] <- do.call("heatmap.DM3", c(list(matrix = exprMat[genes,]), options_list))
        options_list["top_annotation"] <- NULL

        ## 1 MODULE FOR 1 HEATMAP
        pdf(paste0(output_folder, "/", marker_name, "_heatmap.pdf"), useDingbats = F)
            draw(do.call("heatmap.DM3", c(list(matrix = exprMat[genes,]), options_list_separately)))
        dev.off()
    }
}

### CONCATENATE HEATMAPS ON VERTICAL AXE
ht <- htList[[names(markers)[1]]]
if (length(names(markers)) > 1)
    for (cluster in names(markers)[2:length(markers)])
        ht <- ht %v% htList[[cluster]]

pdf(paste0(output_folder, "/heatmap.pdf"), useDingbats = F)
    draw(ht, cluster_columns = T, cluster_column_slices = F, clustering_method_columns="ward.D2", gap = unit(0.5, "mm"), clustering_distance_columns="pearson")
dev.off()

pdf("file_to_remove")
### SAVE GENES ORDER
genes_order <- list()
for(marker in names(htList))
    genes_order[[marker]] <- intersect(markers[[marker]],rownames(exprMat))[row_order(htList[[marker]])]

write(toJSON(genes_order), paste0(output_folder, "/genes_order.json"))
if (length(names(markers)) > 1){
    df_genes_order <- listToDf(genes_order)
    if (!(settings$TF_path == "None")){
        df_genes_order[df_genes_order$gene %in% TF[,1], "TF"] <- TRUE
        df_genes_order[!(df_genes_order$gene %in% TF[,1]), "TF"] <- FALSE
    }
    writeTable(df_genes_order, paste0(output_folder, "/genes_order.tsv"))
}

### SAVE SAMPLES ORDER
samples_order <- list()
samples_order_index <- column_order(htList[[names(markers)[1]]])
if (settings$group_cells == "None"){
   samples_order[["Gene"]] <- colnames(exprMat)[samples_order_index ]
}else{
    for (cluster_annot in names(samples_order_index))
        samples_order[[cluster_annot]] <- colnames(exprMat)[samples_order_index[[cluster_annot]]]
}

write(toJSON(samples_order), paste0(output_folder, "/samples_order.json"))
if (settings$group_cells == "None")
    writeTable(listToDf(samples_order), paste0(output_folder, "/samples_order.tsv"))

dev.off()
save.image(file=paste0(output_folder, "/env.Rdata"))
