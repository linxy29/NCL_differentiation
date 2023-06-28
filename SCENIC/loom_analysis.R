# Read file https://rdrr.io/github/aertslab/SCENIC/f/vignettes/importing_pySCENIC.Rmd
# PACKAGES
##install.packages("devtools")
##devtools::install_github("aertslab/SCopeLoomR")
#devtools::install_github("aertslab/SCENIC", ref="v1.1.2")
library(SCopeLoomR)
library(SummarizedExperiment)
library(SCENIC)
library(AUCell)
library(rjson)
# FUNCTIONS
source("https://gitlab.univ-nantes.fr/E198672Y/embryo-functions/-/raw/master/functions.R?inline=false")

# VARIABLES
settings <- fromJSON(file = "/SCRATCH-BIRD/users/nantunes/embryo-scenic/loom_analysis_settings.json")

#settings <- list(
#  "WD"="/SCRATCH-BIRD/users/nantunes/embryo-scenic/output",
#  "loom"="scRNAseq_SCENIC.loom",
#  "output_dir"="result"
#)

# WD
setwd(settings$WD)

dir.create(settings$output_dir)

pyScenicLoomFile <- file.path(getwd(), settings$loom)
loom <- open_loom(pyScenicLoomFile, mode="r+")

cellInfo <- get_cell_annotation(loom)
regulons_incidMat <- get_regulons(loom, column.attr.name="Regulons") # as incid matrix
regulons_motif <- regulonsToGeneLists(regulons_incidMat) # convert to list

# Regulon AUC and thresholds
regulonsAUC <- as.matrix(assay(get_regulons_AUC(loom, column.attr.name="RegulonsAUC")))
regulonsAucThresholds <- get_regulon_thresholds(loom)

# cellInfo <- get_cellAnnotation(loom) # will also contain AUC values, etc... you can filter them out
clusterings <- get_clusterings_with_name(loom)

# SAVING
## RDS
saveRDS(regulons_motif, paste0(settings$output_dir, "/regulons_motif.rds"))
saveRDS(regulonsAUC, paste0(settings$output_dir, "/regulonsAUC.rds"))
saveRDS(regulonsAucThresholds, paste0(settings$output_dir, "/regulonsAucThresholds.rds"))
saveRDS(clusterings, paste0(settings$output_dir, "/clusterings.rds"))

## JSON AND TSV
regulons_motif_df <- listToDf(regulons_motif, "regulons", "genes")
writeTable(regulons_motif_df, paste0(settings$output_dir, "/regulons_motif.tsv"))
exportJSON <- toJSON(regulons_motif)
write(exportJSON, paste0(settings$output_dir, "/regulons_motif.json"))
df <- data.frame(Regulons=regulonsAucThresholds)
rownames(df) <- df$Regulons
df$Thresholds <- names(regulonsAucThresholds)
writeTable(df, paste0(settings$output_dir, "/regulonsAucThresholds.tsv"))
binary_activity <- list()
for (reg in rownames(regulonsAUC))
 binary_activity[[reg]] <- regulonsAUC[reg, ] > df[reg, "Thresholds"]
binary_activity_regulons <- as.data.frame(binary_activity)
writeTable(binary_activity_regulons, paste0(settings$output_dir, "/binary_activityaryRegulonsCell.tsv"))
write(exportJSON, paste0(settings$output_dir, "/embeddings.json"))
writeTable(clusterings, paste0(settings$output_dir, "/clusterings.tsv"))
temp <- strsplit(names(regulons_motif), "_")
names(regulons_motif) <- sapply(temp,"[[",1)
exportJSON <- toJSON(regulons_motif)
write(exportJSON, paste0(settings$output_dir, "/regulons_motif_splitted_names.json"))

