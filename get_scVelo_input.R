library(Seurat)
library(SeuratDisk)
library(SeuratWrappers)
library(tidyverse)


pathSPC = "/storage/holab/linxy/iPSC/veloAE/"  ## for server
sample = "cfs"
#pathSPC = "/mnt/d/Data"   ## for WSL
#pathSPC = "D:/Data"    ## for windows
options(future.globals.maxSize = 40000 * 1024^2)

## eusample
seuratObj = readRDS("/storage/holab/linxy/iPSC/seuratObj/Dataset_EU_D7.RDS")
## change the annotation
#cell.type=c("Neuromesodermal progenitors","Lateral/paraxial mesoderm progenitors",
            "Axial progenitors", "Primitive streak (cluster3)","Primitive streak (cluster4)",
            "Cardiac mesoderm","Notochord", "Endoderm", "Neuroectoderm")
#Idents(seuratObj)=seuratObj$seurat_clusters0.4
#names(x = cell.type) <- levels(seuratObj)
#seuratObj <- RenameIdents(seuratObj,cell.type)
#saveRDS(seuratObj, "/storage/holab/linxy/iPSC/seuratObj/Dataset_EU_D7.RDS")

#integrateObj = readRDS("/storage/holab/linxy/iPSC/seuratObj/PBMC_TC_clustered_annotated.RDS")
DimPlot(seuratObj)
#Idents(seuratObj) = seuratObj$annotation
ggsave(str_c("/storage/holab/linxy/iPSC/seuratObj/", sample, "_UMAP.png"))

## ncc sample
seuratObj = readRDS("/usersdata/share/linxy/2023-03-17_NCC_clustered_annotated.RDS")
#cell.type=c("Sclerotome","Somitic mesoderm","lateral/paraxial mesoderm", "axial skeleton system","notochord")
#Idents(seuratObj)=seuratObj$SCT.CC.reg_snn_res.0.2
#names(x = cell.type) <- levels(seuratObj)
#seuratObj <- RenameIdents(seuratObj,cell.type)
#saveRDS(seuratObj, "/storage/holab/linxy/iPSC/seuratObj/NCC_230330.RDS")
#Idents(seuratObj) = seuratObj$annotation
DimPlot(seuratObj)
ggsave(str_c("/storage/holab/linxy/iPSC/seuratObj/", sample, "_UMAP.png"))

## cfs sample
seuratObj = readRDS("/usersdata/share/linxy/2023-03-17_CFSEU_annotated.RDS")


## cfseu integration
seuratObj = readRDS("/storage/holab/linxy/iPSC/seuratObj/2023-03-27-Integration.by.Seurat.CFS.EU.SCT.CC.clustered.rds")
cell.type=c("neuromesodermal progenitors","primitive streak", "cardiac", "paraxial mesoderm progenitors","notochord",
"axial progenitors","neuroectoderm", "endoderm")
Idents(seuratObj)=seuratObj$integrated_snn_res.0.3
names(x = cell.type) <- levels(seuratObj)
seuratObj <- RenameIdents(seuratObj,cell.type)
saveRDS(seuratObj, "/storage/holab/linxy/iPSC/seuratObj/cfseu_230415.RDS")

## Plot and save
DimPlot(seuratObj)
head(Cells(seuratObj))
write.csv(Cells(seuratObj), file = str_c(pathSPC, sample, "_cellID_obs.csv"), row.names = FALSE)
head(Embeddings(seuratObj, reduction = "umap")[,1:2])
write.csv(Embeddings(seuratObj, reduction = "umap")[,1:2], file = str_c(pathSPC, sample, "_cell_embeddings.csv"))
head(Idents(seuratObj))
write.csv(Idents(seuratObj), file = str_c(pathSPC, sample,"_clusters.csv"))





