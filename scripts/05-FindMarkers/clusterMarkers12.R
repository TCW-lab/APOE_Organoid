library(tools)
library(Seurat)
library(data.table)
library(ggplot2)
library(SingleCellExperiment)
library(fishpond)
library(scater)
library(org.Hs.eg.db)
library(Matrix)
library(biomaRt)
library(harmony)
library(ggthemes)
setwd("/projectnb/tcwlab/LabMember/akg/APOE_Jorganoid_project/outputs/03-CellTypeAnnotation")

dir <- "/projectnb/tcwlab/LabMember/akg/APOE_Jorganoid_project/outputs/03-CellTypeAnnotation"

organoid <- readRDS('/projectnb/tcwlab/LabMember/akg/APOE_Jorganoid_project/outputs/03-CellTypeAnnotation/organoid_post_harmony.rds')

library(scRNAseq)
library(SingleR)
library(tidyverse)

ident = 12
filename <- paste0("/projectnb/tcwlab/LabMember/akg/APOE_Jorganoid_project/outputs/05-FindMarkers/Markers_Cluster_",ident,".csv")


markers <- FindMarkers(organoid, ident.1 = ident)
write.csv(markers, file = filename, row.names = TRUE)