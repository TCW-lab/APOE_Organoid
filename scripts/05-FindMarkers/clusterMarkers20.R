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
setwd("/projectnb/tcwlab/LabMember/akg/APOE_Jorganoid_project/outputs/05-FindMarkers")

dir <- "/projectnb/tcwlab/LabMember/akg/APOE_Jorganoid_project/outputs/05-FindMarkers"

organoid <- readRDS('/projectnb/tcwlab/LabMember/akg/APOE_Jorganoid_project/outputs/03-CellTypeAnnotation/organoid_post_harmony.rds')

organoid_APOE33Ch <- subset(organoid, subset = genotype == 'APOE33Ch')

library(scRNAseq)
library(SingleR)
library(tidyverse)

ident = 20
filename <- paste0("Markers_Cluster_",ident,".csv")


markers <- FindMarkers(organoid_APOE33Ch, ident.1 = ident)
write.csv(markers, file = filename, row.names = TRUE)