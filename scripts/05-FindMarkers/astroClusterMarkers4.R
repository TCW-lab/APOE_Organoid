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

# organoid_APOE33Ch <- subset(organoid, subset = genotype == 'APOE33Ch')

library(scRNAseq)
library(SingleR)
library(tidyverse)

astrocytes <- c('3', '4', '12', '13', #this is to be a reminder
                '14', '15', '17') # of all of the astrocyte clusters

ident = 4
# rename all the other Astrocyte clusters as Astrocyte
organoid <- RenameIdents(organoid, '3' = "Astrocyte")
organoid <- RenameIdents(organoid, '12' = "Astrocyte")
organoid <- RenameIdents(organoid, '13' = "Astrocyte")
organoid <- RenameIdents(organoid, '14' = "Astrocyte")
organoid <- RenameIdents(organoid, '15' = "Astrocyte")
organoid <- RenameIdents(organoid, '17' = "Astrocyte")

#compare these Astrocytes to the rest of the population
markers <- FindMarkers(organoid, ident.1 = ident, ident.2 = "Astrocyte")

filename <- paste0("Astrocyte_Markers_Cluster_",ident,".csv")
write.csv(markers, file = filename, row.names = TRUE)