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


dir <- "/projectnb/tcwlab/LabMember/akg/APOE_Jorganoid_project/outputs/03-CellTypeAnnotation"


library(scRNAseq)
library(SingleR)
library(tidyverse)

S1 <- readRDS('/projectnb/tcwlab/LabMember/akg/APOE_Jorganoid_project/outputs/03-CellTypeAnnotation/organoid_SingleR_sample1.rds')

organoid <- readRDS('/projectnb/tcwlab/LabMember/akg/APOE_Jorganoid_project/outputs/03-CellTypeAnnotation/organoid_post_harmony.rds')

DimPlot(organoid, reduction = 'humap', group.by = "sample")

DimPlot(organoid, reduction = "umap", group.by = "sample")




organoid <- AddMetaData(organoid, S1$labels)

organoid@meta.data$V1
organoid
DimPlot(subset(organoid, subset = sample == "sample1"), group.by = "V1", reduction = "humap")
FeaturePlot(organoid, features =c("GAD1", "GAD2"), reduction = 'humap')
        