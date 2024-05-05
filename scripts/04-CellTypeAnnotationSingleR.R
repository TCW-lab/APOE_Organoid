#### Script for Cell Type Annotation
#### 02-14-2024 Andrew Gjelsteen 

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

organoid <- readRDS('/projectnb/tcwlab/LabMember/akg/APOE_Jorganoid_project/outputs/03-CellTypeAnnotation/organoid_post_harmony.rds')

library(scRNAseq)
library(SingleR)
library(tidyverse)
?FindMarkers
DimPlot(organoid, reduction = "humap", group.by = "genotype")

table(organoid$seurat_clusters)
markers <- FindMarkers(organoid, ident.1 = 20)

# remotes::install_version("SingleR", version = "1.0.1")
# remotes::install_github('dviraran/SingleR')

#Reference Atlas: Polioudakis et. al.

ref_dir <- "/projectnb/tcwlab/LabMember/akg/APOE_Jorganoid_project/ref_data"
organoid <- readRDS('/projectnb/tcwlab/LabMember/akg/APOE_Jorganoid_project/outputs/03-CellTypeAnnotation/organoid_post_harmony.rds')
#load count matrix
load(paste(ref_dir, "/sc_dev_cortex_geschwind/raw_counts_mat.rdata", sep = ""))

raw.counts.mat <- as.matrix(raw_counts_mat)
reffy <- CreateSeuratObject(raw.counts.mat)

organoid@meta.data$nCount_RNA



#load metadata
ref_metadata <- read.csv(file = paste(ref_dir, "/sc_dev_cortex_geschwind/cell_metadata.csv", sep = ""))
rownames(ref_metadata) <- ref_metadata[,1]
ref_metadata[,1] <- NULL
ref_metadata[,1]
#append metadata to reference object
reffy <- AddMetaData(reffy, ref_metadata)
reffy <- reffy[,!is.na(reffy$Cluster)]

#convert to a single cell experiment
reffo <- as.SingleCellExperiment(reffy)
metadata(reffo)<- ref_metadata

#normalize data
library(scuttle)
reffo <- logNormCounts(reffo)

str(reffo@metadata$Cluster)
### Need to work with this as.matrix() function call.
### Seems like this is what is causing the issue at hand.
# reffo@metadata$Cluster <- as.matrix(reffo@metadata$Cluster)
colData(reffo)$Cluster <- reffo@metadata$Cluster

samples <- unique(organoid@meta.data$sample)

for (sample_id in samples[18:24]) {
  sce <- as.SingleCellExperiment(subset(organoid, subset = sample == sample_id))
  keep_rows <- !grepl("^ENSG", rownames(sce))
  # keep only gene symbols.
  sce <- sce[keep_rows, ]
  sce_counts <- counts(sce)
  reffo_counts <- counts(reffo)
  
  # intersect names of sce / reffo
  common_row_names <- intersect(rownames(sce), rownames(reffo))
  
  # subset sce to only be rownames in common w reffo
  sce_filtered <- sce[common_row_names, ]
  
  # Check to ensure 'sce_filtered' contains only the common rows
  rownames(sce_filtered) 
  
  sce_counts_dense <- as.matrix(sce_counts)
  reffo_counts_dense <- as.matrix(reffo_counts)
  # Run SingleR
  singleR_results <- SingleR(sc_data = sce_counts_dense,
                             ref_data = reffo_counts_dense, types = factor(reffo@colData$Cluster),
                             numCores = 14)
  filename = paste0('organoid_SingleR_',sample_id,'.rds')
  saveRDS(singleR_results, file.path(dir, filename))
}

# ?SingleR
#### use maxcutoff parameter in FeaturePlots (to adjust scaling)
### ALD1H1



