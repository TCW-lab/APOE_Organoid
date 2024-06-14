library(dplyr)
library(patchwork)
library(ggplot2)
library(gridExtra)
library(Matrix)
library(scales)
library(cowplot)
library(Seurat) #Seurat v.3.0
library(RColorBrewer)
library(BiocManager) #v 3.17
library(tibble)
library(patchwork)
# for miQC / flexmix:
library(flexmix)
library(miQC)
library(singleCellTK)
library(fgsea)
library(DropletUtils)
library(SeuratData)
library(tidyr)
library(SeuratWrappers)
library(data.table)
library(SoupX)
library(fishpond)
library(Matrix)
setwd("/projectnb/tcwlab/LabMember/akg/projects/APOE_Organoid")

dir <- "/projectnb/tcwlab/LabMember/akg/projects/APOE_Organoid"
out <- '/outputs/04-Clustering_and_Batch_Correction'
organoid <- readRDS(paste0(dir, '/outputs/03-miQC_and_SoupX/organoid4.rds'))

##### Step 1. Read in RDS object and consolidate counts slots ##################

names(organoid@assays$RNA@layers)  # This will show individual count layers for each sample

#### Here we use the JoinLayers function in order to rejoin all of the counts
#### and data layers into just two layers (currently it is in 24 for each).
organoid <- JoinLayers(organoid)


# saveRDS(organoid, paste0(dir, out, '/organoid.rds')) 




# Step 2. Perform UMAP Clustering and Visualization ############################
### This is if you have not already performed SCTransform, or prefer not to use
### SCTransform due to wanting to define the parameters you use to normalize,
### FindVariableFeats, Scale, etc.
# Preprocessing steps
# organoid <- NormalizeData(organoid, verbose = FALSE)
organoid <- FindVariableFeatures(organoid, selection.method = "vst", nfeatures = 2000)
organoid <- ScaleData(organoid, features = rownames(organoid), verbose = TRUE)

# Dimensionality reduction
organoid <- RunPCA(organoid, features = VariableFeatures(object = organoid), verbose = TRUE)
organoid <- RunUMAP(organoid, dims = 1:50)

# Clustering
organoid <- FindNeighbors(organoid, dims = 1:50, k.param = 30, verbose = TRUE)
organoid <- FindClusters(organoid, resolution = 0.6, verbose = TRUE)

DimPlot(organoid, reduction = 'umap', group.by = 'seurat_clusters')


saveRDS(organoid, file = paste0(dir, '/outputs/04-Clustering_and_Batch_Correction/organoid1.rds')) # 06-12 left off here.
organoid <- readRDS(paste0(dir,out,'/organoid.rds'))

organoid@assays$RNA@layers$data

### Inspect the organoid object here.
DimPlot(organoid, reduction = 'umap', group.by = 'seurat_clusters')

organoid <- readRDS(paste0(dir, out, '/organoid1.rds'))

#### Step 3. Perform Harmony Reduction for Batch Effect Correction  ############
library(harmony)

organoid <- RunHarmony(organoid, group.by.vars = "sample")

# After Harmony integration, proceed with re-clustering
organoid <- FindNeighbors(organoid, dims = 1:50, k.param = 30, reduction = "harmony")
organoid <- FindClusters(organoid, resolution = 0.6, verbose = TRUE)

organoid <- RunUMAP(organoid, reduction = "harmony", dims = 1:50, reduction.name='humap')
DimPlot(organoid, reduction = "humap", group.by = "seurat_clusters")



saveRDS(organoid, file = paste0(dir, '/outputs/05-Doublet_Finder/organoid2.rds'))
