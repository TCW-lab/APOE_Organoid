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

setwd("/projectnb/tcwlab/LabMember/akg/projects/APOE_Organoid")

dir <- "/projectnb/tcwlab/LabMember/akg/projects/APOE_Organoid"
out <- '/outputs/04-Clustering_and_Batch_Correction'
organoid <- readRDS(paste0(dir, '/outputs/03-miQC_and_SoupX/organoid4.rds'))

##### Step 1. Read in RDS object and consolidate counts slots ##################
##### This is being done in order to reduce Seurat object size (currently 20GB 
##### for me).

organoid[["RNA"]] <- as(object = organoid[["SoupX_counts"]], Class = "Assay5")
# Remove the SoupX_counts assay
organoid[["SoupX_counts"]] <- NULL

# Verify that only the RNA assay remains
sum(organoid@assays$RNA@layers$counts)
# This is equal to my output for SoupX_Counts from the previous script.
# So, this is good.

# saveRDS(organoid, paste0(dir, out, '/organoid.rds'))



# Step 2. Perform UMAP Clustering and Visualization ############################

# Preprocessing steps
organoid <- NormalizeData(organoid, verbose = FALSE)
organoid <- FindVariableFeatures(organoid, selection.method = "vst", nfeatures = 2000)
organoid <- ScaleData(organoid, features = rownames(organoid), verbose = TRUE)

# Dimensionality reduction
organoid <- RunPCA(organoid, features = VariableFeatures(object = organoid), verbose = TRUE)
organoid <- RunUMAP(organoid, dims = 1:50)

# Clustering
organoid <- FindNeighbors(organoid, dims = 1:50, k.param = 30, verbose = TRUE)
organoid <- FindClusters(organoid, resolution = 0.6, verbose = TRUE)

DimPlot(organoid, reduction = 'umap', group.by = 'seurat_clusters')

saveRDS(organoid, file = paste0(dir, '/', 'outputs/04-Clustering_and_Batch_Correction/organoid1.rds'))

organoid <- readRDS(paste0(dir, '/outputs/04-Clustering_and_Batch_Correction/organoid1.rds'))



### Inspect the organoid object here.
DimPlot(organoid, reduction = 'umap', group.by = 'seurat_clusters')

organoid[["RNA"]] <- as(object = organoid[["SoupX_counts"]], Class = "Assay5")
# Remove the SoupX_counts assay
organoid[["SoupX_counts"]] <- NULL


#### Step 3. Perform Harmony Reduction for Batch Effect Correction  ############
library(harmony)

organoid <- RunHarmony(organoid, group.by.vars = "sample")

# After Harmony integration, proceed with re-clustering
organoid <- FindNeighbors(organoid, dims = 1:50, k.param = 30, reduction = "harmony")
organoid <- FindClusters(organoid, resolution = 0.6, verbose = TRUE)

organoid <- RunUMAP(organoid, reduction = "harmony", dims = 1:50, reduction.name='humap')
DimPlot(organoid, reduction = "humap", group.by = "seurat_clusters")


organoid <- readRDS(paste0(dir, '/outputs/04-Clustering_and_Batch_Correction/organoid1.rds'))

saveRDS(organoid, file = paste0(dir, '/outputs/05-Cluster_Annotation/organoid.rds'))
        
#### Step 4. Run DoubletFinder to remove possible doublets #####################

library(DoubletFinder)
## The following code is from 
## Chris McGinnis (UCSF): https://github.com/chris-mcginnis-ucsf/DoubletFinder

# SCTransform organoid
library(sctransform)
SCTransform(organoid, useNames = TRUE)

## pK Identification (no ground-truth) ---------------------------------------------------------------------------------------
sweep.res.list_organoid <- paramSweep(organoid, PCs = 1:10, sct = FALSE)
sweep.stats_organoid <- summarizeSweep(sweep.res.list_organoid, GT = FALSE)
bcmvn_organoid <- find.pK(sweep.stats_organoid)

## Homotypic Doublet Proportion Estimate -------------------------------------------------------------------------------------
homotypic.prop <- modelHomotypic(annotations)           ## ex: annotations <- organoid@meta.data$ClusteringResults
nExp_poi <- round(0.075*nrow(organoid@meta.data))  ## Assuming 7.5% doublet formation rate - tailor for your dataset
nExp_poi.adj <- round(nExp_poi*(1-homotypic.prop))

## Run DoubletFinder with varying classification stringencies ----------------------------------------------------------------

organoid <- doubletFinder(organoid, PCs = 1:10, pN = 0.25, pK = 0.09, 
                            nExp = nExp_poi.adj, sct = FALSE)

### Now save the doublet-removed organoid object.
saveRDS(organoid, file = paste0(dir, '/outputs/05-Cluster_Annotation/organoid1.rds'))
