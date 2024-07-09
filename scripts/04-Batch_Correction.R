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

# This will show individual count layers for each sample
organoid@assays$SoupX_counts

organoid@assays$SoupX_counts
#### Here we use the JoinLayers function in order to rejoin all of the counts
#### and data layers into just two layers (currently it is in 24 for each).
organoid <- JoinLayers(organoid, assay = "RNA")


n_APOE22 <- sum(organoid@meta.data$genotype == 'APOE22')
n_APOE33 <- sum(organoid@meta.data$genotype == 'APOE33')
n_APOE33Ch <- sum(organoid@meta.data$genotype == 'APOE33Ch')
n_APOE44 <- sum(organoid@meta.data$genotype == 'APOE44')
# Display the result of the miQC filtering:
p <- VlnPlot(organoid, pt.size = 0, group.by = 'sample', features = 'nFeature_RNA', log = TRUE) +
  theme(legend.position = "none")
print(p)


# Define the color and point size
point_color <- "lightcoral"  # a light shade of red
point_size <- 0.1  # small point size

# Create the first FeatureScatter plot with custom color and small point size
plot1 <- FeatureScatter(organoid, feature1 = "nCount_RNA", feature2 = "percent.mt") +
  geom_point(color = point_color, size = point_size) +
  theme(legend.position = 'none')

# Create the second FeatureScatter plot with custom color and small point size
plot2 <- FeatureScatter(organoid, feature1 = "nCount_RNA", feature2 = "nFeature_RNA") +
  geom_point(color = point_color, size = point_size) +
  theme(legend.position = 'none')

# Combine and display the plots
p <- plot1 + plot2
print(p)


rm(plot1, plot2, p)


# saveRDS(organoid, paste0(dir, out, '/organoid.rds')) 




# Step 2. Perform UMAP Clustering and Visualization ############################
### This is if you have not already performed SCTransform, or prefer not to use
### SCTransform due to wanting to define the parameters you use to normalize,
### FindVariableFeats, Scale, etc.

DefaultAssay(organoid) <- 'SoupX_counts'
organoid@active.assay
# Dimensionality reduction
organoid <- NormalizeData(organoid, verbose = TRUE)
organoid <- FindVariableFeatures(organoid, selection.method = "vst",
                                 nfeatures = 2000, verbose = TRUE)
organoid <- ScaleData(organoid, verbose = TRUE)

organoid <- RunPCA(organoid, features = VariableFeatures(object = organoid), 
                   verbose = TRUE)
organoid <- RunUMAP(organoid, dims = 1:30, verbose = TRUE)
organoid <- FindNeighbors(organoid, dims = 1:30, k.param = 30)

# Find clusters with different resolution values
resolutions <- seq(0.5, 1.5, by = 0.05)  # Adjust the range as needed
for (res in resolutions) {
  organoid <- FindClusters(organoid, resolution = res)
  num_clusters <- length(unique(organoid$seurat_clusters))
  cat("Resolution:", res, "- Number of clusters:", num_clusters, "\n")
  if (num_clusters == 24) {
    cat("Found 24 clusters with resolution:", res, "\n")
    break
  }
}

unique(organoid$seurat_clusters)

DimPlot(organoid, reduction = 'umap', group.by = 'seurat_clusters')


saveRDS(organoid, file = paste0(dir, '/outputs/04-Clustering_and_Batch_Correction/organoid1.rds'))
##### 07-07-2024 HERE
organoid <- readRDS(paste0(dir, '/outputs/04-Clustering_and_Batch_Correction/organoid1.rds'))

organoid@assays$RNA@layers$data

### Inspect the organoid object here.
# DimPlot(organoid, reduction = 'umap', group.by = 'seurat_clusters')

#### Step 3. Perform Harmony Reduction for Batch Effect Correction ############
library(harmony)

# Run Harmony for batch effect correction
organoid <- RunHarmony(organoid, group.by.vars = "sample", assay = "SoupX_counts")

# Find neighbors and clusters after Harmony integration
organoid <- FindNeighbors(organoid, dims = 1:30, k.param = 30, reduction = "harmony")

# Try different resolutions to find 24 clusters (being done to match previous iterations)
resolutions <- seq(0.55, 1.5, by = 0.05)  # Adjust the range as needed
for (res in resolutions) {
  organoid <- FindClusters(organoid, resolution = res)
  num_clusters <- length(unique(organoid$seurat_clusters))
  cat("Resolution:", res, "- Number of clusters:", num_clusters, "\n")
  if (num_clusters == 23 || num_clusters == 24) {
    cat("Found 23 clusters with resolution:", res, "\n")
    break
  }
}

# Proceed with UMAP visualization
organoid <- RunUMAP(organoid, reduction = "harmony", dims = 1:30, reduction.name = 'humap')
DimPlot(organoid, reduction = 'humap', group.by = "seurat_clusters")
DimPlot(organoid, reduction = 'umap', group.by = 'seurat_clusters')

DimPlot(organoid, reduction = 'humap', group.by = "seurat_clusters")
DimPlot(organoid, reduction = 'umap', group.by = 'seurat_clusters')


saveRDS(organoid, file = paste0(dir, '/outputs/05-Cluster_Annotation/organoid.rds'))
