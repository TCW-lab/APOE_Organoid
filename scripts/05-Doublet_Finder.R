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
library(sctransform)
library(DoubletFinder)

setwd("/projectnb/tcwlab/LabMember/akg/projects/APOE_Organoid")

dir <- "/projectnb/tcwlab/LabMember/akg/projects/APOE_Organoid"
out <- '/outputs/05-Doublet_Finder'


# Load the Seurat object
organoid <- readRDS(paste0(dir, '/outputs/05-Doublet_Finder/organoid.rds'))

## The following code is from 
## Chris McGinnis (UCSF): https://github.com/chris-mcginnis-ucsf/DoubletFinder





# Split the Seurat object into smaller subsets
subset_list <- SplitObject(organoid, split.by = "RNA_snn_res.0.6")

# Save each subset to a file
for (i in 1:length(subset_list)) {
  saveRDS(subset_list[[i]], file = paste0(dir, '/outputs/05-Doublet_Finder/subset_', i, '.rds'))
}

# Function to run DoubletFinder on a subset
run_doublet_finder <- function(subset_index) {
  # Load the subset
  subset <- readRDS(paste0(dir, out, "/subset_", subset_index, ".rds"))
  
  # DoubletFinder pipeline
  ## pK Identification (no ground-truth) ---------------------------------------
  sweep.res.list_subset <- paramSweep(subset, PCs = 1:10, sct = FALSE)
  sweep.stats_subset <- summarizeSweep(sweep.res.list_subset, GT = FALSE)
  bcmvn_subset <- find.pK(sweep.stats_subset)
  
  ## Homotypic Doublet Proportion Estimate -------------------------------------
  annotations <- subset@meta.data$seurat_clusters
  homotypic.prop <- modelHomotypic(annotations)
  nExp_poi <- round(0.075 * nrow(subset@meta.data))  ## Assuming 7.5% doublet formation rate
  nExp_poi.adj <- round(nExp_poi * (1 - homotypic.prop))
  
  ## Run DoubletFinder with varying classification stringencies ----------------
  subset <- doubletFinder_v3(subset, PCs = 1:10, pN = 0.25, pK = 0.09, 
                             nExp = nExp_poi.adj, reuse.pANN = FALSE, sct = FALSE)
  
  # Save the processed subset
  saveRDS(subset, paste0(dir, out, "/processed_subset_", subset_index, ".rds"))
}

# Run DoubletFinder on each subset
for (i in 1:length(subset_list)) {
  run_doublet_finder(i)
}

# Load all processed subsets
processed_list <- list()
for (i in 1:length(subset_list)) {
  processed_list[[i]] <- readRDS(paste0("processed_subset_", i, ".rds"))
}

# Merge the processed subsets into a single Seurat object
combined_organoid <- Reduce(function(x, y) merge(x, y), processed_list)

# Save the combined Seurat object
saveRDS(combined_organoid, file = paste0(dir, out, '/organoid.rds'))
