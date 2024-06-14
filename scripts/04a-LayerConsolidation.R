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



# extract individual counts and data layers
individual_count_layers <- lapply(paste0("counts.", 1:24), function(layer_name) organoid@assays$RNA@layers[[layer_name]])
individual_data_layers <- lapply(paste0("data.", 1:24), function(layer_name) organoid@assays$RNA@layers[[layer_name]])

# merge individual counts and data layers into a single layer each
combined_counts <- JoinLayers(individual_count_layers)
combined_data <- JoinLayers(individual_data_layers)

# update counts and data layers in the Seurat object
organoid@assays$RNA@layers$counts <- combined_counts
organoid@assays$RNA@layers$data <- combined_data

# remove individual count and data layers
layers_to_remove <- c(paste0("counts.", 1:24), paste0("data.", 1:24))
organoid@assays$RNA@layers <- organoid@assays$RNA@layers[!names(organoid@assays$RNA@layers) %in% layers_to_remove]



saveRDS(organoid, file = paste0(dir,out,'/organoid.rds'))