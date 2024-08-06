### Edited 07-30-2024 akg

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
### In this particular case, I am NOT using SCTransform, Instead, I am doing each
### step manually.

organoid <- readRDS(paste0(dir, out, '/organoid1.rds'))
## we're going to exclude sample17 aka APOE44_6 from further analysis:
unique(organoid@meta.data$sample)

library(dplyr)

organoid@meta.data <- organoid@meta.data %>%
  mutate(sample = case_when(
    # APOE33
    sample == "sample1" ~ "APOE33_1",
    sample == "sample14" ~ "APOE33_2",
    sample == "sample16" ~ "APOE33_3",
    sample == "sample19" ~ "APOE33_4",
    sample == "sample23" ~ "APOE33_5",
    sample == "sample8" ~ "APOE33_6",
    
    # APOE22
    sample == "sample12" ~ "APOE22_1",
    sample == "sample15" ~ "APOE22_2",
    sample == "sample21" ~ "APOE22_3",
    sample == "sample25" ~ "APOE22_4",
    sample == "sample3" ~ "APOE22_5",
    sample == "sample9" ~ "APOE22_6",
    
    # APOE44
    sample == "sample4" ~ "APOE44_1",
    sample == "sample10" ~ "APOE44_2",
    sample == "sample11" ~ "APOE44_3",
    sample == "sample22" ~ "APOE44_4",
    sample == "sample26" ~ "APOE44_5",
    sample == "sample17" ~ "APOE44_6",
    
    # APOE33Ch
    sample == "sample2" ~ "APOE33Ch_1",
    sample == "sample6" ~ "APOE33Ch_2",
    sample == "sample20" ~ "APOE33Ch_3",
    sample == "sample24" ~ "APOE33Ch_4",
    sample == "sample13" ~ "APOE33Ch_5",
    sample == "sample18" ~ "APOE33Ch_6",
    
    # Keep existing value if no match is found
    TRUE ~ sample  
  ))

# Check for duplicates in mapping
duplicated_samples <- organoid@meta.data %>% 
  group_by(sample) %>% 
  filter(n() > 1) %>% 
  pull(sample)

if (length(duplicated_samples) > 0) {
  print("Warning: Duplicate sample names found after renaming:")
  print(duplicated_samples)
} else {
  print("Renaming completed successfully without duplicates.")
}


### Plot the distribution of nFeature_RNA and nCount_RNA across samples:
library(ggplot2)
library(dplyr)
library(scales)
# Extract metadata
meta_data <- organoid@meta.data

# Identify any NA entries in the sample column
na_samples <- is.na(meta_data$sample)

# Remove NA entries if they exist
if (any(na_samples)) {
  message("Found NA entries in sample column, removing them.")
  meta_data <- meta_data[!na_samples, ]
}

# Define the expected order of samples
expected_samples <- c(
  paste0("APOE22_", 1:6),
  paste0("APOE33_", 1:6),
  paste0("APOE44_", 1:6),
  paste0("APOE33Ch_", 1:6)
)

# Ensure sample is a factor with the specified order
meta_data$sample <- factor(meta_data$sample, levels = expected_samples)

# Plot nFeature_RNA
# Remove non-positive values
meta_data <- meta_data[meta_data$nFeature_RNA > 0, ]
non_positive_features <- meta_data$nFeature_RNA <= 0
sum(non_positive_features)
# Recreate the plot p1
# Plot nFeature_RNA with explicit y-axis limits
y_breaks <- c(1e3, 1e4, 1e5)
y_labels <- scales::scientific(y_breaks)

# Plot nFeature_RNA with specific y-axis breaks and labels
p1 <- ggplot(meta_data, aes(x = sample, y = nFeature_RNA, fill = sample)) +
  geom_violin() +
  scale_y_log10(breaks = y_breaks, labels = y_labels) +
  theme_minimal() +
  ggtitle("nFeature_RNA Distribution by Sample (after QC)") +
  scale_fill_manual(values = scales::hue_pal()(length(expected_samples))) +
  theme(
    axis.title.x = element_text(size = 14, face = "bold"),   # x-axis label size and bold
    axis.title.y = element_text(size = 14, face = "bold"),   # y-axis label size and bold
    axis.text.x = element_text(size = 14, angle = 45, hjust = 1),  # x-axis tick label size
    axis.text.y = element_text(size = 14),    # y-axis tick label size
    legend.position = "none",
    plot.title = element_text(size = 16, hjust = 0.5, face = "bold")
  )


# Plot nCount_RNA
p2 <- ggplot(meta_data, aes(x = sample, y = nCount_RNA, fill = sample)) +
  geom_violin() +
  scale_y_log10() +
  theme_minimal() +
  ggtitle("nCount_RNA Distribution by Sample (after QC)") +
  scale_fill_manual(values = scales::hue_pal()(length(expected_samples))) +
  theme(
    axis.title.x = element_text(size = 14),   # x-axis label size
    axis.title.y = element_text(size = 14),   # y-axis label size
    axis.text.x = element_text(size = 14, angle = 45, hjust = 1),  # x-axis tick label size
    axis.text.y = element_text(size = 14),# y-axis tick label size
    legend.position = "none",
    plot.title = element_text(size = 16, hjust = 0.5, face = "bold")
  )

# Print the plots
print(p1)
print(p2)


# Plot nCount_RNA
p3 <- ggplot(meta_data, aes(x = sample, y = percent.mt, fill = sample)) +
  geom_violin() +
  scale_y_log10() +
  theme_minimal() +
  ggtitle("percent.mt Distribution by Sample (after QC)") +
  scale_fill_manual(values = scales::hue_pal()(length(expected_samples))) +
  theme(
    axis.title.x = element_text(size = 14),   # x-axis label size
    axis.title.y = element_text(size = 14),   # y-axis label size
    axis.text.x = element_text(size = 14, angle = 45, hjust = 1),  # x-axis tick label size
    axis.text.y = element_text(size = 14),# y-axis tick label size
    legend.position = "none",
    plot.title = element_text(size = 16, hjust = 0.5, face = "bold")
  )

print(p3)
organoid$DF.classifications_0.25_0.3_484
# Plot nCount_RNA
p3 <- ggplot(meta_data, aes(x = sample, y = percent.mt, fill = sample)) +
  geom_violin() +
  scale_y_log10() +
  theme_minimal() +
  ggtitle("percent.mt Distribution by Sample (after QC)") +
  scale_fill_manual(values = scales::hue_pal()(length(expected_samples))) +
  theme(
    axis.title.x = element_text(size = 14),   # x-axis label size
    axis.title.y = element_text(size = 14),   # y-axis label size
    axis.text.x = element_text(size = 14, angle = 45, hjust = 1),  # x-axis tick label size
    axis.text.y = element_text(size = 14),# y-axis tick label size
    legend.position = "none",
    plot.title = element_text(size = 16, hjust = 0.5, face = "bold")
  )








meta_data <- organoid@meta.data

# Filter the meta.data for cells with genotype 'APOE44' and count them
num_APOE44_cells <- sum(meta_data$genotype == 'APOE44')

# Print the result
num_APOE44_cells



unique(organoid@meta.data$sample)
# Remove observations corresponding to sample17 (aka APOE44_6)
organoid <- subset(organoid, subset = sample != "APOE44_6")
sum(organoid@meta.data$genotype == 'APOE44')

# saveRDS(organoid, file = paste0(dir, out, '/organoid2.rds')) #saved 07-30
organoid <- readRDS(paste0(dir, out, '/organoid2.rds'))

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

organoid <- FindClusters(organoid, resolution = 0.60)



## Optional: if you want to pre-define a certain number of clusters you can do this:
# Find clusters with different resolution values
# resolutions <- seq(0.5, 1.5, by = 0.05)  # Adjust the range as needed
# for (res in resolutions) {
#   organoid <- FindClusters(organoid, resolution = res)
#   num_clusters <- length(unique(organoid$seurat_clusters))
#   cat("Resolution:", res, "- Number of clusters:", num_clusters, "\n")
#   if (num_clusters == 24) {
#     cat("Found 24 clusters with resolution:", res, "\n")
#     break
#   }
# }

unique(organoid$seurat_clusters)

DimPlot(organoid, reduction = 'umap', group.by = 'seurat_clusters')


saveRDS(organoid, file = paste0(dir, '/outputs/04-Clustering_and_Batch_Correction/organoid1.rds'))
##### 07-07-2024 HERE
organoid <- readRDS(paste0(dir, '/outputs/04-Clustering_and_Batch_Correction/organoid1.rds'))

organoid@assays$RNA@layers$data

### Inspect the organoid object here.
# DimPlot(organoid, reduction = 'umap', group.by = 'seurat_clusters')



#### Step 3. Perform Harmony Reduction for Batch Effect Correction  ############
library(harmony)


# Run Harmony for batch effect correction
organoid <- RunHarmony(organoid, group.by.vars = "sample", assay = "SoupX_counts")

DefaultAssay(organoid)
# Find neighbors and clusters after Harmony integration
organoid <- FindNeighbors(organoid, dims = 1:30, k.param = 30, reduction = "harmony")

# Try different resolutions to find 24 clusters (being done to match previous iterations)
organoid <- FindClusters(organoid, resolution = 0.60)
# resolutions <- seq(0.55, 1.5, by = 0.05)  # Adjust the range as needed
# for (res in resolutions) {
#   organoid <- FindClusters(organoid, resolution = res)
#   num_clusters <- length(unique(organoid$seurat_clusters))
#   cat("Resolution:", res, "- Number of clusters:", num_clusters, "\n")
#   if (num_clusters == 23 || num_clusters == 24) {
#     cat("Found 23 clusters with resolution:", res, "\n")
#     break
#   }
# }

# Proceed with UMAP visualization
organoid <- RunUMAP(organoid, reduction = "harmony", dims = 1:30, reduction.name = 'humap')
DimPlot(organoid, reduction = 'humap', group.by = "seurat_clusters", label = FALSE)
DimPlot(organoid, reduction = 'umap', group.by = 'seurat_clusters')

DimPlot(organoid, reduction = 'humap', group.by = "seurat_clusters")
DimPlot(organoid, reduction = 'umap', group.by = 'seurat_clusters')


saveRDS(organoid, file = paste0(dir, '/outputs/05-Cluster_Annotation/organoid_unannotated.rds'))




