#### 04-08-2024 Script 
library(Seurat)
library(data.table)
library(dplyr)
library(Matrix)

setwd('/projectnb/tcwlab/LabMember/akg/APOE_Jorganoid_project/')


#Idents:
# c("GABAergic", "GABAergic", "Glutamatergic", "Astrocyte", "Astrocyte", "Glutamatergic", "OPC/Oligo", "GABAergic", "GABAergic", "Glutamatergic",
#   "NPCs - cycling", "NPCs - cycling", "Astrocyte", "Astrocyte", "Astrocyte", "Astrocyte", "Glutamatergic", "Astrocyte", "NPCs - cycling", "VLMC",
#   20, "NPCs - cycling", "Glutamatergic", "Pigmented Epithelial", 24)
# organoid <- RenameIdents(organoid, "5" = "Glutamatergic")

# Save the annotated .rds object #
# saveRDS(organoid, 'outputs/organoid_annotated.rds')
organoid <- readRDS('outputs/organoid_annotated.rds')
organoid0 <- readRDS('/projectnb/tcwlab/LabMember/akg/APOE_Jorganoid_project/outputs/miQC/organoid0.rds')

# 
# organoid@meta.data$orig.ident
# 
# organoid0@meta.data$percent.mt <- PercentageFeatureSet(organoid0, features = mt_gene_names)
# organoid0@meta.data$percent.mt



# Filter genes expressed in less than 3 cells
raw_counts_mat <- organoid@assays$RNA@counts
genes_to_keep <- colSums(raw_counts_mat != 0) >= 3
organoid <- subset(organoid, features = rownames(organoid)[genes_to_keep])

# Calculate upper limit, set at 99.5 percentile
upper_limit <- quantile(organoid@meta.data$nCount_RNA, probs = 0.995)

# Here I am assuming you might have a method to deal with thresholding, as CalculateBarcodeInflections and SubsetByBarcodeInflections 
# do not exist by default in Seurat. You would typically use a custom function or adjust based on other Seurat functions:
# This example sets a simpler approach using subset
organoid <- subset(organoid, subset = nCount_RNA <= upper_limit)

# Optionally, also apply a lower threshold
organoid <- subset(organoid, subset = nCount_RNA >= 1000)
organoid <- subset(organoid, subset = percent.mt < 5)
organoid <- subset(organoid, subset = nFeature_RNA > 200) #number of genes detected in each cell
#7,669 cells and 24,255 genes

#Only include cells with at least 500 molecules expressed
organoid <- subset(organoid, subset = nCount_RNA > 500)


# For QC metrics:
vps <- VlnPlot(object = organoid, 
               features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), 
               ncol = 3, 
               slot = "counts", 
               pt.size = 0, 
               group.by = 'orig.ident')  # Forces all data into a single group



print(vps)



# humap by genotype function
plot_genotype <- function(seurat_object, genotype, point_size = 1, alpha_other = 0.1) {
  # Create a new column for highlighting
  seurat_object[["highlight"]] <- ifelse(seurat_object@meta.data$genotype == genotype, genotype, "Other")
  
  # Define colors with vibrant highlight and more opaque others
  highlight_colors <- c("Other" = scales::alpha("lightgray", alpha_other), genotype = "dodgerblue")
  
  # Adjust point sizes dynamically
  point_sizes <- ifelse(seurat_object@meta.data$genotype == genotype, 3, 1)  # Larger size for highlighted genotype
  
  # Generate the UMAP plot with dynamic point sizes
  p <- DimPlot(seurat_object, reduction = "humap", group.by = "highlight", cols = highlight_colors, pt.size = point_sizes) +
    ggtitle(paste("hUMAP Highlighting", genotype))
  
  return(p)
}




# List of genotypes
umap_data <- FetchData(organoid, vars = c("humap_1", "humap_2", "genotype"))

# Iterate over genotypes and create plots
plots <- lapply(genotypes, function(specific_genotype) {
  # Create a new column for highlighting in the DataFrame
  umap_data$highlight <- ifelse(umap_data$genotype == specific_genotype, "highlighted", "other")
  
  # Define colors and sizes based on the highlight status
  umap_data$color <- ifelse(umap_data$highlight == "highlighted", "steelblue3", "lightgray")
  umap_data$size <- ifelse(umap_data$highlight == "highlighted", 1, 0.75)  # Adjust sizes to be more distinct
  
  # Plot using ggplot
  p <- ggplot(umap_data, aes(x = humap_1, y = humap_2, color = color, size = size)) +
    geom_point(alpha = ifelse(umap_data$highlight == "highlighted", 1, 0.2)) +  # Correct alpha and make non-highlighted more transparent
    scale_color_identity() +
    scale_size_identity() +
    labs(title = paste("hUMAP Plot Highlighting", specific_genotype)) +
    theme_minimal() +
    theme(legend.position = "none",  # Remove legend to clean up the plot
          axis.title = element_blank(),  # Remove axis titles
          axis.text = element_blank(),  # Remove axis text
          axis.ticks = element_blank())  # Remove axis ticks
  
  # Return the plot
  return(p)
})


if ("patchwork" %in% rownames(installed.packages())) {
  library(patchwork)
  plot_grid <- wrap_plots(plots, ncol = 2)
  print(plot_grid)
} else {
  lapply(plots, print)
}

library(dplyr)



# Ensure you've loaded the required libraries
library(Seurat)
library(dplyr)
library(ggplot2)

# Convert the metadata to a data frame
meta_data_df <- as.data.frame(organoid@meta.data)

# Confirm that 'old.ident' and 'genotype' are column names in the dataframe
print(colnames(meta_data_df))

# Now perform the operation with the data frame
genotype_by_cluster <- meta_data_df %>%
  dplyr::count(old.ident, genotype) %>%
  dplyr::arrange(old.ident, genotype)

# Check the resulting data frame to ensure it's as expected
print(head(genotype_by_cluster))


# Proceed with plotting
p <- ggplot(genotype_by_cluster, aes(x = old.ident, y = n, fill = genotype)) +
  geom_col(position = position_dodge()) +
  scale_fill_brewer(palette = "Set1") +
  labs(title = "Genotype Counts Across Clusters", x = "Cluster", y = "Count") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        legend.position = "bottom")

# Print the plot
print(p)







unique(organoid$sample)
# Double-check to see our annotations:
DimPlot(organoid, reduction = "humap", group.by = "ident", label = TRUE)
colnames(organoid@assays$SoupX_counts)

###### Displaying the proportion of each genotype beneath each cluster name ####

library(ggplot2)
library(ggtext)

organoid[["old.ident"]] <- Idents(object = organoid)
organoid <- StashIdent(object = organoid, save.name = "old.ident")

# Calculate proportions of each genotype within each cluster
genotype_proportions <- organoid@meta.data %>%
  group_by(old.ident, genotype) %>%
  summarise(count = n(), .groups = 'drop') %>%
  group_by(old.ident) %>%
  mutate(total = sum(count), proportion = count / total) %>%
  summarise(label = paste(old.ident, paste(genotype, sprintf("%.2f%%", proportion * 100), sep=": ", collapse=",\n"), sep="\n"), .groups = 'drop')


# Check the resulting labels
print(genotype_proportions)



# plot the UMAP
p <- DimPlot(organoid, reduction = "humap", group.by = "old.ident", label = FALSE)

# get UMAP coordinates
humap_coordinates <- Embeddings(organoid, "humap")

# convert to a data frame
humap_df <- data.frame(umap_1 = humap_coordinates[,1], umap_2 = humap_coordinates[,2])

# add old.ident information for grouping
humap_df$old.ident <- organoid@meta.data$old.ident

# calculate centroids for each cluster
cluster_centroids <- humap_df %>%
  group_by(old.ident) %>%
  summarise(umap_1 = mean(umap_1), umap_2 = mean(umap_2))


# Merge centroids with labels for direct plotting
plot_data <- merge(genotype_proportions, cluster_centroids, by = "old.ident")

# Add labels to the plot
p <- DimPlot(organoid, reduction = "humap", group.by = "old.ident", label = FALSE) +
  geom_text(data = plot_data, aes(x = umap_1, y = umap_2, label = label),
            hjust = 0.5, vjust = 0.5,
            check_overlap = FALSE, size = 3.5, fontface = "bold") +
            ggtitle("Proportion of Genotypes by hUMAP Cluster")

# Print the plot
print(p)



#### Subsetting each cell type and creating a separate Seurat object from that #
########## Astrocyte ###########################################################
# subset to only include GABAergic cells
Astrocyte <- subset(organoid, idents = "Astrocyte")



saveRDS(Astrocyte, 'outputs/01-Get-Pseudobulk/Astrocyte.rds')

Astrocyte@meta.data$
############# GABAergic ########################################################
# subset to only include GABAergic cells
GABAergic <- subset(organoid, idents = "GABAergic")
saveRDS(GABAergic, 'outputs/01-Get-Pseudobulk/GABAergic.rds')


################## Glutamatergic ###############################################
# subset to only include glutamatergic cells
glutamatergic <- subset(organoid, idents = "Glutamatergic")
saveRDS(glutamatergic, 'outputs/01-Get-Pseudobulk/glutamatergic.rds')

###################### OPC/Oligo ###############################################
# subset to only include OPC cells
OPC <- subset(organoid, idents = "OPC")
saveRDS(OPC, 'outputs/01-Get-Pseudobulk/OPC.rds')

######################### NPCs (cycling) #######################################
# subset to only include NPCs
NPCs <- subset(organoid, idents = "NPCs")
saveRDS(NPCs, 'outputs/01-Get-Pseudobulk/NPCs.rds')

################################### VLMC #######################################
# subset to only include VLMC
VLMC <- subset(organoid, idents = "VLMC")
saveRDS(VLMC, 'outputs/01-Get-Pseudobulk/VLMC.rds')

########################################## Pig. Epi. ###########################
 
# subset to only include Pigmented epithelial
Epithelial <- subset(organoid, idents = "Epithelial")
saveRDS(Epithelial, 'outputs/01-Get-Pseudobulk/Epithelial.rds')

############################################# Cluster 20 #######################
# subset to only include Cluster 20 (unknown)
unknown <- subset(organoid, idents = "20")
saveRDS(unknown, 'outputs/01-Get-Pseudobulk/unknown.rds')




