### 06-Get_Pseudobulk (last edit: 07-02-2024 by akg)

library(Seurat)
library(data.table)
library(tibble)
library(Matrix)


setwd('/projectnb/tcwlab/LabMember/akg/projects/APOE_Organoid/')
dir <- '/projectnb/tcwlab/LabMember/akg/projects/APOE_Organoid'
out <- '/outputs/06-Get_Pseudobulk'

organoid <- readRDS(paste0(dir,'/outputs/05-Cluster_Annotation/organoid1.rds'))



fp<-function(...)file.path(...)

#we create a list containing the original cell type name 
cellTypes <- c('Glutamatergic', 'GABAergic', 'Astrocyte', 'ImmatureAstrocyte', 'OPC', 
               'NPCs', 'VLMC', 'PericyteProgenitor', 'Pericyte', 'unknown', '23')
# cellType <- 'VLMC'
# dir w cellTypes .rds files
inputDir <- './outputs/05-Cluster_Annotation/'
outputDir <- './outputs/06-Get_Pseudobulk'

# seuratObject <- readRDS('./outputs/01-Get-Pseudobulk/Seurat/Astrocyte.rds')

seurat_Object <- readRDS(paste0(inputDir, 'unknown.rds'))
# for each cell type
# Load necessary libraries
library(Seurat)
library(Matrix)

# Define input directory and cell types
inputDir <- "path/to/input/directory/"
cellTypes <- c("cellType1", "cellType2", "cellType3")  # Replace with your actual cell types

# Load necessary libraries
library(Seurat)
library(Matrix)
library(data.table)

# Define input directory and cell types
inputDir <- "path/to/input/directory/"
outputDir <- "path/to/output/directory/"
cellTypes <- c("cellType1", "cellType2", "cellType3")  # Replace with your actual cell types

for (cellType in cellTypes) {
  # Seurat file for the current cell type
  seuratPath <- paste0(inputDir, cellType, '.rds')
  
  # Read Seurat file
  seuratObject <- readRDS(seuratPath)
  
  # Get assay data
  rna_data <- GetAssayData(seuratObject, assay = "RNA")
  
  # Calculate genes detected and UMI counts
  genesDet <- rowSums(rna_data)
  UMI_counts <- colSums(rna_data)
  
  # Ensure the cell names match
  seurat_cell_names <- colnames(seuratObject)
  umi_counts_names <- names(UMI_counts)
  
  # Print the first few cell names for inspection
  cat("Seurat cell names:\n")
  print(head(seurat_cell_names))
  cat("UMI counts names:\n")
  print(head(umi_counts_names))
  
  # Ensure the lengths match
  if (length(seurat_cell_names) != length(umi_counts_names)) {
    stop("Number of cells in Seurat object does not match the number in UMI counts.")
  }
  
  # Create a mapping between simplified and original names
  name_mapping <- data.frame(
    seurat_name = seurat_cell_names,
    umi_name = umi_counts_names,
    stringsAsFactors = FALSE
  )
  
  # Print the mapping for inspection
  cat("Name mapping:\n")
  print(head(name_mapping))
  
  # Reorder UMI_counts to match the Seurat object cell names
  UMI_counts <- UMI_counts[name_mapping$umi_name]
  
  # Add metadata to Seurat object
  seuratObject <- AddMetaData(seuratObject, metadata = UMI_counts, col.name = "median_UMIs")
  seuratObject$median_UMIs <- UMI_counts
  seuratObject@meta.data$median_genes_det <- genesDet[rownames(seuratObject@meta.data)]
  
  # Aggregate expression data
  pseudo_mat <- AggregateExpression(seuratObject, assays = 'RNA', group.by = c('sample'), return.seurat = FALSE)
  
  # Output aggregated expression into .csv
  pseudo_mat_df <- as.data.frame(pseudo_mat$RNA)
  pseudo_mat_df <- pseudo_mat_df[, colSums(is.na(pseudo_mat_df)) < nrow(pseudo_mat_df)]
  
  fwrite(data.table(pseudo_mat_df, keep.rownames = 'gene_id'), file = file.path(outputDir, paste0(cellType, '_pseudobulk.csv.gz')))
  
  # Calculate metrics for meta.data
  mtd <- seuratObject@meta.data
  setDT(mtd, keep.rownames = TRUE)
  
  mtd[, n.cells := .N, by = .(sample)]
  mtd[, n.cells.genotype := .N, by = .(genotype)]
  mtd[, genotype_total := .N, by = genotype]
  mtd[, prop.cells := n.cells / genotype_total, by = .(genotype)]
  mtd[, avg.pct.mt.per.cell := mean(`percent.mt`), by = c('sample')]
  mtd[, pass.threshold.n.cells := n.cells > 20]
  mtd[, outlier.n.cells := !pass.threshold.n.cells, by = .(sample)]
  
  mtsc <- unique(mtd, by = c('sample', 'genotype'))
  mtsc <- mtsc[, c("sample", "genotype", "median_UMIs", "n.cells", "pass.threshold.n.cells", "avg.pct.mt.per.cell", "prop.cells")]
  
  # Output sample level metadata into .csv
  fwrite(mtsc, file.path(outputDir, paste0(cellType, '_sample_level_metadata.csv.gz')), row.names = TRUE)
}




# # humap by genotype function
# plot_genotype <- function(seurat_object, genotype, point_size = 1, alpha_other = 0.1) {
#   # Create a new column for highlighting
#   seurat_object[["highlight"]] <- ifelse(seurat_object@meta.data$genotype == genotype, genotype, "Other")
#   
#   # Define colors with vibrant highlight and more opaque others
#   highlight_colors <- c("Other" = scales::alpha("lightgray", alpha_other), genotype = "dodgerblue")
#   
#   # Adjust point sizes dynamically
#   point_sizes <- ifelse(seurat_object@meta.data$genotype == genotype, 3, 1)  # Larger size for highlighted genotype
#   
#   # Generate the UMAP plot with dynamic point sizes
#   p <- DimPlot(seurat_object, reduction = "humap", group.by = "highlight", cols = highlight_colors, pt.size = point_sizes) +
#     ggtitle(paste("hUMAP Highlighting", genotype))
#   
#   return(p)
# }
# 
# plot_genotype(organoid, 'APOE33Ch')

DimPlot(organoid, reduction = 'humap', group.by = 'genotype')





