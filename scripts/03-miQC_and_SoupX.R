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
setwd("/projectnb/tcwlab/LabMember/akg/projects/APOE_Organoid")

dir <- "/projectnb/tcwlab/LabMember/akg/projects/APOE_Organoid"
out <- 'outputs/03-miQC_and_SoupX'

organoid <- readRDS('/projectnb/tcwlab/LabMember/akg/projects/APOE_Organoid/outputs/02-Quality-Control/organoid1.rds')

#Run miQC
#note: posterior cutoff = the posterior probability of a cell being part of the compromised distribution, a number between 0 and 1.
#Any cells below the appointed cutoff will be marked to keep. Defaults to 0.75.

organoid1 <- RunMiQC(organoid, percent.mt = "percent.mt", nFeature_RNA = "nFeature_RNA", posterior.cutoff = 0.97, model.slot = "flexmix_model") 

pdf(paste(dir, "/miQC/stringent_miQC_probability_plot.pdf", sep = ""), width = 15, height = 10)
PlotMiQC(organoid1, color.by = "miQC.probability") + ggplot2::scale_color_gradient(low = "grey", high = "purple")
dev.off()

pdf(paste(dir, "/miQC/miQC_keep_probability_plot.pdf", sep = ""), width = 15, height = 10)
PlotMiQC(organoid1, color.by = "miQC.keep")
dev.off()

plot1 <- FeatureScatter(organoid1, feature1 = "nCount_RNA", feature2 = "percent.mt") + theme(legend.position = 'none')
plot2 <- FeatureScatter(organoid1, feature1 = "nCount_RNA", feature2 = "nFeature_RNA") + theme(legend.position = 'none')
p <- plot1 + plot2
print(p)

# Perform MiQC Filtering
organoid2 <- subset(organoid1, miQC.keep == "keep") #118,688 cells and 36,601 genes

plot1 <- FeatureScatter(organoid2, feature1 = "nCount_RNA", feature2 = "percent.mt") + theme(legend.position = 'none')
plot2 <- FeatureScatter(organoid2, feature1 = "nCount_RNA", feature2 = "nFeature_RNA") + theme(legend.position = 'none')
p <- plot1 + plot2
print(p)

saveRDS(organoid2, file = paste0(dir, '/outputs/03-miQC_and_SoupX/organoid2.rds'))
organoid2 <- readRDS(paste0(dir, '/outputs/03-miQC_and_SoupX/organoid2.rds'))


# Step 2. SoupX Ambient RNA Correction ####
# This step involves 1. reading the alevin fry quant files for each sample
# then 2. running SoupX which identifies the 'ambient' RNAs, and marks them for 
# removal. Then 3. We update our Seurat object.

# In order to accomplish this, we'll split our organoid object by sample, and then
# run SoupX across each sample. We'll also update the rownames (gene names) in these
# sce objects in the same manner as in 01-Get_Seurat_Object.

## First, we need to cluster our Seurat object. This is necessary for SoupX
## But before we do this, we will integrate all our individual counts layers into one

# Function to merge individual Seurat objects based on sample IDs and combine their counts
# merge_seurat_objects_by_sample <- function(seurat_obj) {
#   # Extract unique sample IDs from meta.data
#   sample_ids <- unique(seurat_obj@meta.data$sample)
#   
#   # List to store individual counts matrices
#   counts_list <- list()
#   
#   # Loop through each sample ID, subset, and extract the counts matrix
#   for (sample_id in sample_ids) {
#     subset_obj <- subset(seurat_obj, subset = sample == sample_id)
#     counts_list[[sample_id]] <- GetAssayData(subset_obj, layer = "counts")
#   }
#   
#   # Get the union of all genes (rows) across the counts matrices
#   all_genes <- unique(unlist(lapply(counts_list, rownames)))
#   
#   # Align all counts matrices to have the same genes (rows)
#   aligned_counts_list <- lapply(counts_list, function(counts) {
#     aligned_counts <- Matrix::Matrix(0, nrow = length(all_genes), ncol = ncol(counts), sparse = TRUE)
#     rownames(aligned_counts) <- all_genes
#     colnames(aligned_counts) <- colnames(counts)
#     common_genes <- intersect(rownames(counts), all_genes)
#     aligned_counts[common_genes, ] <- counts[common_genes, ]
#     return(aligned_counts)
#   })
#   
#   # Combine the aligned counts matrices into one
#   combined_counts <- do.call(cbind, aligned_counts_list)
#   
#   # Create a new Seurat object with the combined counts matrix
#   new_seurat_obj <- CreateSeuratObject(counts = combined_counts, meta.data = seurat_obj@meta.data)
#   
#   # Normalize and find variable features for the new Seurat object
#   new_seurat_obj <- NormalizeData(new_seurat_obj)
#   new_seurat_obj <- FindVariableFeatures(new_seurat_obj)
#   new_seurat_obj <- ScaleData(new_seurat_obj)
#   
#   return(new_seurat_obj)
# }

# Apply the function to your combined Seurat object
# organoid2_merged <- merge_seurat_objects_by_sample(organoid2)


#merge the counts layers into one:



organoid2[["RNA"]] <- JoinLayers(organoid2[["RNA"]])




organoid2 <- NormalizeData(organoid2)
organoid2 <- FindVariableFeatures(organoid2)
organoid2 <- ScaleData(organoid2)
organoid2 <- RunPCA(organoid2)
organoid2 <- FindNeighbors(organoid2)
organoid2 <- FindClusters(organoid2, resolution = 0.5)

# saveRDS(organoid2, file = paste0(dir, '/outputs/03-miQC_and_SoupX/organoid3.rds'))

# Directories for raw counts of each sample
raw_counts_directories <- c("1_S16/af_quant", "2_S8/af_quant", "3_S9/af_quant", "4_S14/af_quant", "6_S23/af_quant",
                            "8_S17/af_quant", "9_S21/af_quant", "10_S7/af_quant", "11_S1/af_quant", "12_S20/af_quant", "13_S11/af_quant",
                            "14_S13/af_quant", "15_S5/af_quant", "16_S24/af_quant", "17_S19/af_quant", "18_S22/af_quant",
                            "19_S3/af_quant", "20_S15/af_quant", "21_S6/af_quant", "22_S2/af_quant", "23_S12/af_quant", "24_S10/af_quant",
                            "25_S18/af_quant", "26_S4/af_quant") # paths to quant files



# Load raw counts
custom_format <- list("counts" = c("U","S","A"))

# Function to apply SoupX correction
apply_soupX_correction <- function(seurat_obj, raw_counts_directory) {
  library(SoupX)
  library(fishpond)
  library(Seurat)
  
  # Load raw counts
  sce <- fishpond::loadFry(raw_counts_directory, outputFormat = "scRNA")
  
  # Map Ensembl IDs to gene symbols
  ensembl_ids <- rownames(sce)
  geneid_to_symbol <- setNames(gene_names$SYMBOL, gene_names$GENEID)
  symbol_names <- geneid_to_symbol[ensembl_ids]
  symbol_names <- ifelse(is.na(symbol_names), ensembl_ids, symbol_names)
  rownames(sce) <- make.unique(symbol_names)
  
  # Extract counts and intersect with Seurat object genes
  raw <- counts(sce)
  common_genes <- intersect(rownames(seurat_obj), rownames(raw))
  raw <- raw[common_genes, , drop = FALSE]
  
  # Access counts based on Seurat version
  if (packageVersion("Seurat") >= "4") {
    seurat_counts <- GetAssayData(seurat_obj, slot = "counts")[common_genes, , drop = FALSE]
  } else {
    seurat_counts <- seurat_obj@assays$RNA@counts[common_genes, , drop = FALSE]
  }
  
  # Ensure the gene order is the same in both matrices
  raw <- raw[order(rownames(raw)), , drop = FALSE]
  seurat_counts <- seurat_counts[order(rownames(seurat_counts)), , drop = FALSE]
  
  # Apply SoupX correction
  sc <- SoupChannel(tod = raw, toc = seurat_counts)
  sc <- setClusters(sc, seurat_obj$initial_seurat_clustering_0.5)
  sc <- autoEstCont(sc, doPlot = FALSE)
  soup_out <- adjustCounts(sc, roundToInt = TRUE)
  
  # Integrate corrected counts into Seurat object
  seurat_obj[["SoupX_counts"]] <- CreateAssayObject(counts = soup_out)
  seurat_obj <- SetAssayData(seurat_obj, slot = "counts", new.data = soup_out, assay = "RNA")
  
  # Return the corrected Seurat object
  return(seurat_obj)
}


mrna_counts <- data.frame(sample = character(), 
                          genotype = character(),
                          mrna_before_correction = numeric(), 
                          mrna_after_correction = numeric()
)

samples <- unique(organoid2@meta.data$samples)

# Loop over each sample to apply SoupX correction and record mRNA sums
for (i in 1:length(samples)) {
  sample_id <- samples[i]
  raw_counts_directory_ending <- raw_counts_directories[i]
  raw_counts_directory <- paste0(dir, '/outputs/01-Simpleaf_Outputs/',raw_counts_directory_ending)
  # Subset Seurat object for the current sample
  indiv_seurat <- subset(organoid2, subset = sample == sample_id)
  
  # Record the sum of mRNA molecules before correction
  mrna_before <- sum(indiv_seurat@meta.data$nCount_RNA)
  
  # Apply SoupX correction
  indiv_seurat <- apply_soupX_correction(indiv_seurat, raw_counts_directory)
  
  # Record the sum of mRNA molecules after correction
  mrna_after <- sum(indiv_seurat@assays$RNA@counts)
  
  # Print statement for debugging
  print(paste0("Sample: ", sample_id, ". mrna_before: ", mrna_before, ". mrna_after: ", mrna_after))
  
  # Append the corrected object onto the organoid_list
  organoid_list[[sample_id]] <- indiv_seurat
  
  # Extract genotype for the current sample
  genotype <- unique(indiv_seurat@meta.data$genotype)
  
  # Append the mRNA count information to the data frame
  mrna_counts <- rbind(mrna_counts, data.frame(sample = sample_id,
                                               genotype = genotype,
                                               mrna_before_correction = mrna_before,
                                               mrna_after_correction = mrna_after))
}

# Merge all corrected Seurat objects into one (all already have unique cell identifiers)
if(length(organoid_list) > 1) {
  # do.call to apply the merge function to a list of Seurat objects
  organoid3 <- do.call(merge, c(list(x = organoid_list[[1]], y = organoid_list[-1]), list(project = "APOE_Organoid")))
} else if(length(organoid_list) == 1) {
  # only one object, simply rename the project if necessary
  organoid3 <- RenameProject(organoid_list[[1]], "APOE_Organoid")
} else {
  # error catch
  message("No data processed successfully.")
}

print(mrna_counts)

# Aggregate mRNA counts by genotype
aggregated_mrna_counts <- aggregate(cbind(mrna_before_correction, mrna_after_correction) ~ genotype, data = mrna_counts, FUN = sum)

# Print the result
print(aggregated_mrna_counts)

