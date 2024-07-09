library(dplyr)
library(patchwork)
library(ggplot2)
library(gridExtra)
library(Matrix)
library(scales)
library(cowplot)
library(Seurat) #Seurat v.5
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
out <- 'outputs/03-miQC_and_SoupX'





organoid <- readRDS(paste0(dir, '/outputs/03-miQC_and_SoupX/organoid3.rds'))



# Directories for raw counts of each sample
raw_counts_directories <- c("1_S16/af_quant", "2_S8/af_quant", "3_S9/af_quant", "4_S14/af_quant", "6_S23/af_quant",
                            "8_S17/af_quant", "9_S21/af_quant", "10_S7/af_quant", "11_S1/af_quant", "12_S20/af_quant", "13_S11/af_quant",
                            "14_S13/af_quant", "15_S5/af_quant", "16_S24/af_quant", "17_S19/af_quant", "18_S22/af_quant",
                            "19_S3/af_quant", "20_S15/af_quant", "21_S6/af_quant", "22_S2/af_quant", "23_S12/af_quant", "24_S10/af_quant",
                            "25_S18/af_quant", "26_S4/af_quant") # paths to quant files


gene_names <- fread('/projectnb/tcwlab/LabMember/adpelle1/projects/APOE_Jorganoid/outputs/CellRangerCount/sample_1/outs/filtered_feature_bc_matrix/features.tsv.gz', 
                    col.names = c('GENEID', 'SYMBOL', 'c3'))
# remove the c3 column from the gene_names data.table df
gene_names <- gene_names[, c3 := NULL]
# set GENEID to rownames attribute, in order to make the following easier
rownames(gene_names) <- gene_names$GENEID

# Load raw counts
custom_format <- list("counts" = c("U","S","A"))

# Function to apply SoupX correction
apply_soupX_correction <- function(seurat_obj, raw_counts_directory) {
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
  if (packageVersion("Seurat") >= "5") {
    seurat_counts <- GetAssayData(seurat_obj, layer = "counts")[common_genes, , drop = FALSE]
  } else {
    print('You need to update this code to work with an earlier (<=4) version of Seurat.')
  }
  # Ensure the gene order is the same in both matrices
  raw <- raw[order(rownames(raw)), , drop = FALSE]
  seurat_counts <- seurat_counts[order(rownames(seurat_counts)), , drop = FALSE]
  
  # Apply SoupX correction
  sc <- SoupChannel(tod = raw, toc = seurat_counts)
  sc <- setClusters(sc, seurat_obj$seurat_clusters)
  sc <- autoEstCont(sc, doPlot = FALSE)
  soup_out <- adjustCounts(sc, roundToInt = TRUE)
  
  # Integrate corrected counts into Seurat object
  seurat_obj[["SoupX_counts"]] <- CreateAssayObject(counts = soup_out)
  #seurat_obj <- SetAssayData(seurat_obj, layer = "counts", new.data = soup_out, assay = "RNA")
  
  
  #### Edit this so that there is a retained version of the original assay data
  # Remove the SoupX_counts assay
  # Return the corrected Seurat object
  return(seurat_obj)
}



mrna_counts <- data.frame(sample = character(), 
                          genotype = character(),
                          mrna_before_correction = numeric(), 
                          mrna_after_correction = numeric()
)

samples <- unique(organoid@meta.data$sample)
organoid_list <- list()

# Loop over each sample to apply SoupX correction and record mRNA sums
for (i in 1:length(samples)) {
  sample_id <- samples[i]
  raw_counts_directory_ending <- raw_counts_directories[i]
  raw_counts_directory <- paste0(dir, '/outputs/01-Simpleaf_Outputs/',raw_counts_directory_ending)
  # Subset Seurat object for the current sample
  indiv_seurat <- subset(organoid, subset = sample == sample_id)
  
  # Record the sum of mRNA molecules before correction
  mrna_before <- sum(indiv_seurat@meta.data$nCount_RNA)
  
  # Apply SoupX correction
  indiv_seurat <- apply_soupX_correction(indiv_seurat, raw_counts_directory)
  
  # Record the sum of mRNA molecules after correctione
  mrna_after <- sum(indiv_seurat@assays$RNA@layers$counts)
  
  # Print statement for debugging
  print(paste0("Sample: ", sample_id, ". mrna_before: ", mrna_before, ". mrna_after: ", mrna_after))
  
  # Append the corrected object onto the organoid_list
  organoid_list[[i]] <- indiv_seurat
  
  # Extract genotype for the current sample
  genotype <- unique(indiv_seurat@meta.data$genotype)
  
  # Explicitly set column names of the data frame before appending
  new_row <- data.frame(sample = sample_id,
                        genotype = genotype,
                        mrna_before_correction = mrna_before,
                        mrna_after_correction = mrna_after)
  colnames(new_row) <- colnames(mrna_counts)
  
  # Append the new row to mrna_counts
  mrna_counts <- rbind(mrna_counts, new_row)
  
}

print(paste0('length of the Seurat_obj list:', length(organoid_list), '.'))

count = 1
# merge all the Seurat objects into one.
tryCatch({
  if(length(organoid_list) > 1) {
    print(paste0('Merging function is processing #', count))
    # do.call to apply the merge function to a list of Seurat objects
    organoid4 <- do.call(merge, c(list(x = organoid_list[[1]], y = organoid_list[-1]), 
                                  list(project = "APOE_Organoid")))
    count = count + 1
  } else if(length(organoid_list) == 1) {
    # only one object, simply rename the project if necessary
    organoid4 <- RenameProject(organoid_list[[1]], "APOE_Organoid")
  } else {
    # error catch
    message("No data processed successfully.")
  }
}, error = function(e) {
  print(paste("Error:", e))
  # Additional error handling if required
})


print(mrna_counts)

# Aggregate mRNA counts by genotype
aggregated_mrna_counts <- aggregate(cbind(mrna_before_correction, mrna_after_correction) ~ genotype, data = mrna_counts, FUN = sum)

# Print the result
print(aggregated_mrna_counts)

write.csv(mrna_counts, file = paste0(dir, '/outputs/03-miQC_and_SoupX/SoupX_counts_by_Sample.csv'))
write.csv(aggregated_mrna_counts, file = paste0(dir, '/outputs/03-miQC_and_SoupX/SoupX_counts_by_Genotype.csv'))

saveRDS(organoid4, file = paste0(dir, '/outputs/03-miQC_and_SoupX/organoid4.rds'))