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

#merge the counts layers into one:
organoid2[["RNA"]] <- JoinLayers(organoid2[["RNA"]])

# SCTransform on the organoid object to do normalization, scaling, etc.
library(sctransform)
organoid3 <- SCTransform(organoid2, vars.to.regress = "percent.mt", verbose = FALSE)

# Dimensionality reduction
organoid3 <- RunPCA(organoid3, features = VariableFeatures(object = organoid3), verbose = TRUE)
organoid3 <- RunUMAP(organoid3, dims = 1:50)

# Clustering
organoid3 <- FindNeighbors(organoid3, dims = 1:50, k.param = 30, verbose = TRUE)
organoid3 <- FindClusters(organoid3, resolution = 0.6, verbose = TRUE)

saveRDS(organoid3, file = paste0(dir, '/outputs/03-miQC_and_SoupX/organoid3.rds'))
# organoid3 <- readRDS(paste0(dir, '/outputs/03-miQC_and_SoupX/organoid3.rds'))

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
  seurat_obj <- SetAssayData(seurat_obj, layer = "counts", new.data = soup_out, assay = "RNA")
  
  seurat_obj[["RNA"]] <- as(object = seurat_obj[["SoupX_counts"]], Class = "Assay5")
  # Remove the SoupX_counts assay
  seurat_obj[["SoupX_counts"]] <- NULL
  # Return the corrected Seurat object
  return(seurat_obj)
}


mrna_counts <- data.frame(sample = character(), 
                          genotype = character(),
                          mrna_before_correction = numeric(), 
                          mrna_after_correction = numeric()
)

samples <- unique(organoid3@meta.data$sample)
organoid_list <- list()
# Loop over each sample to apply SoupX correction and record mRNA sums
for (i in 1:length(samples)) {
  sample_id <- samples[i]
  raw_counts_directory_ending <- raw_counts_directories[i]
  raw_counts_directory <- paste0(dir, '/outputs/01-Simpleaf_Outputs/',raw_counts_directory_ending)
  # Subset Seurat object for the current sample
  indiv_seurat <- subset(organoid3, subset = sample == sample_id)
  
  # Record the sum of mRNA molecules before correction
  mrna_before <- sum(indiv_seurat@meta.data$nCount_RNA)
  
  # Apply SoupX correction
  indiv_seurat <- apply_soupX_correction(indiv_seurat, raw_counts_directory)
  
  # Record the sum of mRNA molecules after correction
  mrna_after <- sum(indiv_seurat@assays$SoupX_counts$counts)
  
  # Print statement for debugging
  print(paste0("Sample: ", sample_id, ". mrna_before: ", mrna_before, ". mrna_after: ", mrna_after))
  
  # Append the corrected object onto the organoid_list
  organoid_list[[sample_id]] <- indiv_seurat
  
  # Extract genotype for the current sample
  genotype <- unique(indiv_seurat@meta.data$genotype)
  
  # Append the mRNA count information to the data frame
  new_row <- data.frame(sample = sample_id,
                        genotype = genotype,
                        mrna_before_correction = mrna_before,
                        mrna_after_correction = mrna_after)
  colnames(new_row) <- colnames(mrna_counts)
  
  # Append the new row to mrna_counts
  mrna_counts <- rbind(mrna_counts, new_row)
}
count = 1
# Merge all corrected Seurat objects into one (all already have unique cell identifiers)
if(length(organoid_list) > 1) {
  
  print(paste0('Now at: ', count))
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

print(mrna_counts)

# Aggregate mRNA counts by genotype
aggregated_mrna_counts <- aggregate(cbind(mrna_before_correction, mrna_after_correction) ~ genotype, data = mrna_counts, FUN = sum)

# Print the result
print(aggregated_mrna_counts)


write.csv(mrna_counts, file = paste0(dir, '/outputs/03-miQC_and_SoupX/SoupX_counts_by_Sample.csv'))
write.csv(aggregated_mrna_counts, file = paste0(dir, '/outputs/03-miQC_and_SoupX/SoupX_counts_by_Genotype.csv'))

saveRDS(organoid4, file = paste0(dir, '/outputs/03-miQC_and_SoupX/organoid4.rds'))
