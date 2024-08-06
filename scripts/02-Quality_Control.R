## 02-Quality_Control.R
## Adapted from a script written by Deepti Murthy for the TCW Lab.

#In this analysis, I remove cells with high mitochondrial percentage and 
## remove outliers, scale, normalize, and cluster.
###
#R version 4.3.1
sessionInfo()

#load libraries 
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
library(flexmix)
library(miQC)
library(singleCellTK)
library(fgsea)
library(singleCellTK)
library(DropletUtils)
library(SeuratData)
library(tidyr)
library(DoubletFinder)
library(sctransform)
library(SeuratWrappers)

# .libPaths()
# source("~/.Rprofile")

# Step 1: Read in the indiv samples' Seurat object files & combine into one ####

###* denotes changes made by AKG
###* Then, combine them into one Seurat object. This object will be used to
###* be a point of comparison for before/after quality control steps
setwd("/projectnb/tcwlab/LabMember/akg/projects/APOE_Organoid")

dir <- "/projectnb/tcwlab/LabMember/akg/projects/APOE_Organoid/"
out<-'outputs/02-Quality_Control/'


S1 <- readRDS(file = "./outputs/01-individual-Seurat-Files/S1.rds")
S2 <- readRDS(file = "./outputs/01-individual-Seurat-Files/S2.rds")
S3 <- readRDS(file = "./outputs/01-individual-Seurat-Files/S3.rds")
S4 <- readRDS(file = "./outputs/01-individual-Seurat-Files/S4.rds")
###* No S 5 (deemed poor quality by researchers)
S6 <- readRDS(file = "./outputs/01-individual-Seurat-Files/S6.rds")
###* No S 7 (deemed poor quality by researchers)
S8 <- readRDS(file = "./outputs/01-individual-Seurat-Files/S8.rds")
S9 <- readRDS(file = "./outputs/01-individual-Seurat-Files/S9.rds")
S10 <- readRDS(file = "./outputs/01-individual-Seurat-Files/S10.rds")
S11 <- readRDS(file = "./outputs/01-individual-Seurat-Files/S11.rds")
S12 <- readRDS(file = "./outputs/01-individual-Seurat-Files/S12.rds")
S13 <- readRDS(file = "./outputs/01-individual-Seurat-Files/S13.rds")
S14 <- readRDS(file = "./outputs/01-individual-Seurat-Files/S14.rds")
S15 <- readRDS(file = "./outputs/01-individual-Seurat-Files/S15.rds")
S16 <- readRDS(file = "./outputs/01-individual-Seurat-Files/S16.rds")
S17 <- readRDS(file = "./outputs/01-individual-Seurat-Files/S17.rds")
S18 <- readRDS(file = "./outputs/01-individual-Seurat-Files/S18.rds")
S19 <- readRDS(file = "./outputs/01-individual-Seurat-Files/S19.rds")
S20 <- readRDS(file = "./outputs/01-individual-Seurat-Files/S20.rds")
S21 <- readRDS(file = "./outputs/01-individual-Seurat-Files/S21.rds")
S22 <- readRDS(file = "./outputs/01-individual-Seurat-Files/S22.rds")
S23 <- readRDS(file = "./outputs/01-individual-Seurat-Files/S23.rds")
S24 <- readRDS(file = "./outputs/01-individual-Seurat-Files/S24.rds")
S25 <- readRDS(file = "./outputs/01-individual-Seurat-Files/S25.rds")
S26 <- readRDS(file = "./outputs/01-individual-Seurat-Files/S26.rds")


# Load the Seurat objects into a list with names
samples_list <- list(
  S1 = S1, S2 = S2, S3 = S3, S4 = S4, S6 = S6, S8 = S8, 
  S9 = S9, S10 = S10, S11 = S11, S12 = S12, S13 = S13, 
  S14 = S14, S15 = S15, S16 = S16, S17 = S17, S18 = S18, 
  S19 = S19, S20 = S20, S21 = S21, S22 = S22, S23 = S23, 
  S24 = S24, S25 = S25, S26 = S26
)

# Define the groupings and corresponding sample names
genotype_groups <- list(
  APOE22 = c("S3", "S9", "S12", "S15", "S21", "S25"),
  APOE33 = c("S1", "S8", "S14", "S16", "S19", "S23"),
  APOE33Ch = c("S2", "S6", "S13", "S18", "S20", "S24"),
  APOE44 = c("S4", "S10", "S11", "S17", "S22", "S26")
)

sample_names_groups <- list(
  APOE22 = c("APOE22_5", "APOE22_6", "APOE22_1", "APOE22_2", "APOE22_3", "APOE22_4"),
  APOE33 = c("APOE33_1", "APOE33_6", "APOE33_2", "APOE33_3", "APOE33_4", "APOE33_5"),
  APOE33Ch = c("APOE33Ch_1", "APOE33Ch_2", "APOE33Ch_5", "APOE33Ch_6", "APOE33Ch_3", "APOE33Ch_4"),
  APOE44 = c("APOE44_1", "APOE44_2", "APOE44_3", "APOE44_6", "APOE44_4", "APOE44_5")
)

# Loop through each genotype group and assign metadata
for (genotype in names(genotype_groups)) {
  sample_names <- sample_names_groups[[genotype]]
  sample_objects <- genotype_groups[[genotype]]
  
  for (i in seq_along(sample_objects)) {
    sample_name <- sample_names[i]
    sample_key <- sample_objects[i]
    
    # Check if sample_obj exists in the list
    if (!sample_key %in% names(samples_list)) {
      warning(paste("Object", sample_key, "not found in list. Skipping."))
      next
    }
    
    # Get the Seurat object
    sample_obj <- samples_list[[sample_key]]
    
    # Add metadata
    sample_obj@meta.data$sample <- sample_name
    sample_obj@meta.data$genotype <- genotype
    
    # Update the list with the modified object
    samples_list[[sample_key]] <- sample_obj
  }
}


# # Check that the metadata was assigned correctly
# S1@meta.data[c("sample", "genotype")]
# S17@meta.data[c("sample", "genotype")]
# Repeat for other objects as needed


################################################################################


mt_gene_names <- c(
  "ENSG00000198888", # ND1
  "ENSG00000198763", # ND2
  "ENSG00000198804", # COX1
  "ENSG00000198712", # COX2
  "ENSG00000228253", # ATP8
  "ENSG00000198899", # ATP6
  "ENSG00000198938", # COX3
  "ENSG00000198840", # ND3
  "ENSG00000212907", # ND4L
  "ENSG00000198886", # ND4
  "ENSG00000198786", # ND5
  "ENSG00000198695", # ND6
  "ENSG00000198727", # CYTB
  "ENSG00000210049", # MT-TF
  "ENSG00000210082", # MT-TV
  "ENSG00000209082", # MT-TL1
  "ENSG00000198888", # MT-TR
  "ENSG00000210100", # MT-TN
  "ENSG00000210107", # MT-TG
  "ENSG00000210112", # MT-TL2
  "ENSG00000210117", # MT-TS1
  "ENSG00000210127", # MT-TV
  "ENSG00000210135", # MT-TE
  "ENSG00000210144", # MT-TS2
  "ENSG00000210151", # MT-TH
  "ENSG00000210154", # MT-TD
  "ENSG00000210156", # MT-TK
  "ENSG00000210164", # MT-TM
  "ENSG00000210174", # MT-TI
  "ENSG00000210184", # MT-TT
  "ENSG00000210191", # MT-TW
  "ENSG00000210194", # MT-TC
  "ENSG00000210211", # MT-TY
  "ENSG00000210228", # MT-TA
  "ENSG00000210243", # MT-TQ
  "ENSG00000198899", # MT-RNR1 (12S rRNA)
  "ENSG00000198763",  # MT-RNR2 (16S rRNA)
  "ND1", "ND2", "COX1", "COX2", "ATP8", "ATP6", "COX3", 
  "ND3", "ND4L", "ND4", "ND5", "ND6", "CYTB",
  "TF", "TV", "TL1", "TR", "TN", "TG", "TL2", "TS1", "TV",
  "TE", "TS2", "TH", "TD", "TK", "TM", "TI", "TT", "TW", 
  "TC", "TY", "TA", "TQ", "RNR1", "RNR2",
  "MT-ND4L", "MT-CYB", "MT-ND2", "MT-ATP6","MT-CO2", "MT-ND5", "MT-CO1", "MT-ATP8", "MT-ND6", "MT-CO3", "MT-ND4", "MT-ND1" 
)

lower_nCount_RNA_cutoffs <- c(
  1000, 1000, 1000, 1000, 1000, 1000, 900, 1000, 1200, 1000, 
  1000, 1500, 1000, 1000, 1000, 1000, 1500, 1000, 900, 1200, 
  1200, 1500, 1400, 1400
)


# Initialize a data frame to store cell counts before and after QC
cell_counts_DF <- data.frame(
  sample = character(), 
  genotype = character(),
  cells_before_miQC = integer(),
  cells_after_miQC = integer(),
  percent_recovery_miQC = numeric(),
  percent_removed_miQC = numeric(),
  cells_before_QC = integer(), 
  cells_after_QC = integer(),
  percent_recovery_QC = numeric(),
  percent_removed_QC = numeric(),
  cells_after_doubletF = integer(),
  percent_recovery_DF = numeric(),
  percent_removed_DF = numeric(),
  stringsAsFactors = FALSE
)

# Initialize an empty list to store sample's individual RDS files
organoid_list <- list()

# Directory to save miQC plots
img_dir <- 'outputs/images/02-Quality_Control/'

# Loop through each sample
for (i in 1:length(samples_list)) { 
  cat("Processing sample", i, "\n")
  indiv_Seurat <- samples_list[[i]]
  
  sample_id <- unique(indiv_Seurat@meta.data$sample)
  genotype <- unique(indiv_Seurat@meta.data$genotype)
  
  cells_before_QC <- ncol(indiv_Seurat)
  cat("Sample ID:", sample_id, "- Cells before QC:", cells_before_QC, "\n")
  
  # Subset mt_gene_names to include only genes present in the Seurat object
  counts_data <- GetAssayData(indiv_Seurat, layer = "counts", assay = "RNA")
  valid_mt_gene_names <- mt_gene_names[mt_gene_names %in% rownames(counts_data)]
  additional_mt_genes <- rownames(counts_data)[grepl("^mt-", rownames(counts_data), ignore.case = TRUE)]
  valid_mt_gene_names <- unique(c(valid_mt_gene_names, additional_mt_genes))
  indiv_Seurat[["percent.mt"]] <- PercentageFeatureSet(indiv_Seurat, features = valid_mt_gene_names)
  
  # Additionally label the ribosomal percentage
  indiv_Seurat[["percent.ribo"]] <- PercentageFeatureSet(indiv_Seurat, pattern = "^RP[SL]")
  
  # Run miQC as the first filtering step
  cells_before_miQC <- ncol(indiv_Seurat)
  indiv_Seurat <- RunMiQC(indiv_Seurat, percent.mt = "percent.mt", nFeature_RNA = "nFeature_RNA", 
                          posterior.cutoff = 0.75, model.slot = "flexmix_model")
  
  # Create miQC probability plot
  miQC_probability_plot <- PlotMiQC(indiv_Seurat, color.by = "miQC.probability") +
    ggplot2::scale_color_gradient(low = "grey", high = "purple")
  
  # Create miQC keep probability plot
  miQC_keep_probability_plot <- PlotMiQC(indiv_Seurat, color.by = "miQC.keep")
  # Combine the miQC plots
  miQC_combined_plot <- gridExtra::grid.arrange(miQC_probability_plot, miQC_keep_probability_plot, ncol = 2)
  
  # Save the combined miQC plots
  ggsave(paste0(img_dir, sample_id, "_miQC_combined_plots.png"), plot = miQC_combined_plot, 
         width = 15, height = 10, units = "in", dpi = 300)
  
  plot1 <- FeatureScatter(indiv_Seurat, feature1 = "nCount_RNA", feature2 = "percent.mt") + theme(legend.position = 'none') 
  plot2 <- FeatureScatter(indiv_Seurat, feature1 = "nCount_RNA", feature2 = "nFeature_RNA") + theme(legend.position = 'none') 
  p <- plot1 + plot2 + plot_annotation(title = paste(sample_id))
  ggsave(paste0(img_dir, sample_id, "_FeatureScatter_plots.png"), plot = p, width = 15, height = 10, units = "in", dpi = 300)
  
  ## Subset based on miQC filtering criteria
  indiv_Seurat <- subset(indiv_Seurat, subset = miQC.probability < 0.97)
  cells_after_miQC <- ncol(indiv_Seurat)
  percent_recovery_miQC <- 100 * (cells_after_miQC / cells_before_miQC)
  percent_removed_miQC <- 100 - percent_recovery_miQC

  
  # Proceed with additional QC steps
  lower_cell_cutoff <- lower_nCount_RNA_cutoffs[i]
  cat("Lower cell cutoff:", lower_cell_cutoff, "\n")
  sorted_counts <- sort(indiv_Seurat@meta.data$nCount_RNA, decreasing = FALSE)
  upper_cell_cutoff_rank <- ceiling(0.99 * length(sorted_counts))
  cat("Upper cell cutoff rank:", upper_cell_cutoff_rank, "\n")
  
  # Subset the Seurat object based on nCount_RNA
  indiv_Seurat <- subset(indiv_Seurat, subset = nCount_RNA > lower_cell_cutoff)
  
  # Subset for nFeature_RNA
  nFeature_RNA_values <- indiv_Seurat@meta.data$nFeature_RNA
  sorted_nFeature_RNA <- sort(nFeature_RNA_values, decreasing = FALSE)
  upper_index_nFeature <- ceiling(0.99 * length(sorted_nFeature_RNA))
  upper_feature_cutoff <- sorted_nFeature_RNA[upper_index_nFeature]
  cat("Upper feature cutoff:", upper_feature_cutoff, "\n")
  
  lower_limit_nFeature <- 200
  indiv_Seurat <- subset(indiv_Seurat, subset = nFeature_RNA > lower_limit_nFeature)
  
  # Recalculate number of cells after QC
  cells_after_QC <- ncol(indiv_Seurat)
  percent_recovery_QC <- 100 * (cells_after_QC / cells_before_QC)
  percent_removed_QC <- 100 - percent_recovery_QC
  cat("Cells after additional QC:", cells_after_QC, "\n")
  
  # Normalize and find variable features
  indiv_Seurat <- NormalizeData(indiv_Seurat)
  indiv_Seurat <- FindVariableFeatures(indiv_Seurat, selection.method = "vst", nfeatures = 2000)
  indiv_Seurat <- ScaleData(indiv_Seurat)
  indiv_Seurat <- RunPCA(indiv_Seurat)
  indiv_Seurat <- RunUMAP(indiv_Seurat, dims = 1:30)
  
  # Subset by barcode inflections
  indiv_Seurat <- CalculateBarcodeInflections(indiv_Seurat, barcode.column = "nCount_RNA", group.column = "sample", threshold.low = NULL, threshold.high = upper_cell_cutoff_rank)
  indiv_Seurat <- SubsetByBarcodeInflections(object = indiv_Seurat)
  
  # Perform DoubletFinder
  indiv_Seurat <- NormalizeData(indiv_Seurat, assay = "RNA", verbose = TRUE)
  indiv_QC <- RunPCA(indiv_Seurat, verbose = FALSE)
  
  pca_results <- Embeddings(indiv_Seurat, reduction = "pca")
  cat("PCA results dimension:", dim(pca_results), "\n")
  
  indiv_QC <- FindNeighbors(indiv_Seurat, dims = 1:30)
  indiv_QC <- FindClusters(indiv_QC, resolution = 0.6)
  annotations <- indiv_QC@meta.data$seurat_clusters
  
  sweep.res.list_organoid <- paramSweep(indiv_QC, sct = FALSE)
  sweep.stats_organoid <- summarizeSweep(sweep.res.list_organoid, GT = FALSE)
  bcmvn_organoid <- find.pK(sweep.stats_organoid)
  
  optimal_pK <- as.numeric(as.character(bcmvn_organoid$pK[which.max(bcmvn_organoid$BCmetric)]))
  homotypic.prop <- modelHomotypic(annotations)
  nExp_poi <- round(0.075 * nrow(indiv_QC@meta.data))
  nExp_poi.adj <- round(nExp_poi * (1 - homotypic.prop))
  
  cat("Expected doublets (nExp_poi.adj):", nExp_poi.adj, "\n")
  
  # Perform DoubletFinder
  indiv_QC <- doubletFinder(indiv_QC, pN = 0.25, PCs = 1:10, pK = optimal_pK, nExp = nExp_poi.adj, sct = FALSE)
  
  # Find the pANN column dynamically
  pANN_colname <- grep("^pANN_", colnames(indiv_QC@meta.data), value = TRUE)
  
  # Ensure that exactly one pANN column is found
  if (length(pANN_colname) != 1) {
    stop("Expected exactly one pANN column. Found: ", length(pANN_colname))
  }
  
  # Assign the pANN column as DoubletScore
  indiv_QC@meta.data$doubletScore <- indiv_QC@meta.data[, pANN_colname]
  
  # Find the DF.classifications column dynamically
  DF_classification_colname <- grep("^DF.classifications_", colnames(indiv_QC@meta.data), value = TRUE)
  
  # Ensure that exactly one DF.classifications column is found
  if (length(DF_classification_colname) != 1) {
    stop("Expected exactly one DF.classifications column. Found: ", length(DF_classification_colname))
  }
  
  # Subset to keep only singlets based on classifications
  singlet_cells <- indiv_QC@meta.data[, DF_classification_colname] == "Singlet"
  seurat_obj_clean <- subset(indiv_QC, cells = which(singlet_cells))
  
  
  
  
  # Calculate the number of cells after DoubletFinder and percentage metrics
  cells_after_DF <- ncol(seurat_obj_clean)
  percent_recovery_DF <- 100 * (cells_after_DF / cells_after_miQC)
  percent_removed_DF <- 100 - percent_recovery_DF
  cat("Sample:", sample_id, "- Cells before:", cells_before_QC, "- Cells after DF:", cells_after_DF, "\n")
  
  
  
  # Example checks after running DoubletFinder
  if (!exists("seurat_obj_clean") || is.null(seurat_obj_clean)) {
    cells_after_DF <- NA
    percent_recovery_DF <- NA
    percent_removed_DF <- NA
    warning(paste("DoubletFinder results are missing for sample:", sample_id))
  } else {
    cells_after_DF <- ncol(seurat_obj_clean)
    percent_recovery_DF <- 100 * (cells_after_DF / cells_after_QC)
    percent_removed_DF <- 100 - percent_recovery_DF
    print(paste('percent_recovery_DF:', percent_recovery_DF, ' percent_removed_DF: ', percent_removed_DF))
  }
  
  # Append results to dataframe only if all necessary data is present
  if (!is.na(cells_before_QC) && !is.na(cells_after_QC) && !is.na(cells_after_DF)) {
    new_row <- data.frame(
      sample = sample_id,
      genotype = genotype,
      cells_before_miQC = cells_before_miQC,
      cells_after_miQC = cells_after_miQC,
      percent_recovery_miQC = percent_recovery_miQC,
      percent_removed_miQC = percent_removed_miQC,
      cells_before_QC = cells_before_QC,
      cells_after_QC = cells_after_QC,
      percent_recovery_QC = percent_recovery_QC,
      percent_removed_QC = percent_removed_QC,
      cells_after_doubletF = cells_after_DF,
      percent_recovery_DF = percent_recovery_DF,
      percent_removed_DF = percent_removed_DF,
      stringsAsFactors = FALSE
    )
    cell_counts_DF <- rbind(cell_counts_DF, new_row)
  } else {
    warning(paste("Incomplete data for sample:", sample_id))
  }
  
  # Store the cleaned Seurat object
  organoid_list[[sample_id]] <- seurat_obj_clean
  rm(indiv_QC, indiv_Seurat, seurat_obj_clean)
}


write.csv(cell_counts_DF, file = paste0(dir, out, "cell_counts_DF.csv"), row.names = FALSE)

rm(S1, S2, S3, S4, S6, S8, S9, S10,
   S11, S12, S13, S14, S15, S16, S17, S18, S19, S20,
   S21, S22, S23, S24, S25, S26)

library(Seurat)

# Check if the organoid_list has more than one object
if (length(organoid_list) > 1) {
  # Rename the project name for all objects in the list
  for (i in 1:length(organoid_list)) {
    organoid_list[[i]]@project.name <- "APOE_Organoid"
  }
  
  # Ensure unique cell names across all objects
  for (i in 1:length(organoid_list)) {
    organoid_list[[i]] <- RenameCells(organoid_list[[i]], add.cell.id = paste0("Sample_", i, "_"))
  }
  
  # Initialize the combined Seurat object with the first object in the list
  organoid_combined <- organoid_list[[1]]
  
  # Merge all remaining Seurat objects in the list into the combined object
  for (i in 2:length(organoid_list)) {
    organoid_combined <- merge(organoid_combined, y = organoid_list[[i]], add.cell.ids = c(paste0("Sample_", 1, "_"), paste0("Sample_", i, "_")))
  }
  
  print('Ran organoid_combined <- merge(...) successfully.')
  
} else if (length(organoid_list) == 1) {
  # Only one object, rename the project and set it as the combined object
  organoid_combined <- organoid_list[[1]]
  organoid_combined@project.name <- "APOE_Organoid"
  print('Only one object in the list, renamed the project successfully.')
  
} else {
  # Error catch for no objects in the list
  stop("No data processed successfully.")
}

# Set the default assay to 'RNA' after merging
DefaultAssay(organoid_combined) <- "RNA"
print('Ran DefaultAssay(...) successfully.')

print('Saving organoid_combined as organoid2.rds')
saveRDS(organoid_combined, file = paste0(dir, out, 'organoid2.rds'))


## Save the cell_counts_DF object

write.csv(cell_counts_DF, file = paste0(dir, out, "cell_counts_DF.csv"), row.names = FALSE)

# Aggregate the dataframe by genotype
aggregated_data <- cell_counts_DF %>%
  group_by(genotype) %>%
  summarise(
    total_cells_before_miQC = sum(cells_before_miQC, na.rm = TRUE),
    total_cells_after_miQC = sum(cells_after_miQC, na.rm = TRUE),
    total_cells_before_QC = sum(cells_before_QC, na.rm = TRUE),
    total_cells_after_QC = sum(cells_after_QC, na.rm = TRUE),
    total_cells_after_doubletF = sum(cells_after_doubletF, na.rm = TRUE),
    .groups = 'drop'
  ) %>%
  mutate(
    percent_recovery_miQC = 100 * (total_cells_after_miQC / total_cells_before_miQC),
    percent_recovery_QC = 100 * (total_cells_after_QC / total_cells_after_miQC),
    percent_recovery_DF = 100 * (total_cells_after_doubletF / total_cells_after_QC),
    percent_removed_DF = 100 * (1 - (total_cells_after_doubletF / total_cells_after_QC))
  )


write.csv(aggregated_data, file = paste0(dir, out, "aggregated_cell_counts_DF.csv"), row.names = FALSE)





##### Step 3. Visualize the results of QC above            #####################
library(stringr)

## Now, we go through the miQC again, but this time with a 0.97 posterior cutoff
## We are running through this a second time because we were originally just doing it 
# 
# 
# p <- VlnPlot(organoid_combined, pt.size = 0, group.by = 'sample', features = 'percent.mt', log = TRUE)
# ggsave(filename = "outputs/images/02-Quality_Control/combined_Vln_plots_before_QC.png", 
#        plot = p, width = 6, height = 4)
# 
# 
# 
# # Step 3. Visualize after-QC metrics                                             ####
# 
# features <- c("nFeature_RNA", "nCount_RNA", "percent.mt")
# ## plot the before-QC metrics
# vps<-VlnPlot(object = organoid_combined, features = features, ncol = 3, layer = "counts", pt.size = 0)
# vps
# vp1<-vps[[1]] +geom_hline(yintercept = 200,colour="black", na.rm = TRUE) +labs(x = "organoid")
# vp2<-vps[[2]]+geom_hline(yintercept = upper_cell_cutoff,colour="black") +geom_hline(yintercept = lower_cell_cutoff,colour="black") +labs(x = "organoid") #+ scale_y_continuous(limits=c(0,100000))#+ scale_y_continuous(trans='log10')
# vp3<-vps[[3]]+geom_hline(yintercept = 5, colour="black") +labs(x = "organoid")
# p <- vp1 +  theme(axis.text.x = element_blank()) +vp2 + theme(axis.text.x = element_blank()) +vp3 + theme(axis.text.x = element_blank())
# print(p)
# 
# ggsave(filename = "outputs/images/02-Quality_Control/After_QC_Vln_Plot.png", 
#        plot = p, width = 6, height = 4)

