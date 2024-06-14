## 02-Quality-Control.R
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
# .libPaths()
# source("~/.Rprofile")

# Step 1: Read in the indiv samples' Seurat object files & combine into one ####

###* denotes changes made by AKG
###* Then, combine them into one Seurat object. This object will be used to
###* be a point of comparison for before/after quality control steps
setwd("/projectnb/tcwlab/LabMember/akg/projects/APOE_Organoid")

dir <- "/projectnb/tcwlab/LabMember/akg/projects/APOE_Organoid/"
out<-'outputs/02-Quality-Control/'


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


###* Proceed by merging all the samples into a single object, called organoid0
APOE22_samples_group <- c("S3", "S9", "S12", "S15", "S21", "S25")
APOE33_samples_group <- c("S1", "S8", "S14", "S16", "S19", "S23")
APOE33Ch_samples_group <- c("S2", "S6", "S13", "S18", "S20", "S24")
APOE44_samples_group <- c("S4", "S10", "S11", "S17", "S22", "S26")

sample_names_group22 <- c("sample3", "sample9", "sample12", "sample15", "sample21", "sample25")
sample_names_group33 <- c("sample1", "sample8", "sample14", "sample16", "sample19", "sample23")
sample_names_group33Ch <- c("sample2", "sample6", "sample13", "sample18", "sample20", "sample24")
sample_names_group44 <- c("sample4", "sample10", "sample11", "sample17", "sample22", "sample26")

# Loop through each sample and add a metadata column for sample labels
# for APOE22
for (i in seq_along(APOE22_samples_group)) {
  object_name <- APOE22_samples_group[i]
  sample_name <- sample_names_group22[i]
  
  eval(parse(text = paste0(object_name, "[['genotype']] <- 'APOE22'")))
}

for (i in seq_along(APOE33_samples_group)) {
  object_name <- APOE33_samples_group[i]
  sample_name <- sample_names_group33[i]
  
  eval(parse(text = paste0(object_name, "[['genotype']] <- 'APOE33'")))
}

for (i in seq_along(APOE33Ch_samples_group)) {
  object_name <- APOE33Ch_samples_group[i]
  sample_name <- sample_names_group33Ch[i]
  
  eval(parse(text = paste0(object_name, "[['genotype']] <- 'APOE33Ch'")))
}

for (i in seq_along(APOE44_samples_group)) {
  object_name <- APOE44_samples_group[i]
  sample_name <- sample_names_group44[i]
  
  eval(parse(text = paste0(object_name, "[['genotype']] <- 'APOE44'")))
}

################################################################################




# Initialize a data frame to store cell counts before and after QC
cell_counts_DF <- data.frame(
  sample = character(), 
  genotype = character(),
  cells_before_QC = integer(), 
  cells_after_QC = integer(),
  percent_recovery_QC = numeric(),
  percent_removed_QC = numeric(),
  cells_after_doubletF = integer(),
  percent_recovery_DF = numeric(),
  percent_removed_DF = numeric(),
  stringsAsFactors = FALSE  # Ensure strings are not converted to factors
)

# init empty list to store sample's individual RDS files
organoid_list <- list()

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

samples_list <- list(S1, S2, S3, S4, S6, S8, S9, S10, S11, S12, S13, S14,
             S15, S16, S17, S18, S19, S20, S21, S22, S23, S24, S25, S26)


# Loop through each sample
for (i in 1:length(samples_list)) {
  cat("Processing sample", i, "\n")
  indiv_Seurat <- samples_list[[i]]
  
  sample_id <- unique(indiv_Seurat@meta.data$sample)
  genotype <- unique(indiv_Seurat@meta.data$genotype)
  
  cells_before_QC <- ncol(indiv_Seurat)
  cat("Sample ID:", sample_id, "- Cells before QC:", cells_before_QC, "\n")
  
  # Subset mt_gene_names to include only genes present in the Seurat object
  counts_data <- GetAssayData(indiv_Seurat, slot = "counts", assay = "RNA")
  valid_mt_gene_names <- mt_gene_names[mt_gene_names %in% rownames(counts_data)]
  additional_mt_genes <- rownames(counts_data)[grepl("^mt-", rownames(counts_data), ignore.case = TRUE)]
  valid_mt_gene_names <- unique(c(valid_mt_gene_names, additional_mt_genes))
  indiv_Seurat[["percent.mt"]] <- PercentageFeatureSet(indiv_Seurat, features = valid_mt_gene_names)
  
  # Additionally label the ribosomal percentage
  indiv_Seurat[["percent.ribo"]] <- PercentageFeatureSet(indiv_Seurat, pattern = "^RP[SL]")
  
  lower_cell_cutoff <- lower_nCount_RNA_cutoffs[i]
  cat("Lower cell cutoff:", lower_cell_cutoff, "\n")
  sorted_counts <- sort(indiv_Seurat@meta.data$nCount_RNA, decreasing = FALSE)
  upper_cell_cutoff_rank <- ceiling(0.99 * length(sorted_counts))
  cat("Upper cell cutoff rank:", upper_cell_cutoff_rank, "\n")
  
  # Subset the Seurat object based on nCount_RNA
  indiv_Seurat <- subset(indiv_Seurat, subset = nCount_RNA > lower_cell_cutoff)
  
  # Perform SCTransform
  indiv_Seurat <- SCTransform(indiv_Seurat, vars.to.regress = "percent.mt", verbose = TRUE)
  
  # Calculate and subset by barcode inflections
  indiv_Seurat <- CalculateBarcodeInflections(
    indiv_Seurat,
    barcode.column = "nCount_RNA",
    group.column = "sample",
    threshold.low = NULL,
    threshold.high = upper_cell_cutoff_rank
  )
  indiv_Seurat <- SubsetByBarcodeInflections(object = indiv_Seurat)
  
  # Filter genes expressed in less than 3 cells for both SCT and RNA assays
  for (assay in c("SCT", "RNA")) {
    raw_counts_mat <- GetAssayData(indiv_Seurat, slot = "counts", assay = assay)
    genes_to_keep <- rowSums(raw_counts_mat > 0) >= 3
    indiv_Seurat <- subset(indiv_Seurat, features = names(genes_to_keep[genes_to_keep]))
  }
  
  # Subset for nFeature_RNA
  nFeature_RNA_values <- indiv_Seurat@meta.data$nFeature_RNA
  sorted_nFeature_RNA <- sort(nFeature_RNA_values, decreasing = FALSE)
  upper_index_nFeature <- ceiling(0.99 * length(sorted_nFeature_RNA))
  upper_feature_cutoff <- sorted_nFeature_RNA[upper_index_nFeature]
  cat("Upper feature cutoff:", upper_feature_cutoff, "\n")
  
  # Set thresholds for nFeature_RNA
  lower_limit_nFeature <- 200
  indiv_Seurat <- subset(indiv_Seurat, subset = nFeature_RNA > lower_limit_nFeature)
  
  cells_after_QC <- ncol(indiv_Seurat)
  percent_recovery_QC <- 100 * (cells_after_QC / cells_before_QC)
  
  # Perform DoubletFinder
  indiv_Seurat <- NormalizeData(indiv_Seurat, assay = "RNA", verbose = TRUE)
  indiv_QC <- RunPCA(indiv_Seurat, verbose = FALSE)
  
  # Verify PCA results
  pca_results <- Embeddings(indiv_QC, reduction = "pca")
  cat("PCA results dimension:", dim(pca_results), "\n")
  
  indiv_QC <- RunUMAP(indiv_QC, dims = 1:10, reduction = "pca")
  indiv_QC <- FindNeighbors(indiv_QC, dims = 1:15, reduction = "pca")
  indiv_QC <- FindClusters(indiv_QC, resolution = 0.6)
  annotations <- indiv_QC@meta.data$seurat_clusters
  
  sweep.res.list_organoid <- paramSweep(indiv_QC, PCs = 1:10, sct = TRUE)
  sweep.stats_organoid <- summarizeSweep(sweep.res.list_organoid, GT = FALSE)
  bcmvn_organoid <- find.pK(sweep.stats_organoid)
  
  optimal_pK <- bcmvn_organoid$pK[which.max(bcmvn_organoid$BCmetric)]
  homotypic.prop <- modelHomotypic(annotations)
  nExp_poi <- round(0.075 * nrow(indiv_QC@meta.data))
  nExp_poi.adj <- round(nExp_poi * (1 - homotypic.prop))
  
  # Check the number of expected doublets
  cat("Expected doublets (nExp_poi.adj):", nExp_poi.adj, "\n")
  
  indiv_QC <- doubletFinder_v3(indiv_QC, PCs = 1:10, pN = 0.25, pK = optimal_pK, nExp = nExp_poi.adj, sct = TRUE)
  
  cells_after_DF <- ncol(indiv_QC)
  percent_recovery_DF <- 100 * (cells_after_DF / cells_after_QC)
  percent_removed_DF <- 100 - percent_recovery_DF
  cat("Sample:", sample_id, "- Cells before:", cells_before_QC, "- Cells after DF:", cells_after_DF, "\n")
  
  # Append results to dataframe
  new_row <- data.frame(
    sample = sample_id,
    genotype = genotype,
    cells_before_QC = cells_before_QC,
    cells_after_QC = cells_after_QC,
    percent_recovery_QC = percent_recovery_QC,
    cells_after_doubletF = cells_after_DF,
    percent_recovery_DF = percent_recovery_DF,
    percent_removed_DF = percent_removed_DF,
    stringsAsFactors = FALSE
  )
  
  cell_counts_df <- rbind(cell_counts_df, new_row)
}


if (length(organoid_list) > 1) {
  # Merge all Seurat objects in the organoid_list
  organoid1 <- Reduce(function(x, y) merge(x, y, project = "APOE_Organoid"), organoid_list)
} else if (length(organoid_list) == 1) {
  # Only one object, rename the project
  organoid1 <- organoid_list[[1]]
  organoid1@project.name <- "APOE_Organoid"
} else {
  # Error catch
  message("No data processed successfully.")
}

saveRDS(organoid1, file = paste0(dir, out, 'organoid1.rds'))


## Save the cell_counts_DF object

write.csv(cell_counts_df, file = paste0(dir, out, "cell_counts_DF.csv"), row.names = FALSE)





##### Step 3. Visualize the results of QC above            #####################
library(stringr)

## Now, we go through the miQC again, but this time with a 0.97 posterior cutoff
## We are running through this a second time because we were originally just doing it 


VlnPlot(organoid1, pt.size = 0, group.by = 'sample', features = 'percent.mt', log = TRUE)
ggsave(filename = "outputs/images/02-Quality-Control/combined_Vln_plots_before_QC.png", 
       plot = p, width = 6, height = 4)



# Step 3. Visualize after-QC metrics                                             ####

features <- c("nFeature_RNA", "nCount_RNA", "percent.mt")
## plot the before-QC metrics
vps<-VlnPlot(object = organoid1, features = features, ncol = 3, layer = "counts", pt.size = 0)
vps
vp1<-vps[[1]] +geom_hline(yintercept = 200,colour="black", na.rm = TRUE) +labs(x = "organoid")
vp2<-vps[[2]]+geom_hline(yintercept = upper_cell_cutoff,colour="black") +geom_hline(yintercept = lower_cell_cutoff,colour="black") +labs(x = "organoid") #+ scale_y_continuous(limits=c(0,100000))#+ scale_y_continuous(trans='log10')
vp3<-vps[[3]]+geom_hline(yintercept = 5, colour="black") +labs(x = "organoid")
p <- vp1 +  theme(axis.text.x = element_blank()) +vp2 + theme(axis.text.x = element_blank()) +vp3 + theme(axis.text.x = element_blank())


output_directory <- paste0(dir, 'outputs/images/02-Quality-Control/After_QC_Vln_Plot.png')

# Save the plot
ggsave(output_file, plot = p, width = 10, height = 8)
