#Final QC

### Note: This is taken from Deepti Murthy and adapted by Andrew Gjelsteen,
### most of the comments are hers ###
### Andrew's comments are denoted by ##*  *##
#TCW Single-Cell RNA Seq Data
#Dataset Source: TCW_14311_run1

#Date: 11/10/22 - 11/13/22
#Date: 1/1/23 (Altering Cluster Resolution)

#Note, prior to this analysis:
#The seurat_object.rds" object was created using RNA and HTO information, including
#only those cells with at least 10 nUMI for any given HTO.
#HTODemux was run using default parameters.
#Across-sample Doublets and dropouts have been removed from the dataset.

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
.libPaths()
source("~/.Rprofile")

###* denotes changes made by AKG
###* Step 1: Read in the individual samples' Seurat object files ####
###* Then, combine them into one Seurat object. This object will be used to
###* be a point of comparison for before/after quality control steps
setwd("/projectnb/tcwlab/LabMember/akg/projects/APOE_Organoid")

dir <- "/projectnb/tcwlab/LabMember/akg/projects/APOE_Organoid"
out<-'outputs/02-Final_QC'


S1 <- readRDS(file = "./outputs/simpleafSeurat/sample1.rds")
S2 <- readRDS(file = "./outputs/simpleafSeurat/sample2.rds")
S3 <- readRDS(file = "./outputs/simpleafSeurat/sample3.rds")
S4 <- readRDS(file = "./outputs/simpleafSeurat/sample4.rds")
###* No Sample 5 (deemed poor quality by researchers)
S6 <- readRDS(file = "./outputs/simpleafSeurat/sample6.rds")
###* No Sample 7 (deemed poor quality by researchers)
S8 <- readRDS(file = "./outputs/simpleafSeurat/sample8.rds")
S9 <- readRDS(file = "./outputs/simpleafSeurat/sample9.rds")
S10 <- readRDS(file = "./outputs/simpleafSeurat/sample10.rds")
S11 <- readRDS(file = "./outputs/simpleafSeurat/sample11.rds")
S12 <- readRDS(file = "./outputs/simpleafSeurat/sample12.rds")
S13 <- readRDS(file = "./outputs/simpleafSeurat/sample13.rds")
S14 <- readRDS(file = "./outputs/simpleafSeurat/sample14.rds")
S15 <- readRDS(file = "./outputs/simpleafSeurat/sample15.rds")
S16 <- readRDS(file = "./outputs/simpleafSeurat/sample16.rds")
S17 <- readRDS(file = "./outputs/simpleafSeurat/sample17.rds")
S18 <- readRDS(file = "./outputs/simpleafSeurat/sample18.rds")
S19 <- readRDS(file = "./outputs/simpleafSeurat/sample19.rds")
S20 <- readRDS(file = "./outputs/simpleafSeurat/sample20.rds")
S21 <- readRDS(file = "./outputs/simpleafSeurat/sample21.rds")
S22 <- readRDS(file = "./outputs/simpleafSeurat/sample22.rds")
S23 <- readRDS(file = "./outputs/simpleafSeurat/sample23.rds")
S24 <- readRDS(file = "./outputs/simpleafSeurat/sample24.rds")
S25 <- readRDS(file = "./outputs/simpleafSeurat/sample25.rds")
S26 <- readRDS(file = "./outputs/simpleafSeurat/sample26.rds")


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
  
  eval(parse(text = paste0(object_name, "[['sample']] <- '", sample_name, "'")))
  eval(parse(text = paste0(object_name, "[['genotype']] <- 'APOE22'")))
}

for (i in seq_along(APOE33_samples_group)) {
  object_name <- APOE33_samples_group[i]
  sample_name <- sample_names_group33[i]
  
  eval(parse(text = paste0(object_name, "[['sample']] <- '", sample_name, "'")))
  eval(parse(text = paste0(object_name, "[['genotype']] <- 'APOE33'")))
}

for (i in seq_along(APOE33Ch_samples_group)) {
  object_name <- APOE33Ch_samples_group[i]
  sample_name <- sample_names_group33Ch[i]
  
  eval(parse(text = paste0(object_name, "[['sample']] <- '", sample_name, "'")))
  eval(parse(text = paste0(object_name, "[['genotype']] <- 'APOE33Ch'")))
}

for (i in seq_along(APOE44_samples_group)) {
  object_name <- APOE44_samples_group[i]
  sample_name <- sample_names_group44[i]
  
  eval(parse(text = paste0(object_name, "[['sample']] <- '", sample_name, "'")))
  eval(parse(text = paste0(object_name, "[['genotype']] <- 'APOE44'")))
}

################################################################################
#########@######################################################################


##* Merge all Seurat objects into one *##
##* Can easily process it as one object *##
##* 
organoid0 <- merge(S1, y = c(S2, S3, S4, S6, S8, S9, S10,
                             S11, S12, S13, S14, S15, S16, S17, S18, S19,
                             S20, S21, S22, S23, S24, S25, S26), 
                   add.cell.ids = c("S1", "S2", "S3", "S4", "S6", "S8", "S9", "S10",
                                    "S11", "S12", "S13", "S14", "S15", "S16", "S17", "S18", "S19",
                                    "S20", "S21", "S22", "S23", "S24", "S25", "S26"),
                   project = "APOE_Organoid"
)


saveRDS(organoid0, file = paste0(dir, '/outputs/organoid0.rds'))
unique(organoid0@meta.data$sample)
organoid0 <- readRDS('./outputs/organoid0.rds')
# rm(S1, S2, S3, S4, S6, S8, S9, S10,
#    S11, S12, S13, S14, S15, S16, S17, S18, S19,
#    S20, S21, S22, S23, S24, S25, S26)

rownames(organoid0)

# Initialize a data frame to store cell counts before and after QC
cell_counts <- data.frame(sample = character(), 
                          cells_before_QC = integer(), 
                          cells_after_QC = integer())

# Get unique sample identifiers
samples <- unique(organoid0@meta.data$sample)

# Define mitochondrial genes (adjust according to your organism)
mt_gene_names <- grep("^MT-", rownames(organoid0), value = TRUE)

organoid_list <- list()

# Loop over each sample
for (sample_id in samples) {
  # Subset Seurat object for the current sample
  T_QC <- subset(organoid0, subset = sample == sample_id)
  
  # Record the number of cells before QC
  cells_before <- ncol(T_QC)
  
  # Filter genes expressed in less than 3 cells
  raw_counts_mat <- T_QC@assays$RNA@counts
  genes_to_keep <- Matrix::colSums(raw_counts_mat != 0) >= 3
  T_QC <- subset(T_QC, features = rownames(T_QC)[genes_to_keep])
  
  # Calculate upper limit, set at 99.5 percentile
  upper_limit <- quantile(T_QC@meta.data$nCount_RNA, probs = 0.995)
  
  T_QC <- CalculateBarcodeInflections(
    T_QC,
    barcode.column = "nCount_RNA",
    group.column = "orig.ident",
    threshold.low = 1000,
    threshold.high = NULL
  )
  T_QC <- SubsetByBarcodeInflections(object = T_QC)
  
  # Record the number of cells after QC
  cells_after <- ncol(T_QC)
  
  # append the QC'd object onto the organoid_list (so we can combine after)
  organoid_list[[sample_id]] <- T_QC
  
  # Append the cell count information to the data frame
  cell_counts <- rbind(cell_counts, data.frame(sample = sample_id,
                                               cells_before_QC = cells_before,
                                               cells_after_QC = cells_after))
}

# merge all QC'd Seurat objects into one (all alrdy have unique cell identifiers)
if(length(organoid_list) > 1) {
  # do.call to apply the merge function to a list of Seurat objects
  organoid1 <- do.call(merge, c(list(x = organoid_list[[1]], y = organoid_list[-1]), list(project = "APOE_Organoid")))
} else if(length(organoid_list) == 1) {
  # Only one object, simply rename the project if necessary
  organoid1 <- RenameProject(organoid_list[[1]], "APOE_Organoid")
} else {
  # No objects were processed successfully
  message("No data processed successfully.")
}



organoid1 <- merge(x = organoid_list[[1]], y = organoid_list[-1], 
                           add.cell.ids = names(organoid_list[-1]), 
                           project = "APOE_Organoid")



# View the cell counts data frame
print(cell_counts)

# Step 2. Label and quantify Mitochrondrial DNA and Ribosomal DNA ####

### Note: A lot of the mitochondrial genes, when undergoing the mapping of 
### ENSEMBL_IDs to gene symbols, somehow got excluded from this process.
### In order to deal with this, we will define a list of all the ENSEMBL IDs AND
### their gene symbols for
### all 37 Human Mitochondrial genes

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
  "TC", "TY", "TA", "TQ", "RNR1", "RNR2"
)




# Step 3. Apply SoupX Correction (Ambient RNA Removal) ####


# Running 
library(tools)
library(Seurat)
library(data.table)
library(ggplot2)
library(SingleCellExperiment)
library(fishpond)
library(scater)
library(org.Hs.eg.db)
library(Matrix)
library(biomaRt)
library(SoupX)

## initialize a dataframe which will contain the before/after SoupX correction
## cell counts
SoupXSummary <- data.frame(
  sampleID = character(24),
  before_SoupX = numeric(24),
  after_SoupX = numeric(24),
  percent_removed = numeric(24),
  stringsAsFactors = FALSE
)


# 'organoid0' is combined Seurat object
sample_ids <- unique(organoid0@meta.data$sample)
# sample_ids
seurat_objects <- SplitObject(organoid0, split.by = "sample")


source <- '/projectnb/tcwlab/LabMember/akg/APOE_Jorganoid_project/outputs/'

### Function to convert ENSEMBL_ID names to gene symbols
getGeneSymbols <- function(ensembl_ids) {
  require(org.Hs.eg.db)
  gene_symbols <- mapIds(org.Hs.eg.db,
                         keys = ensembl_ids,
                         column = "SYMBOL",
                         keytype = "ENSEMBL",
                         multiVals = "first")
  return(gene_symbols)
}

# Load raw counts
custom_format <- list("counts" = c("U","S","A"))

# Function to apply SoupX correction
apply_soupX_correction <- function(seurat_obj, raw_counts_directory) {
  
  # # set a "raw counts" count for each of the individual Seurat samples' objects
  # seurat_obj[["original_counts"]] <- CreateAssayObject(counts = seurat_obj@assays$RNA@counts) #raw counts
  # 
  sce <- fishpond::loadFry(raw_counts_directory,
                           outputFormat = custom_format)
  # replace Ensembl IDs in row names w gene symbols
  ensembl_ids_sce <- rownames(sce)
  gene_symbols_sce <- getGeneSymbols(ensembl_ids_sce)
  
  # If no gene symbol is found, keep the original Ensembl ID
  updated_rownames <- ifelse(is.na(gene_symbols_sce[ensembl_ids_sce]), ensembl_ids_sce, gene_symbols_sce[ensembl_ids_sce])
  rownames(sce) <- updated_rownames
  
  #create the seurat object by filtering for selected cells
  S1Counts <- counts(sce)
  
  raw <- S1Counts
  #SoupX Correction
  #modify genes in raw
  genes <- intersect(rownames(seurat_obj), rownames(raw))
  length(genes)
  raw <- raw[genes,]

  #run SoupX algo (this package is what I ended up using instead of DecontX)
  sc <- SoupChannel(raw, seurat_obj@assays$RNA@counts)
  sc <- setClusters(sc, seurat_obj$initial_seurat_clustering_0.5)
  sc <- autoEstCont(sc, doPlot = FALSE)
  soup_out <- adjustCounts(sc, roundToInt = TRUE) #adjusted counts matrix
  
  seurat_obj[["SoupX_counts"]] <- CreateAssayObject(counts = soup_out)
  seurat_obj@assays$RNA@counts <- soup_out #USED
  
  #calculate percent of mRNA removed
  percent_removed_soup <- (1 - (sum(seurat_obj@assays$RNA@counts) / sum(seurat_obj@assays$original_counts@counts)))*100
  
  
  return(seurat_obj)
}

# Directories for raw counts of each sample
raw_counts_directories <- c("1_S16/af_quant/", "2_S8/af_quant/", "3_S9/af_quant/", "4_S14/af_quant/", "6_S23/af_quant/",
                            "8_S17/af_quant/", "9_S21/af_quant/", "10_S7/af_quant/", "11_S1/af_quant/", "12_S20/af_quant/", "13_S11/af_quant/",
                            "14_S13/af_quant/", "15_S5/af_quant/", "16_S24/af_quant/", "17_S19/af_quant/", "18_S22/af_quant/",
                            "19_S3/af_quant/", "20_S15/af_quant/", "21_S6/af_quant/", "22_S2/af_quant/", "23_S12/af_quant/", "24_S10/af_quant/",
                            "25_S18/af_quant/", "26_S4/af_quant/") # paths to quant files

organoid1 <- readRDS('outputs/simpleafSeurat/organoid1.rds')

# Process each Seurat object
for (i in seq_along(seurat_objects)) {
  sampleID <- names(seurat_objects)[i]
  print(sampleID)
  raw_counts_dir <- paste(source, raw_counts_directories[i], sep = '') 
  print(raw_counts_dir)
  original_cell_count <- ncol(seurat_objects[[sampleID]]@assays$RNA@counts)
  
  seurat_objects[[sampleID]] <- apply_soupX_correction(seurat_objects[[sampleID]], raw_counts_dir)
  
  after_cell_count <- ncol(seurat_objects[[sampleID]]@assays$RNA@counts)
  percent_removed <- (1 - (after_cell_count / original_cell_count)) * 100
  
  samples_summary[i, ] <- c(sampleID, original_cell_count, after_cell_count, percent_removed)
}


raw_counts_dir <- paste(source, "1_S16/af_quant/", sep = '') 
# # set a "raw counts" count for each of the individual Seurat samples' objects
# seurat_obj[["original_counts"]] <- CreateAssayObject(counts = seurat_obj@assays$RNA@counts) #raw counts
# 
sce <- fishpond::loadFry(raw_counts_dir,
                         outputFormat = custom_format)
# Replace Ensembl IDs in row names with gene symbols
ensembl_ids_sce <- rownames(sce)
gene_symbols_sce <- getGeneSymbols(ensembl_ids_sce)

# if no gene symbol is found, keep the original Ensembl ID
updated_rownames <- ifelse(is.na(gene_symbols_sce[ensembl_ids_sce]), 
                           ensembl_ids_sce, gene_symbols_sce[ensembl_ids_sce])
rownames(sce) <- updated_rownames

#create the seurat object by filtering for selected cells
S1Counts <- counts(sce)


seurat_obj <- seurat_objects[["sample1"]]
colnames(seurat_obj)
colnames(S1Counts)

colnames(S1Counts) <- paste0("S1_", colnames(S1Counts))

raw <- S1Counts
#SoupX Correction
#modify genes in raw
genes <- intersect(rownames(seurat_obj), rownames(raw))
length(genes)
raw <- raw[genes,]

#normalize data
seurat_obj <- NormalizeData(seurat_obj)
seurat_obj <- FindVariableFeatures



#scale all
seurat_obj <- ScaleData(seurat_obj, assay = "RNA", verbose = FALSE)


seurat_obj <- RunPCA(seurat_obj, verbose = FALSE)
seurat_obj <- RunUMAP(seurat_obj, dims = 1:50, verbose = FALSE)
seurat_obj <- FindNeighbors(seurat_obj, dims = 1:50, verbose = FALSE)
seurat_obj <- FindClusters(seurat_obj, resolution = 0.5, verbose = FALSE)



#run SoupX algo (this package is what I ended up using instead of DecontX)
sc <- SoupChannel(raw, seurat_obj@assays$RNA@counts[genes,])
seurat_obj@meta.data$

sc <- setClusters(sc, seurat_obj$initial_seurat_clustering_0.5) #### Ask Deepti (mention that this clustering is called for before the cluster is actually performed)
sc <- autoEstCont(sc, doPlot = FALSE)
soup_out <- adjustCounts(sc, roundToInt = TRUE) #adjusted counts matrix

seurat_obj[["SoupX_counts"]] <- CreateAssayObject(counts = soup_out)
seurat_obj@assays$RNA@counts <- soup_out #USED

#calculate percent of mRNA removed
percent_removed_soup <- (1 - (sum(seurat_obj@assays$RNA@counts) / sum(seurat_obj@assays$original_counts@counts)))*100


return(seurat_obj)

##### Deepti's code to set original idents on objects #####
##### to keep track of how many mt / ribo genes were removed #####

gene_names <- rownames(organoid0@assays$RNA@counts)
gene_names
###* Filter gene names based on mitochondrial ENSEMBL IDs/gene symbols list *###
mt_genes <- gene_names[gene_names %in% mt_gene_names]

# Print the mitochondrial genes
print(mt_genes)
###* Assign percent.mt and percent.ribo to the organoid0 object:            *###
###* 
# Set Idents, calculate MT and Ribosomal percentages
Idents(organoid0) <- organoid0$orig.ident

# Subset mt_gene_names to include only genes present in the Seurat object
valid_mt_gene_names <- mt_gene_names[mt_gene_names %in% rownames(organoid0[["RNA"]]@counts)]

# Recalculate the percentage of mitochondrial genes
organoid0[["percent.mt"]] <- PercentageFeatureSet(organoid0, features = valid_mt_gene_names)
organoid0[["percent.ribo"]] <- PercentageFeatureSet(organoid0, pattern = "^RP[SL]")

###* Let's now plot the violin plots for percent.mt, percent.ribo, etc.     *###


# organoid0 is the Seurat object and 'genotype' is a column in meta.data
unique_genotypes <- unique(organoid0@meta.data$genotype)

# Output directory
out <- "/projectnb/tcwlab/LabMember/akg/APOE_Jorganoid_project/outputs/images/" # dir path

# Loop through each genotype
library(Seurat)
library(ggplot2)
library(reshape2)

rm(S1, S2, S3, S4, S6, S8, S9, S10, 
   S11, S12, S13, S14, S15, S16, S17, S18, S19, S20, 
   S21, S22, S23, S24, S25, S26) #clear recent memory

# Define the features to plot
features <- c("nFeature_RNA", "nCount_RNA", "percent.mt")

# Loop through each genotype and make violin plots for each APOE genotype
for (genotype in unique_genotypes) {
  # Subset Seurat object by genotype
  organoid_subset <- subset(organoid0, subset = genotype == genotype)
  
  # Create violin plots
  vps <- VlnPlot(object = organoid_subset, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3, slot = "counts", pt.size = 0)
  
  # Modify the plots to spread them out evenly
  vps <- vps + plot_layout(guides = 'collect') & theme(legend.position = 'none')
  
  # Add the genotype as the overall figure title
  combined_plot <- vps + plot_annotation(title = genotype, theme = theme(plot.title = element_text(hjust = 0.5)))
  
  # Construct the full path for the output file
  output_file_path <- file.path(out, paste0("before_QC_violin_plots_", genotype, ".png"))
  
  # Save the combined plot to the specified directory
  ggsave(output_file_path, plot = combined_plot, width = 12, height = 6, units = "in")
}

out <- "/projectnb/tcwlab/LabMember/akg/APOE_Jorganoid_project/outputs/miQC"
filename = 'organoid0.rds'
###* Saving our orgnaoid0 object. Good idea to do so, so we can restart R   *###
###* prior to next step and load in the object                              *###
saveRDS(organoid0, file.path(out, filename))

################################################################################
#### Following Deepti's Pipeline here: #########################################
#### Step 2: Remove Low Quality Cells  #########################################
################################################################################

organoid0 <- readRDS('outputs/simpleafSeurat/organoid0.rds')

organoid0@active.assay
#set assay to RNA
organoid0@active.assay <- "RNA" #8,149 cells & 33,546 genes

### Before deciding cutoffs, use vln plots to visualize the number of detected
### genes expressed, number of UMIs detected, and pt.mitochondrial genes:

vps<-VlnPlot(object = organoid0, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3, slot = "counts", pt.size = 0)
vp1<-vps[[1]] +geom_hline(yintercept = 200,colour="red") +labs(x = "organoid0")
vp2<-vps[[2]]#+geom_hline(yintercept = 75000,colour="red") +geom_hline(yintercept = 500,colour="red") +labs(x = "organoid0")
vp3<-vps[[3]]#+geom_hline(yintercept = 14,colour="red") +labs(x = "organoid0")

vps <- VlnPlot(object = organoid0, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3, slot = "counts", pt.size = 0)
vp1 <- vps[[1]] + labs(x = "organoid0")
vp2 <- vps[[2]] + labs(x = "organoid0")
vp3 <- vps[[3]] + labs(x = "organoid0")




pdf(paste(dir, "/QC_metric_VP_before_cutoffs.pdf", sep = ""), width = 10, height = 10)
p <- vp1 +  theme(axis.text.x = element_blank()) +vp2 + theme(axis.text.x = element_blank()) +vp3 + theme(axis.text.x = element_blank()) 
plot(p)
dev.off()

### We can also visualize a scatter plot of nFeature_RNA vs. nCount_RNA
### which ought to have a linear relationship:

#FeatureScatter is used to visualize feature-feature relationships
#Ensure linear relationship between nFeature_RNA and nCount_RNA
plot1 <- FeatureScatter(organoid0, feature1 = "nCount_RNA", feature2 = "percent.mt") + labs(legend = "TCW_14311_run1") & NoLegend() 
plot2 <- FeatureScatter(organoid0, feature1 = "nCount_RNA", feature2 = "nFeature_RNA") +labs(legend = "TCW_14311_run1") & NoLegend() 
pdf(paste(dir, "/images/QC_metrics_FP_before_cutoffs.pdf", sep = ""), width = 10, height = 5)
plot1 + plot2 
dev.off()

# Expect this to not be linear prior to following QC steps

#### Remove low quality cells based on miQC threshold ##########################

#BiocManager::install("flexmix")
organoid1 <- readRDS('outputs/simpleafSeurat/organoid0.rds')
library(miQC)
library(Seurat)
library(flexmix)
library(SeuratWrappers)
setwd("/projectnb/tcwlab/LabMember/akg/APOE_Jorganoid_project/")

dir <- "/projectnb/tcwlab/LabMember/akg/APOE_Jorganoid_project/outputs"
out <- 'outputs/miQC'

#Run miQC
#note: posterior cutoff = the posterior probability of a cell being part of the compromised distribution, a number between 0 and 1.
#Any cells below the appointed cutoff will be marked to keep. Defaults to 0.75.

organoid1 <- RunMiQC(organoid1, percent.mt = "percent.mt", nFeature_RNA = "nFeature_RNA", posterior.cutoff = 0.97, model.slot = "flexmix_model") 

pdf(paste(dir, "/miQC/stringent_miQC_probability_plot.pdf", sep = ""), width = 15, height = 10)
PlotMiQC(organoid1, color.by = "miQC.probability") + ggplot2::scale_color_gradient(low = "grey", high = "purple")
dev.off()

pdf(paste(dir, "/miQC/miQC_keep_probability_plot.pdf", sep = ""), width = 15, height = 10)
PlotMiQC(organoid1, color.by = "miQC.keep")
dev.off()



#Run MiQC Filtering
organoid1 <- subset(organoid1, miQC.keep == "keep") #118,688 cells and 36,601 genes

nrow(organoid1@assays$RNA@counts) #genes in organoid1
ncol(organoid1@assays$RNA@counts) #cells in organoid1
nrow(organoid0@assays$RNA@counts) #genes in organoid0
ncol(organoid0@assays$RNA@counts) #cells in organoid0


############ Remove Additional Outliers ########################################

#Remove Additional Outliers
#Only include cells with at least 200 genes expressed
organoid2 <- subset(organoid1, subset = nFeature_RNA > 200) #number of genes detected in each cell
#7,669 cells and 24,255 genes

#Only include cells with at least 500 molecules expressed
organoid2 <- subset(organoid2, subset = nCount_RNA > 500) #number of genes detected in each cell
#7,669 cells and 24,255 genes

#Remove the upper tail observed on the "nCount_RNA" violin plot (99.5% percentile)
ub <- quantile(organoid2[["nCount_RNA"]]$nCount_RNA, probs = 0.995) 
organoid2 <- organoid2[, organoid2[["nCount_RNA"]] < ub] #7,592 cells and 24255 genes

nrow(organoid2@assays$RNA@counts) #genes in organoid2 (36,601)
ncol(organoid2@assays$RNA@counts) #cells in organoid2 (118,094)


######## Visualize cell expression profile after QC ############################

#FeatureScatter is used to visualize feature-feature relationships
#Ensure linear relationship between nFeature_RNA and nCount_RNA
plot1 <- FeatureScatter(organoid2, feature1 = "nCount_RNA", feature2 = "percent.mt") + labs(legend = "TCW_14311_run1") & NoLegend() 
plot2 <- FeatureScatter(organoid2, feature1 = "nCount_RNA", feature2 = "nFeature_RNA") +labs(legend = "TCW_14311_run1") & NoLegend() 
pdf(paste(dir, "/miQC/QC_metrics_FP_after_cutoffs.pdf", sep = ""), width = 10, height = 5)
plot1 + plot2 
dev.off()


###### Now, let's examine how many cells were recovered: #######################
# Get the list of unique genotypes
genotypes <- unique(organoid2@meta.data$genotype)

# Initialize a data frame to store the counts
cell_counts_table <- data.frame(
  Genotype = character(),
  Cells_Before_QC = integer(),
  Cells_After_QC = integer(),
  Percent_Recovered = double(),
  stringsAsFactors = FALSE
)

# Loop over each genotype to get gene counts before and after QC
for (geno in genotypes) {
  # Subset for the current genotype
  subset_organoid0 <- subset(organoid0, subset = genotype == geno)
  subset_organoid2 <- subset(organoid2, subset = genotype == geno)
  
  # Count the number of cells before QC
  cells_before_QC <- length(subset_organoid0@meta.data$nFeature_RNA)
  
  # Count the number of genes after QC
  cells_after_QC <- length(subset_organoid2@meta.data$nFeature_RNA)
  
  # Calculate percent recovered
  percent_recovered <- ifelse(cells_before_QC > 0, (cells_after_QC / cells_before_QC) * 100, 0)
  
  # Add the counts to the data frame
  cell_counts_table <- rbind(
    cell_counts_table,
    data.frame(
      Genotype = geno,
      Cells_Before_QC = cells_before_QC,
      Cells_After_QC = cells_after_QC,
      Percent_Recovered = percent_recovered
    )
  )
}

# Print the table
print(cell_counts_table)

######### save organoid2 object before redeclaring it as organoid3.rds #########
saveRDS(organoid2,file.path(out,'organoid3.rds'))

################################################################################
### Remove Ambient RNA #########################################################
################################################################################

organoid3 <- readRDS('outputs/miQC/organoid3.rds')
organoid4 <- organoid3
rm(organoid4)
### Run SoupX Here ####
#SoupX 
library(tools)
library(Seurat)
library(data.table)
library(ggplot2)
library(SingleCellExperiment)
library(fishpond)
library(scater)
library(org.Hs.eg.db)
library(Matrix)
library(biomaRt)
library(SoupX)
organoid3@meta.data$
#### Scale, Normalize, Cluster data:
#normalize data
TCW <- CellCycleScoring(organoid3, s.features = cc.genes.updated.2019$s.genes,
                        g2m.features = cc.genes.updated.2019$g2m.genes,
                        set.ident = TRUE,
                        search=TRUE)

TCW[["CC.Difference"]] <- TCW$G2M.Score - TCW$S.Score

all.genes <- rownames(TCW)
TCW <- NormalizeData(TCW, normalization.method = "LogNormalize", scale.factor = 10000, assay = "RNA")

#find variable features
TCW <- FindVariableFeatures(TCW, selection.method = "vst", nfeatures = 2000)
# 
# #scale only mt
# tcw <- ScaleData(TCW, assay = "RNA", vars.to.regress = "percent.mt") #visualization scaling
# tcw  <- RunPCA(tcw, features = VariableFeatures(object = tcw))
### This is for visualizing the before/after cell cycle sorting.
pdf(paste(dir, "/pca_dim_1&2_after_first_scaling.pdf", sep = ""), width = 10, height = 10)
VizDimLoadings(tcw, dims = 1:2, reduction = "pca")
dev.off()

pdf(paste(dir, "/dim_heatmap_after_first_scaling.pdf", sep = ""), width = 10, height = 10)
DimHeatmap(tcw, dims = 1, cells = 500, balanced = TRUE)
dev.off()

pdf(paste(dir, "/ridge_plot_cell_cycle_before_scaling.pdf", sep = ""), width = 10, height = 10)
RidgePlot(tcw, features = c("PCNA", "TOP2A", "MCM6", "MKI67"), ncol = 2)
dev.off()

pdf(paste(dir, "/pca_plot_cell_cycle_before_scaling.pdf", sep = ""), width = 10, height = 10)
tcw <- RunPCA(tcw, features = c(cc.genes.updated.2019$s.genes, cc.genes.updated.2019$g2m.genes))
DimPlot(tcw)
dev.off()

#scale all
# TCW <- ScaleData(TCW, features = all.genes, assay = "RNA", vars.to.regress = c("percent.mt", "nCount_RNA", "CC.Difference"), verbose = FALSE)

#BiocManager::install("glmGamPoi")
#TCW <- SCTransform(TCW, method = "glmGamPoi", vars.to.regress = c("percent.mt", "nCount_RNA", "CC.Difference"), verbose = FALSE)

TCW <- ScaleData(TCW, assay = "RNA", verbose = FALSE)
TCW <- RunPCA(TCW, verbose = FALSE)
TCW <- RunUMAP(TCW, dims = 1:50, verbose = FALSE)
TCW <- FindNeighbors(TCW, dims = 1:50, verbose = FALSE)
TCW <- FindClusters(TCW, resolution = 0.7, verbose = FALSE)



#scale all
# seurat_obj <- ScaleData(seurat_obj, assay = "RNA", verbose = FALSE)
# 
# 
# seurat_obj <- RunPCA(seurat_obj, verbose = FALSE)
# seurat_obj <- RunUMAP(seurat_obj, dims = 1:50, verbose = FALSE)
# seurat_obj <- FindNeighbors(seurat_obj, dims = 1:50, verbose = FALSE)
# seurat_obj <- FindClusters(seurat_obj, resolution = 0.5, verbose = FALSE)





####### Beginning of SoupX Correction Code HERE ################################
####### Begin by reading in RDS file            ################################



## initialize a dataframe which will contain the before/after SoupX correction
## cell counts
SoupXSummary <- data.frame(
  sampleID = character(24),
  before_SoupX = numeric(24),
  after_SoupX = numeric(24),
  percent_removed = numeric(24),
  stringsAsFactors = FALSE
)

organoid3 <- TCW
##### DEBUGGING
sample_ids <- unique(TCW@meta.data$sample)
# sample_ids
seurat_objects <- SplitObject(TCW, split.by = "sample")



rm(TCW)
source <- '/projectnb/tcwlab/LabMember/akg/APOE_Jorganoid_project/outputs/'

### Function to convert ENSEMBL_ID names to gene symbols
getGeneSymbols <- function(ensembl_ids) {
  require(org.Hs.eg.db)
  gene_symbols <- mapIds(org.Hs.eg.db,
                         keys = ensembl_ids,
                         column = "SYMBOL",
                         keytype = "ENSEMBL",
                         multiVals = "first")
  return(gene_symbols)
}


# Load raw counts
custom_format <- list("counts" = c("U","S","A"))

# Function to apply SoupX correction
apply_soupX_correction <- function(seurat_obj, raw_counts_directory) {
  
  # # set a "raw counts" count for each of the individual Seurat samples' objects
  # seurat_obj[["original_counts"]] <- CreateAssayObject(counts = seurat_obj@assays$RNA@counts) #raw counts
  # 
  sce <- fishpond::loadFry(raw_counts_directory,
                           outputFormat = custom_format)
  # Replace Ensembl IDs in row names with gene symbols
  ensembl_ids_sce <- rownames(sce)
  gene_symbols_sce <- getGeneSymbols(ensembl_ids_sce)
  
  # if no gene symbol found, keep original ensembl id
  updated_rownames <- ifelse(is.na(gene_symbols_sce[ensembl_ids_sce]), ensembl_ids_sce, gene_symbols_sce[ensembl_ids_sce])
  rownames(sce) <- updated_rownames
  
  #create the seurat object by filtering for selected cells
  S1Counts <- counts(sce)
  
  raw <- S1Counts
  
  seurat_obj[["original_counts"]] <- CreateAssayObject(counts = seurat_obj@assays$RNA@counts) #raw counts
  #SoupX Correction
  #modify genes in raw
  genes <- intersect(rownames(seurat_obj), rownames(raw))
  length(genes)
  raw <- raw[genes,]
  
  #run SoupX algo (this package is what I ended up using instead of DecontX)
  sc <- SoupChannel(raw, seurat_obj@assays$RNA@counts[genes,])
  sc <- setClusters(sc, seurat_obj$seurat_clusters)
  sc <- autoEstCont(sc, doPlot = FALSE)
  soup_out <- adjustCounts(sc, roundToInt = TRUE) #adjusted counts matrix
  
  seurat_obj[["SoupX_counts"]] <- CreateAssayObject(counts = soup_out)
  seurat_obj@assays$RNA@counts <- soup_out #USED
  
  #calculate percent of mRNA removed
  percent_removed_soup <- (1 - (sum(seurat_obj@assays$RNA@counts) / sum(seurat_obj@assays$original_counts@counts)))*100
  seurat_obj[["percent_mRNA_removed"]] <- percent_removed_soup
  
  return(seurat_obj)
}







# Directories for raw counts of each sample
raw_counts_directories <- c("1_S16/af_quant/", "2_S8/af_quant/", "3_S9/af_quant/", "4_S14/af_quant/", "6_S23/af_quant/",
                            "8_S17/af_quant/", "9_S21/af_quant/", "10_S7/af_quant/", "11_S1/af_quant/", "12_S20/af_quant/", "13_S11/af_quant/",
                            "14_S13/af_quant/", "15_S5/af_quant/", "16_S24/af_quant/", "17_S19/af_quant/", "18_S22/af_quant/",
                            "19_S3/af_quant/", "20_S15/af_quant/", "21_S6/af_quant/", "22_S2/af_quant/", "23_S12/af_quant/", "24_S10/af_quant/",
                            "25_S18/af_quant/", "26_S4/af_quant/") # paths to quant files


# Process each Seurat object
for (i in seq_along(seurat_objects)) {
  sampleID <- names(seurat_objects)[i]
  print(sampleID)
  raw_counts_dir <- paste(source, raw_counts_directories[i], sep = '') 
  print(raw_counts_dir)
  original_cell_count <- ncol(seurat_objects[[sampleID]]@assays$RNA@counts)
  
  seurat_objects[[sampleID]] <- apply_soupX_correction(seurat_objects[[sampleID]], raw_counts_dir)
  
  after_cell_count <- ncol(seurat_objects[[sampleID]]@assays$RNA@counts)
  percent_removed <- (1 - (after_cell_count / original_cell_count)) * 100
  
  SoupXSummary[i, ] <- c(sampleID, original_cell_count, after_cell_count, percent_removed)
}


### this is to get the total mRNAs, removed mRNAs, percent removed mRNAs from each sample
for (i in seq_along(seurat_objects)) {
  sampleID <- names(seurat_objects)[i]
  print(sampleID)
  print(unique(seurat_objects[[i]]@meta.data$percent_mRNA_removed))
  
  print(sum(seurat_objects[[i]]@assays$original_counts@counts))
  print((sum(seurat_objects[[i]]@assays$RNA@counts)))
}


TCW <- merge(seurat_objects[[1]], y = c(seurat_objects[2:length(seurat_objects)]),
                   project = "APOE_Jorganoid"
)

rm(TCW, x, organoid3)

saveRDS(organoid5,file.path(out,'organoid5.rds'))



########## Everything below this point is from Deepti's original code #########
######### As I progress through this, I will be deleting and editing parts of ##
### this code.



#preserve original counts prior to ambient RNA removal
organoid3[["original_counts"]] <- CreateAssayObject(counts = organoid3@assays$RNA@counts) #raw counts

dir <- "/projectnb/tcwlab/LabMember/akg/APOE_Jorganoid_project/outputs/02-Final_QC/"

### Removal of Homotypic Doublets ##############################################
### First, Scale and Normalize the Data ########################################
TCW <- CellCycleScoring(TCW, s.features = cc.genes.updated.2019$s.genes,
                        g2m.features = cc.genes.updated.2019$g2m.genes,
                        set.ident = TRUE,
                        search=TRUE)

TCW[["CC.Difference"]] <- TCW$G2M.Score - TCW$S.Score

all.genes <- rownames(TCW)
TCW <- NormalizeData(TCW, normalization.method = "LogNormalize", scale.factor = 10000, assay = "RNA")

#find variable features
TCW <- FindVariableFeatures(TCW, selection.method = "vst", nfeatures = 2000)

#scale only mt
tcw <- ScaleData(TCW, assay = "RNA", vars.to.regress = "percent.mt") #visualization scaling
tcw  <- RunPCA(tcw, features = VariableFeatures(object = tcw))

pdf(paste(dir, "/pca_dim_1&2_after_first_scaling.pdf", sep = ""), width = 10, height = 10)
VizDimLoadings(tcw, dims = 1:2, reduction = "pca")
dev.off()

pdf(paste(dir, "/dim_heatmap_after_first_scaling.pdf", sep = ""), width = 10, height = 10)
DimHeatmap(tcw, dims = 1, cells = 500, balanced = TRUE)
dev.off()

pdf(paste(dir, "/ridge_plot_cell_cycle_before_scaling.pdf", sep = ""), width = 10, height = 10)
RidgePlot(tcw, features = c("PCNA", "TOP2A", "MCM6", "MKI67"), ncol = 2)
dev.off()

pdf(paste(dir, "/pca_plot_cell_cycle_before_scaling.pdf", sep = ""), width = 10, height = 10)
tcw <- RunPCA(tcw, features = c(cc.genes.updated.2019$s.genes, cc.genes.updated.2019$g2m.genes))
DimPlot(tcw)
dev.off()

saveRDS(TCW,file.path(out,'organoid5.rds'))
TCW <- readRDS('outputs/simpleafSeurat/organoid5.rds')

# HERE, this gave me trouble, so I'm separating out the functions into three separate regressions
#scale all
# TCW <- ScaleData(TCW, features = all.genes, assay = "RNA", vars.to.regress = c("percent.mt", "nCount_RNA", "CC.Difference"), verbose = FALSE)
all.genes <- rownames(TCW)
# Scale data regressing out 'percent.mt'
TCW <- ScaleData(TCW, features = all.genes, assay = "RNA", vars.to.regress = "percent.mt", verbose = TRUE)

# Scale data regressing out 'nCount_RNA'
TCW <- ScaleData(TCW, features = all.genes, assay = "RNA", vars.to.regress = "nCount_RNA", verbose = TRUE)

# Scale data regressing out 'CC.Difference'
TCW <- ScaleData(TCW, features = all.genes, assay = "RNA", vars.to.regress = "CC.Difference", verbose = TRUE)

####### Run it with these 3 functions having been run first thing after lab meeting #####

TCW <- ScaleData(TCW)

RunUMAP(TCW, dims = 1:50)

rm(organoid_subset, plot1, plot2, sce)

DimPlot(TCW, reduction = "umap", group.by = c("sample", "genotype"))

TCW@assays$

####### Batch Correction necessary? 02-02-2024 #################################
####### Batch correction deeemed unnecessary ###################################
#Harmony
install.packages("ggthemes")
install.packages("harmony") #devtools::install_github("immunogenomics/harmony", build_vignettes=TRUE)
library(harmony)
library(ggthemes)
#perform Harmony integration: remove bias by individual
TCW <- RunHarmony(TCW, group.by.vars = "Individual", assay.use = "SCT", reduction.save = "harmony")
#Note: for some reason, "RunHarmony" doesn't work on SCC -> run locally & download harmonized obj.#converged after 2 iterations

nn_list = c(10, 20, 30, 40, 50)
res_list = c(0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7)

for (x in nn_list) {
  for (y in res_list) {
    
    #cluster
    TCW <- TCW %>%
      RunUMAP(reduction = "harmony", n.neighbors = x, dims = 1:50) %>% 
      FindNeighbors(reduction = "harmony", dims = 1:50) %>%
      FindClusters(resolution = y) 
    
    #umap visualizations
    after <- DimPlot(TCW, reduction = 'umap', group.by ='APOE_Genotype', cols = APOE_colors)
    
    pdf(file = paste(dir, "/Harmony_Clustering/UMAP_by_APOE_genotype_after_harmony_nn_",x,"_res_",y,"_.pdf", sep = ""))
    print(after)
    dev.off()
    
    after <- DimPlot(TCW, reduction = 'umap', group.by ='Individual', cols = Individual_colors)
    
    pdf(file = paste(dir, "/Harmony_Clustering/UMAP_by_Individual_after_harmony_nn_",x,"_res_",y,"_.pdf", sep = ""))
    print(after)
    dev.off()
    
    after <- DimPlot(TCW, reduction = 'umap', group.by ='seurat_clusters')
    
    pdf(file = paste(dir, "/Harmony_Clustering/UMAP_by_clusters_after_harmony_nn_",x,"_res_",y,"_.pdf", sep = ""))
    print(after)
    dev.off()
    
  }
}

#USED: res = 0.5, nn = 30.








##############################


#run sctransform
#BiocManager::install("glmGamPoi") ### H E R E is where we stopped for now.
library(glmGamPoi)

TCW <- SCTransform(TCW, method = "glmGamPoi", vars.to.regress = c("percent.mt", "nCount_RNA", "CC.Difference"), verbose = FALSE)


TCW <- FindVariableFeatures(TCW)
TCW <- ScaleData(TCW)
TCW <- RunPCA(TCW, verbose = FALSE)
TCW <- RunUMAP(TCW, dims = 1:35, verbose = FALSE)
TCW <- FindNeighbors(TCW, dims = 1:35, verbose = FALSE)
TCW <- FindClusters(TCW, resolution = 0.7, verbose = FALSE)


# split object (by samples)


rm(reffy, reffo, ref_assay, ref_metadata)

###### Removal the doublets ####################################################
#DoubletFinder
remotes::install_github('chris-mcginnis-ucsf/DoubletFinder')
library(DoubletFinder)


##### Fix this ...

#pK identification
sweep.res.list_TCW <- paramSweep(TCW, PCs = 1:50, sct = FALSE)
sweep.stats_TCW <- summarizeSweep(sweep.res.list_TCW, GT = FALSE)
bcmvn_TCW <- find.pK(sweep.stats_TCW)

pdf(paste(dir, "/DoubletFinder_pk_param_ident_plot.pdf", sep = ""), width = 10, height = 10) 
ggplot(bcmvn_TCW, aes(pK, BCmetric, group = 1)) +
  geom_point() +
  geom_line()
dev.off()

pK <- bcmvn_TCW %>% # select the pK that corresponds to max bcmvn to optimize doublet detection
  filter(BCmetric == max(BCmetric)) %>%
  select(pK) 
pK <- as.numeric(as.character(pK[[1]])) #optimal pK = 0.23

#pN = 0.25 (default)
#exp = 12% (high throughput) or 24% (standard) #USED 12

#Homotypic Doublet Estimation
annotations <- TCW@meta.data$SoupX_seurat_clustering_0.6
homotypic.prop <- modelHomotypic(annotations)           ## ex: annotations <- seu_kidney@meta.data$ClusteringResults
nExp_poi <- round(0.12*nrow(TCW@meta.data))  ## 0.12 doublet rate?
nExp_poi.adj <- round(nExp_poi*(1-homotypic.prop))

#run doublet finder
TCW_doubletfinder <- doubletFinder_v3(TCW, 
                                      PCs = 1:50, 
                                      pN = 0.25, 
                                      pK = pK, 
                                      nExp = nExp_poi.adj,
                                      reuse.pANN = FALSE, sct = FALSE)

TCW$DF.classifications_0.25_0.23_841 <- TCW_doubletfinder$DF.classifications_0.25_0.23_841

DF_doublets <- subset(TCW_doubletfinder, DF.classifications_0.25_0.23_841 %in% "Doublet") #841




