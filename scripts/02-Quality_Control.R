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

unique(organoid0@meta.data$genotype)

# saveRDS(organoid0, file = paste0(dir, '/outputs/02-Quality-Control/organoid0.rds'))
unique(organoid0@meta.data$sample)

# Step 2. Iterate through organoid object & remove outliers                 ####
# By this, I mean subset for genes found in >=3 genes, 
# 1000 < number of genes in a cell < 99.5th percentile,
# Additionally, subset to remove cells with > 5% percent.mt DNA


rm(S1, S2, S3, S4, S6, S8, S9, S10,
   S11, S12, S13, S14, S15, S16, S17, S18, S19,
   S20, S21, S22, S23, S24, S25, S26)

saveRDS('.outputs/02-Quality-Control/organoid0.rds')


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


# Subset mt_gene_names to include only genes present in the Seurat object
valid_mt_gene_names <- mt_gene_names[mt_gene_names %in% rownames(organoid0[["RNA"]]$counts)]

# Add genes that start with "mt-" or "MT-" to the list
additional_mt_genes <- rownames(organoid0[["RNA"]]$counts)[grepl("^mt-", rownames(organoid0[["RNA"]]$counts), ignore.case = TRUE)]

# Combine the two lists and remove duplicates
valid_mt_gene_names <- unique(c(valid_mt_gene_names, additional_mt_genes))

# Recalculate the percentage of mitochondrial genes
organoid0[["percent.mt"]] <- PercentageFeatureSet(organoid0, features = valid_mt_gene_names)

### Now that we have our percent.mt labeled, take a snapshot of what the data
## looks like across all 24 samples, in terms of nCount_RNA, nFeature_RNA, and
## percent.mt

## Additionally label the ribosomal percentage.
organoid0[["percent.ribo"]] <- PercentageFeatureSet(organoid0, pattern = "^RP[SL]")




# saveRDS(organoid0, file = paste0(dir, out, 'organoid0.rds'))
# organoid0 <- readRDS(paste0(dir, out, 'organoid0.rds'))

## 
# First, for plotting reasons, calculate the cutoffs for nCount_RNA:
# This time, we'll do a more-stringent upper_cell_cutoff.
# Note that the reason for the last one being 0.997 was due to my having run through
# it as 0.997 the first time, then recovering and replicating the same results at that cutoff.
# We will do a 0.99 cutoff and use that this time around.
upper_cutoff_rank <- quantile(organoid0@meta.data$nCount_RNA, probs = 0.99)

lower_cell_cutoff <- 900 #minimum of the lower_cell_cutoffs prev. used.
# Manually set the lower cutoff

cat('lower_cell_cutoff:', lower_cell_cutoff, '\n')
sorted_counts <- sort(organoid0@meta.data$nCount_RNA, decreasing = FALSE)

# Calculate the index for the 99th percentile
index <- ceiling(0.97 * length(sorted_counts))

# Retrieve the exact nCount_RNA value at the 99.5th percentile
upper_cell_cutoff <- sorted_counts[index]

cat('upper_cell_cutoff:', upper_cell_cutoff, '\n')

# Extract the nFeature_RNA values from the metadata
nFeature_RNA_values <- organoid0@meta.data$nFeature_RNA

# Sort the nFeature_RNA values in ascending order
sorted_nFeature_RNA <- sort(nFeature_RNA_values, decreasing = FALSE)

# Calculate the index for the 99th percentile
index_nFeature <- ceiling(0.97 * length(sorted_nFeature_RNA))

# Retrieve the exact nFeature_RNA value at the 99th percentile
upper_feature_cutoff <- sorted_nFeature_RNA[index_nFeature]

# Print the upper cutoff for nFeature_RNA
cat('upper_feature_cutoff:', upper_feature_cutoff, '\n')


features <- c("nFeature_RNA", "nCount_RNA", "percent.mt")
## plot the before-QC metrics
vps<-VlnPlot(object = organoid0, features = features, ncol = 3, layer = "counts", pt.size = 0)
vps
vp1<-vps[[1]] +geom_hline(yintercept = 500,colour="black", na.rm = TRUE) + geom_hline(yintercept = upper_gene_cutoff, color = 'black', na.rm = TRUE) +labs(x = "organoid")
vp2<-vps[[2]]+geom_hline(yintercept = upper_cell_cutoff,colour="black") +geom_hline(yintercept = lower_cell_cutoff,colour="black", na.rm = TRUE) +labs(x = "organoid") #+ scale_y_continuous(limits=c(0,100000))#+ scale_y_continuous(trans='log10')
vp3<-vps[[3]]+geom_hline(yintercept = 12.5,colour="black") +labs(x = "organoid")
p <- vp1 +  theme(axis.text.x = element_blank()) +vp2 + theme(axis.text.x = element_blank()) +vp3 + theme(axis.text.x = element_blank())
plot(p)
ggsave(filename = "outputs/images/02-Quality-Control/combined_Vln_plots_before_QC.png", 
       plot = p, width = 6, height = 4)

## Now, we go through the miQC again, but this time with a 0.97 posterior cutoff
## We are running through this a second time because we were originally just doing it 


  

  


# Initialize a data frame to store cell counts before and after QC
cell_counts <- data.frame(sample = character(), 
                          genotype = character(),
                          cells_before_QC = integer(), 
                          cells_after_QC = integer()
                          )

# Get unique sample identifiers
samples <- unique(organoid0@meta.data$sample)
# unique(samples)
organoid_list <- list()

VlnPlot(organoid0, pt.size = 0, group.by = 'sample', features = 'nFeature_RNA', log = TRUE)

# organoid0 <- readRDS('outputs/organoid0.rds')
print(unique(organoid0@meta.data$sample))

#### Count initial sample's cell counts here.
# Load required libraries


### Then do the QC here on the organoid0 object

# Loop over each sample
for (sample_id in samples) {
  # Subset Seurat object for the current sample
  indiv_QC <- subset(organoid0, subset = sample == sample_id)
  
  # Record the number of cells before QC
  cells_before <- ncol(indiv_QC)
  
  # Add cells_before_QC to the meta.data of the original Seurat object
  organoid0@meta.data$cells_before_QC[organoid0@meta.data$sample == sample_id] <- cells_before
  
  # Set thresholds for nFeature_RNA
  lower_limit_nFeature <- 500
  upper_limit_nFeature <- quantile(indiv_QC@meta.data$nFeature_RNA, probs = 0.97)
  
  # Subset cells based on nFeature_RNA thresholds
  indiv_QC <- subset(indiv_QC, subset = nFeature_RNA > lower_limit_nFeature & nFeature_RNA < upper_limit_nFeature)
  
  num_umis_before <- sum(indiv_QC[["nCount_RNA"]] < 1000)
  
  # Filter genes expressed in less than 3 cells
  raw_counts_mat <- GetAssayData(indiv_QC, slot = "counts", assay = "RNA")
  genes_to_keep <- rowSums(raw_counts_mat > 0) >= 3
  indiv_QC <- subset(indiv_QC, features = rownames(indiv_QC)[genes_to_keep])
  
  # Set thresholds for nCount_RNA
  upper_limit_nCount <- quantile(indiv_QC@meta.data$nCount_RNA, probs = 0.97)
  lower_limit_nCount <- 500  # Example lower limit
  
  # Subset cells based on nCount_RNA thresholds
  indiv_QC <- subset(indiv_QC, subset = nCount_RNA > lower_limit_nCount & nCount_RNA < upper_limit_nCount)
  
  indiv_QC <- subset(indiv_QC, subset = percent.mt < 10)
  
  # Record the number of cells after QC
  cells_after <- ncol(indiv_QC)
  
  # Add cells_after_QC to the meta.data of the original Seurat object
  organoid0@meta.data$cells_after_QC[organoid0@meta.data$sample == sample_id] <- cells_after
  
  # Print statement for debugging
  print(paste0("Sample: ", sample_id, ". cells_before:", cells_before, ". cells_after:", cells_after))
  
  # Append the QC'd object onto the organoid_list
  organoid_list[[sample_id]] <- indiv_QC
  
  # Extract genotype for the current sample
  genotype <- unique(indiv_QC@meta.data$genotype)
  
  # Append the cell count information to the data frame
  cell_counts <- rbind(cell_counts, data.frame(sample = sample_id,
                                               genotype = genotype,
                                               cells_before_QC = cells_before,
                                               cells_after_QC = cells_after))
  
  num_umis_after <- sum(indiv_QC[["nCount_RNA"]] < 1000)
  
  # Print statement for debugging
  print(paste("Number of UMIs per cell before (under 1000):", num_umis_before, 
              ". Number of UMIs per cell after (under 1000):", num_umis_after))
}


### Count the after-QC cell counts here.

# merge all QC'd Seurat objects into one (all alrdy have unique cell identifiers)
if(length(organoid_list) > 1) {
  # do.call to apply the merge function to a list of Seurat objects
  organoid1 <- do.call(merge, c(list(x = organoid_list[[1]], y = organoid_list[-1]), list(project = "APOE_Organoid")))
} else if(length(organoid_list) == 1) {
  # only one object, simply rename the project if necessary
  organoid1 <- RenameProject(organoid_list[[1]], "APOE_Organoid")
} else {
  # error catch
  message("No data processed successfully.")
}

print(cell_counts)

# Aggregate cell counts by genotype
aggregated_cell_counts <- aggregate(cbind(cells_before_QC, cells_after_QC) ~ genotype, data = cell_counts, FUN = sum)

# Print the result
print(aggregated_cell_counts)

###



library(stringr)
sum(str_detect(rownames(organoid0), '^ENSG'))

# View the cell counts data frame
print(cell_counts)

saveRDS(organoid1, file = paste0(dir, "/", out, "/organoid1.rds"))

# Step 3. Visualize after-QC metrics                                             ####

organoid1 <- readRDS('/projectnb/tcwlab/LabMember/akg/projects/APOE_Organoid/outputs/02-Quality-Control/organoid1.rds')

features <- c("nFeature_RNA", "nCount_RNA", "percent.mt")
## plot the before-QC metrics
vps<-VlnPlot(object = organoid1, features = features, ncol = 3, layer = "counts", pt.size = 0)
vps
vp1<-vps[[1]] +geom_hline(yintercept = 200,colour="black", na.rm = TRUE) +labs(x = "organoid")
vp2<-vps[[2]]+geom_hline(yintercept = upper_cell_cutoff,colour="black") +geom_hline(yintercept = lower_cell_cutoff,colour="black") +labs(x = "organoid") #+ scale_y_continuous(limits=c(0,100000))#+ scale_y_continuous(trans='log10')
vp3<-vps[[3]]+geom_hline(yintercept = 10,colour="black") +labs(x = "organoid")
p <- vp1 +  theme(axis.text.x = element_blank()) +vp2 + theme(axis.text.x = element_blank()) +vp3 + theme(axis.text.x = element_blank())
plot(p)


### This below is just to re-define the quantile and upper_cell and lower_cell
### cutoffs for the sake of visuali
# We will do a 0.99 cutoff and use that this time around.
upper_cutoff_rank <- quantile(organoid1@meta.data$nCount_RNA, probs = 0.99)

lower_cell_cutoff <- 900 #minimum of the lower_cell_cutoffs prev. used.
# Manually set the lower cutoff

cat('lower_cell_cutoff:', lower_cell_cutoff, '\n')
sorted_counts <- sort(organoid1@meta.data$nCount_RNA, decreasing = FALSE)

# Calculate the index for the 99th percentile
index <- ceiling(0.99 * length(sorted_counts))

# Retrieve the exact nCount_RNA value at the 99.5th percentile
upper_cell_cutoff <- sorted_counts[index]

cat('upper_cell_cutoff:', upper_cell_cutoff, '\n')

organoid1 <- subset(organoid1, subset = nCount_RNA < upper_cell_cutoff)

VlnPlot(organoid1, pt.size = 0, group.by = 'sample', features = 'nCount_RNA', log = TRUE)
