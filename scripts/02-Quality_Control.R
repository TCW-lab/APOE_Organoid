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
.libPaths()
source("~/.Rprofile")

# Step 1: Read in the indiv samples' Seurat object files & combine into one ####

###* denotes changes made by AKG
###* Then, combine them into one Seurat object. This object will be used to
###* be a point of comparison for before/after quality control steps
setwd("/projectnb/tcwlab/LabMember/akg/projects/APOE_Organoid")

dir <- "/projectnb/tcwlab/LabMember/akg/projects/APOE_Organoid"
out<-'outputs/02-Quality-Control'


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

saveRDS(organoid0, file = paste0(dir, '/outputs/02-Quality-Control/organoid0.rds'))
unique(organoid0@meta.data$sample)

# Step 2. Iterate through organoid object & remove outliers                 ####
# By this, I mean subset for genes found in >=3 genes, 
# 1000 < number of genes in a cell < 99.5th percentile,
# Additionally, subset to remove cells with > 5% percent.mt DNA

organoid0 <- readRDS('./outputs/organoid0.rds')



# rm(S1, S2, S3, S4, S6, S8, S9, S10,
#    S11, S12, S13, S14, S15, S16, S17, S18, S19,
#    S20, S21, S22, S23, S24, S25, S26)

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

features <- c("nFeature_RNA", "nCount_RNA", "percent.mt")


## 
# First, for plotting reasons, calculate the cutoffs for nCount_RNA:
# This time, we'll do a more-stringent upper_cell_cutoff.
# Note that the reason for the last one being 0.997 was due to my having run through
# it as 0.997 the first time, then recovering and replicating the same results at that cutoff.
# We will do a 0.97 cutoff and use that this time around.
upper_cutoff_rank <- quantile(organoid0@meta.data$nCount_RNA, probs = 0.97)

lower_cell_cutoff <- 900 #minimum of the lower_cell_cutoffs prev. used.
# Manually set the lower cutoff

cat('lower_cell_cutoff:', lower_cell_cutoff, '\n')
sorted_counts <- sort(organoid0@meta.data$nCount_RNA, decreasing = FALSE)

# Calculate the index for the 97th percentile
index <- ceiling(0.97 * length(sorted_counts))

# Retrieve the exact nCount_RNA value at the 99.5th percentile
upper_cell_cutoff <- sorted_counts[index]

cat('upper_cell_cutoff:', upper_cell_cutoff, '\n')



## plot the before-QC metrics
vps<-VlnPlot(object = organoid0, features = features, ncol = 3, layer = "counts", pt.size = 0)
vps
vp1<-vps[[1]] +geom_hline(yintercept = 200,colour="black", na.rm = TRUE) +labs(x = "organoid")
vp2<-vps[[2]]+geom_hline(yintercept = upper_cell_cutoff,colour="black") +geom_hline(yintercept = lower_cell_cutoff,colour="black") +labs(x = "organoid") + scale_y_continuous(limits=c(0,100000))#+ scale_y_continuous(trans='log10')
vp3<-vps[[3]]+geom_hline(yintercept = 5,colour="black") +labs(x = "organoid")
p <- vp1 +  theme(axis.text.x = element_blank()) +vp2 + theme(axis.text.x = element_blank()) +vp3 + theme(axis.text.x = element_blank())
plot(p)

organoid0[["percent.ribo"]] <- PercentageFeatureSet(organoid0, pattern = "^RP[SL]")

## Do the subsetting to remove cells with >5% percent.mt:

organoid0 <- subset(organoid0, subset = percent.mt < 5)

## Additionally label the ribosomal percentage.
organoid0[["percent.ribo"]] <- PercentageFeatureSet(organoid0, pattern = "^RP[SL]")

## Now, we go through the miQC again, but this time with a 0.97 posterior cutoff
## We are running through this a second time because we were originally just doing it 


  
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
  


# Initialize a data frame to store cell counts before and after QC
cell_counts <- data.frame(sample = character(), 
                          cells_before_QC = integer(), 
                          cells_after_QC = integer())

# Get unique sample identifiers
samples <- unique(organoid0@meta.data$sample)
# unique(samples)
organoid_list <- list()

VlnPlot(organoid0, pt.size = 0, group.by = 'sample', features = 'nCount_RNA', log = TRUE)




# Loop over each sample
for (sample_id in samples) { # samples vector defined at beginning of code
  # Subset Seurat object for the current sample
  indiv_QC <- subset(organoid0, subset = sample == sample_id)
  
  # record the number of cells before QC
  cells_before <- ncol(indiv_QC)
  
  # only include cells with at least 200 genes expressed
  indiv_QC <- subset(indiv_QC, subset = nFeature_RNA > 200)
  num_umis_before <- sum(indiv_QC[["nCount_RNA"]] < 1000)
  
  # Filter genes expressed in less than 3 cells
  raw_counts_mat <- indiv_QC@assays$RNA@counts
  genes_to_keep <- Matrix::colSums(raw_counts_mat) >= 3
  indiv_QC <- subset(indiv_QC, features = rownames(indiv_QC)[genes_to_keep])
  
  
  ncol(indiv_QC)
  # Calculate upper limit, set at 99.5 percentile
  head(sort(unique(indiv_QC@meta.data$nCount_RNA), decreasing = TRUE))
  upper_limit <- quantile(indiv_QC@meta.data$nCount_RNA, probs = 0.995)
  print(upper_limit)
  indiv_QC1 <- CalculateBarcodeInflections(
    indiv_QC,
    barcode.column = "nCount_RNA",
    group.column = "sample",
    threshold.low = NULL, #this is barcode rank, not absolute counts.
    threshold.high = upper_limit
  )
  indiv_QC2 <- SubsetByBarcodeInflections(object = indiv_QC1)
  
  # Record the number of cells after QC
  cells_after <- ncol(indiv_QC2)
  print(cells_before)
  print(cells_after)
  # append the QC'd object onto the organoid_list (so we can combine after)
  organoid_list[[sample_id]] <- indiv_QC2
  
  # Append the cell count information to the data frame
  cell_counts <- rbind(cell_counts, data.frame(sample = sample_id,
                                               cells_before_QC = cells_before,
                                               cells_after_QC = cells_after))
  
  num_umis_after <- sum(indiv_QC2[["nCount_RNA"]] < 1000)
  
  # Print statement
  print(paste("Number of UMIs per cell before (under 1000):", num_umis_before, 
              ". Number of UMIs per cell after (under 1000):", num_umis_after))
  
}


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

aggregated_cell_counts <- aggregate(. ~ genotype, data = cell_counts, FUN = sum)

# View the cell counts data frame
#print(cell_counts)
print(aggregated_cell_counts)





organoid1 <- merge(x = organoid_list[[1]], y = organoid_list[-1], 
                   add.cell.ids = names(organoid_list[-1]), 
                   project = "APOE_Organoid")

library(stringr)
sum(str_detect(rownames(organoid0), '^ENSG'))

# View the cell counts data frame
print(cell_counts)
