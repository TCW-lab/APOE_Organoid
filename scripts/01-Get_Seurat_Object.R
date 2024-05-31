##############################################
## 12-01-23 akg adjusted script to separate ##
## un-spliced from spliced data             ##
##############################################
setwd('/projectnb/tcwlab/LabMember/akg/projects/APOE_Organoid')

out<-'outputs/01-Simpleaf_Outputs/Seurat'
# dir.create(out)

library(tools)
library(Seurat)
library(data.table)
library(ggplot2)
library(SingleCellExperiment)
library(fishpond)
library(scater)
library(loomR)
library(biomaRt)
library(ensembldb)
library(Matrix)
library(AnnotationDbi)
library(dplyr)
library(DropletQC)
library(DropletUtils)
library(ggplot2)
# 1. Begin by reading in mat and gene names (10X genomics changes as per Alexandre's pipeline) ####

###* Below is the QC pipeline, performed on the outputted matrices of Simpleaf
###* This was done on all 24 samples.



### Read Alexandre's cellRanger output to get gene names map in GENEID / SYMBOL format
gene_names <- fread('/projectnb/tcwlab/LabMember/adpelle1/projects/APOE_Jorganoid/outputs/CellRangerCount/sample_1/outs/filtered_feature_bc_matrix/features.tsv.gz', 
                    col.names = c('GENEID', 'SYMBOL', 'c3'))
# remove the c3 column from the gene_names data.table df
gene_names <- gene_names[, c3 := NULL]
# set GENEID to rownames attribute, in order to make the following easier
rownames(gene_names) <- gene_names$GENEID

### 2. Now that we have handled the gene_names list, we create our SCE object ####
#### Use this to open up the files, then once you have your matrix, 
#### counts(sce[,cells]) gives you the matrix. With this matrix you just change rownames
#### and then you can convert to Seurat object. 

custom_format <- list("counts" = c("U","S","A"))

sce <- fishpond::loadFry("outputs/01-Simpleaf_Outputs/1_S16/af_quant/",
                         outputFormat = custom_format)


head(colnames(sce))

length(unique(colnames(sce)))

##### Replace Ensembl IDs in row names with gene symbols ####

# Extract the rownames (ENSEMBL IDs) from the SCE object
ensembl_ids <- rownames(sce)

# Create a named vector for mapping
geneid_to_symbol <- setNames(gene_names$SYMBOL, gene_names$GENEID)

# Map the rownames (ENSEMBL IDs) to SYMBOL values
symbol_names <- geneid_to_symbol[ensembl_ids]

# If there is no matching symbol, keep the original ENSEMBL ID
symbol_names <- ifelse(is.na(symbol_names), ensembl_ids, symbol_names)

# Set the rownames of the SCE object to the mapped SYMBOL values
rownames(sce) <- make.unique(symbol_names)

# View the result
# print(rownames(sce))

#### 3. Now that we have our SCE object, convert it to a Seurat object ####


# Calculate the number of UMIs per cell
cell_counts <- colSums(counts(sce))
cells_count <- data.table(cell = names(cell_counts), count = cell_counts)
cells_count[, cell_rank := rank(-count)]



# Define cell cutoff
lower_cell_cutoff <- 900
selected_cells <- names(cell_counts[cell_counts > lower_cell_cutoff])

# Display the number of cells after filtering
filtered_cell_count <- length(selected_cells)
cat("Number of cells after filtering:", filtered_cell_count, "\n")

# Plot the 'knee' plot
knee_plot <- ggplot(cells_count, aes(x = cell_rank, y = count)) + 
  geom_line() +
  scale_x_log10() +
  scale_y_log10() +
  theme_bw() +
  geom_hline(yintercept = lower_cell_cutoff) + 
  ggtitle('sample1')
print(knee_plot)
# Save the plot
ggsave(filename = "outputs/images/01-before-QC-Knee-Plots/kp_sample1_before_QC.png", 
       plot = knee_plot, width = 8, height = 6)

# Create Seurat object with filtered cells
seurat_obj <- CreateSeuratObject(counts = counts(sce[, selected_cells]), project = 'APOE_Organoid')

length(colnames(seurat_obj))




## 4. Assign percent.mt to each cell & do before-QC plotting                ####
# First, for plotting reasons, calculate the cutoffs for nCount_RNA:
upper_cutoff_rank <- quantile(seurat_obj@meta.data$nCount_RNA, probs = 0.995)

# Manually set the lower cutoff
cat('lower_cell_cutoff:', lower_cell_cutoff, '\n')
sorted_counts <- sort(seurat_obj@meta.data$nCount_RNA, decreasing = FALSE)

# Calculate the index for the 97th percentile
index <- ceiling(0.995 * length(sorted_counts))

# Retrieve the exact nCount_RNA value at the 99.5th percentile
upper_cell_cutoff <- sorted_counts[index]

cat('upper_cell_cutoff:', upper_cell_cutoff, '\n')



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

seurat_obj@assays$RNA$counts

# Subset mt_gene_names to include only genes present in the Seurat object
valid_mt_gene_names <- mt_gene_names[mt_gene_names %in% rownames(seurat_obj[["RNA"]]$counts)]

# Recalculate the percentage of mitochondrial genes
seurat_obj[["percent.mt"]] <- PercentageFeatureSet(seurat_obj, features = valid_mt_gene_names)
# 
# # Subset the Seurat object to remove cells with percent.mt >= 5%
# seurat_obj <- subset(seurat_obj, subset = percent.mt < 5)
ncol(seurat_obj)

# Visualize the nCount_RNA, nFeature_RNA, and percent.mt for the seurat_obj pre-QC
# Define the features to plot
features <- c("nFeature_RNA", "nCount_RNA", "percent.mt")

## plot the before-QC metrics

vps<-VlnPlot(object = seurat_obj, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3, layer = "counts", pt.size = 0)
vps
vp1<-vps[[1]] +geom_hline(yintercept = 200,colour="black", na.rm = TRUE) +labs(x = "sample1")
vp2<-vps[[2]]+geom_hline(yintercept = upper_cell_cutoff,colour="black") +geom_hline(yintercept = lower_cell_cutoff,colour="black") +labs(x = "sample1") + scale_y_continuous(limits=c(0,100000))#+ scale_y_continuous(trans='log10')
vp3<-vps[[3]]+geom_hline(yintercept = 5,colour="black") +labs(x = "sample1")
p <- vp1 +  theme(axis.text.x = element_blank()) +vp2 + theme(axis.text.x = element_blank()) +vp3 + theme(axis.text.x = element_blank())
plot(p)

ggsave(filename = "outputs/images/01-before-Vln-Plots/sample1_combined_violin_plot.png", plot = p, width = 8, height = 4)


dev.off()




#### 5. Perform rudimentary QC (nCount_RNA filtering) ####
# first, Filter by nCount_RNA :

library(Seurat)

# Define upper nCount_RNA cutoff as 99.5th percentile
upper_cutoff <- upper_cell_cutoff

# Manually set the lower cutoff
lower_cutoff <- lower_cell_cutoff

# Filter cells based on the calculated inflection points
selected_cells <- WhichCells(seurat_obj, expression = nCount_RNA > lower_cutoff & nCount_RNA < upper_cutoff)

# Subset the Seurat object to keep only the selected cells
seurat_obj <- subset(seurat_obj, cells = selected_cells)

# Display the cutoffs
cat("Lower cutoff:", lower_cell_cutoff, "\n")
cat("Upper cutoff:", upper_cell_cutoff, "\n")

# Display the number of cells after filtering
cat("Number of cells after filtering (nCount_RNA filtering):", ncol(seurat_obj), "\n")



length(colnames(seurat_obj))

seurat_obj[['sample']] <- 'sample1'

saveRDS(seurat_obj, 'outputs/01-individual-Seurat-Files/S1.rds')
















#### 6. Use miQC to perform filtering based on percent.mt, nFeature_RNA     ####





library(miQC)
library(flexmix)
library(SeuratWrappers)
seurat_obj1 <- RunMiQC(seurat_obj, percent.mt = "percent.mt", nFeature_RNA = "nFeature_RNA", 
                       posterior.cutoff = 0.97, model.slot = "flexmix_model")
seurat_obj1 <- subset(seurat_obj1, miQC.keep == "keep")

ncol(seurat_obj1)




## Plot the knee plot of the nCount_RNA barcode rank (after filtering)    ##

# Extract nCount_RNA values
nCount_RNA_values <- seurat_obj$nCount_RNA

# Create a data frame for the plot
cells_count <- data.frame(cell = names(nCount_RNA_values), count = nCount_RNA_values)
cells_count$cell_rank <- rank(-cells_count$count)

# Manually define lower and upper cutoffs
lower_cutoff <- cell_cutoff
# Example value for upper_cutoff, adjust based on your data
upper_cutoff <- upper_cutoff

# Create the barcode rank plot with log10 transformation and custom y-axis title
knee_plot <- ggplot(cells_count, aes(x = cell_rank, y = count)) + 
  geom_line() +
  scale_x_log10() +
  scale_y_log10() +
  theme_bw() +
  geom_hline(yintercept = lower_cutoff, linetype = "dashed", color = "blue") +
  geom_hline(yintercept = upper_cutoff, linetype = "dashed", color = "red") +
  labs(y = "log10(nCount_RNA)", x = "Rank") +
  ggtitle('S1')
print(knee_plot)

# Save the plot
ggsave(filename = "outputs/images/01-after-QC-Knee-Plots/kp_sample1_after_QC.png", plot = knee_plot, width = 8, height = 6)

## Recorded numberCells at this point



# 7. Perform filtering based on nFeature_RNA now (>200)                    #####

seurat_obj <- subset(seurat_obj1, subset = nFeature_RNA > 200)
ncol(seurat_obj1)

#### At this point, proceed with the nFeature_RNA and then percent.mt filtering:
#### Use DropletQC on nFeature, also 200 as lower cutoff.
#### Use 5% at upper percent.mt cutoff.




seurat_obj[["percent.ribo"]] <- PercentageFeatureSet(seurat_obj, pattern = "^RP[SL]")

# Identify mitochondrial genes that start with "mt-" or "MT-"
mt_genes <- grep("^mt-|^MT-", rownames(seurat_obj), value = TRUE)

# Verify the identified mitochondrial genes
print(mt_genes)

seurat_obj[["percent.mt"]] <- PercentageFeatureSet(seurat_obj, features = mt_genes)





#here test different threshold to find the best one that is between the 2 knees, or at the middle the first knee
#we expecrt cutoff between 100 and 2000 UMIS (can be more if the sequencing depth or cell RNA amount is bigger)
#we expect to have between 2000 and 25k cells

sum(rowSums(mat)>cell_cutoff) #number of cells 

?rowSums
cells<-rownames(sce)[rowSums(sce@assays@data$counts)>cell_cutoff] 
#cellS5<-cells_count[count>cell_cutoff]$cell #same
matf<-mat[cells,]



#create the seurat object by filtering for selected cells
S1<-CreateSeuratObject(counts(sce[,cells]),project = 'APOE_Organoid')
head(rownames(S1))
# An object of class Seurat 
# 58219 features across 8979 samples within 1 assay 
# Active assay: RNA (58219 features, 0 variable features)
# 1 layer present: counts

saveRDS(S1,file.path(out,'sample1.rds'))


nrow(mat) #191 363
colnames(mat)
#get the number of umis per cell
cells_count<-data.table(cell=rownames(mat),
                        count=rowSums(mat))
cells_count[,cell_rank:=rank(-count)]
#plot the 'knee' plot
ggplot(cells_count,aes(x=cell_rank,y=count))+geom_line()+
  scale_x_log10()+scale_y_log10()+theme_bw()+
  geom_hline(yintercept = 1000) + ggtitle('sample1')



###############################################################

# Calculate QC metrics
sce <- scater::addPerCellQC(sce, subsets = list(Mito = grep("^MT-", rownames(sce))))

class(sce)
# Define QC filtering criteria
min_genes <- 1000
max_genes <- 10000
max_mito <- 5
class(sce)

####




S1<-CreateSeuratObject(t(as.matrix(matf)),project = 'Organoid')
# An object of class Seurat 
# 174657 features across 5460 samples within 1 assay 
# Active assay: RNA (174657 features, 0 variable features)
# 1 layer present: counts
S5
saveRDS(S1,file.path(out,'S1.rds'))

S1$sample<-'S1'

#did for all


#read in and then merge the Seurat objects
S1 <- readRDS(file = "./outputs/01-get_seurat_object/S1.rds")
S2 <- readRDS(file = "./outputs/01-get_seurat_object/S2.rds")
S3 <- readRDS(file = "./outputs/01-get_seurat_object/S3.rds")
S4 <- readRDS(file = "./outputs/01-get_seurat_object/S4.rds")
S5 <- readRDS(file = "./outputs/01-get_seurat_object/S5.rds")
S6 <- readRDS(file = "./outputs/01-get_seurat_object/S6.rds")
S7 <- readRDS(file = "./outputs/01-get_seurat_object/S7.rds")
S8 <- readRDS(file = "./outputs/01-get_seurat_object/S8.rds")
S9 <- readRDS(file = "./outputs/01-get_seurat_object/S9.rds")
S10 <- readRDS(file = "./outputs/01-get_seurat_object/S10.rds")
S11 <- readRDS(file = "./outputs/01-get_seurat_object/S11.rds")
S12 <- readRDS(file = "./outputs/01-get_seurat_object/S12.rds")
S13 <- readRDS(file = "./outputs/01-get_seurat_object/S13.rds")
S14 <- readRDS(file = "./outputs/01-get_seurat_object/S14.rds")
S15 <- readRDS(file = "./outputs/01-get_seurat_object/S15.rds")
S16 <- readRDS(file = "./outputs/01-get_seurat_object/S16.rds")
S17 <- readRDS(file = "./outputs/01-get_seurat_object/S17.rds")
S18 <- readRDS(file = "./outputs/01-get_seurat_object/S18.rds")
S19 <- readRDS(file = "./outputs/01-get_seurat_object/S19.rds")
S20 <- readRDS(file = "./outputs/01-get_seurat_object/S20.rds")
S21 <- readRDS(file = "./outputs/01-get_seurat_object/S21.rds")
S22 <- readRDS(file = "./outputs/01-get_seurat_object/S22.rds")
S23 <- readRDS(file = "./outputs/01-get_seurat_object/S23.rds")
S24 <- readRDS(file = "./outputs/01-get_seurat_object/S24.rds")

# install.packages('SeuratData')
library(remotes)
# remotes::install_github("hoxo-m/githubinstall")
library(githubinstall)

# githubinstall("SatijaLab/SeuratData")
library(SeuratData)


# unloadNamespace("fastmap")

# detach("package:fastmap", unload = TRUE)
# install.packages("fastmap")
# devtools::install_github('satijalab/seurat-data')



APOE_Jorganoid <- merge(S1, y = c(S2, S3, S4, S5, S6, S7, S8, S9, S10,
                            S11, S12, S13, S14, S15, S16, S17, S18, S19,
                            S20, S21, S22, S23, S24), 
                  add.cell.ids = c("S1", "S2", "S3", "S4", "S5", "S6", "S7", "S8", "S9", "S10",
                                  "S11", "S12", "S13", "S14", "S15", "S16", "S17", "S18", "S19",
                                  "S20", "S21", "S22", "S23", "S24"),
                  project = "APOE_Jorganoid"
                  )
APOE_Jorganoid <- NormalizeData(APOE_Jorganoid)
APOE_Jorganoid <- FindVariableFeatures(APOE_Jorganoid)
APOE_Jorganoid <- ScaleData(APOE_Jorganoid)

APOE_Jorganoid <- RunPCA(APOE_Jorganoid, features = VariableFeatures(object = APOE_Jorganoid))
DimPlot(APOE_Jorganoid, reduction = "pca")
ElbowPlot(object = APOE_Jorganoid)


APOE_Jorganoid <- RunUMAP(APOE_Jorganoid, dims = 1:10)
DimPlot(APOE_Jorganoid, reduction = "umap")

saveRDS(APOE_Jorganoid,file.path(out,'APOE_Organoid.rds'))

#---- end of data pre-processing ----#
# APOE_Jorganoid <- readRDS(file = "./outputs/01-get_seurat_object/APOE_Jorganoid.rds")

ggplot(APOE_Jorganoid@reductions$umap@cell.embeddings, aes(x=UMAP_1, y=UMAP_2, color=orig.ident)) + geom_point()

#--- Data Analysis for CellRanger pipeline ---#
# load in CellRanger Data


#--- end of CellRanger ---#
#--- Beginning of downstream data analysis ---#


ggplot()
