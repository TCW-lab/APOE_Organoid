##############################################
## 12-01-23 akg adjusted script to separate ##
## un-spliced from spliced data             ##
##############################################
setwd('/projectnb/tcwlab/LabMember/akg/projects/APOE_Organoid_project/')

out<-'outputs/01-Simpleaf_Outputs'
# dir.create(out)

library(tools)
library(Seurat)
library(data.table)
library(ggplot2)
library(SingleCellExperiment)
library(fishpond)
library(scater)

if (!requireNamespace("biomaRt", quietly = TRUE))
  install.packages("biomaRt")
library(biomaRt)

library(org.Hs.eg.db)
library(Matrix)

getGeneSymbols <- function(ensembl_ids) {
  require(org.Hs.eg.db)
  gene_symbols <- mapIds(org.Hs.eg.db,
                         keys = ensembl_ids,
                         column = "SYMBOL",
                         keytype = "ENSEMBL",
                         multiVals = "first")
  return(gene_symbols)
}
 ####* 05-10-2024####
counts.mat <- organoid0@assays$RNA@counts

cr_features <- fread('/projectnb/tcwlab/LabMember/adpelle1/projects/APOE_Jorganoid/outputs/CellRangerCount/sample_1/outs/filtered_feature_bc_matrix/features.tsv.gz', col.names = c('ENSEMBL', 'symbol', 'expression'))

rownames(counts.mat) = cr_features[rownames(counts.mat), on = 'ENSEMBL']$symbol

rownames(counts.mat)


ensembl_ids_sce <- rownames(sce)
gene_symbols_sce <- getGeneSymbols(ensembl_ids_sce)

# Replace Ensembl IDs with gene symbols in row names
# If no gene symbol is found, retain the original Ensembl ID
updated_rownames <- ifelse(is.na(gene_symbols_sce[ensembl_ids_sce]), ensembl_ids_sce, gene_symbols_sce[ensembl_ids_sce])
rownames(sce) <- updated_rownames

micromamba create -n simpleaf -y -c bioconda -c conda-forge simpleaf piscem
micromamba activate simpleaf

micromamba create -n simpleaf -y -c bioconda -c conda-forge simpleaf piscem
micromamba activate simpleaf



###* 10X genomics changes as per Alexandre's pipeline              *############

###* Below is the QC pipeline, performed on the outputted matrices of Simpleaf
###* This was done on all 24 samples. The 
# Read the matrix
mat <- ReadMtx('outputs/01-Simpleaf_Outputs/12_S20/af_quant/alevin/quants_mat.mtx',
               cells = 'outputs/01-Simpleaf_Outputs/12_S20/af_quant/alevin/quants_mat_cols.txt',
               features = 'outputs/01-Simpleaf_Outputs/12_S20/af_quant/alevin/quants_mat_rows.txt',
               feature.column = 1)

# Use getGeneSymbols function to map Ensembl IDs to gene symbols
ensembl_ids <- colnames(mat)
gene_symbols <- getGeneSymbols(ensembl_ids)

# Replace Ensembl IDs with gene symbols in column names
colnames(mat) <- unname(gene_symbols[ensembl_ids])

# Check the updated column names
head(colnames(mat))

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

################################################ delete this below after converting all 
############################################### ENSEMBL to symbols
cell_cutoff=1000
cells<-rownames(mat)[rowSums(mat)>cell_cutoff] 
custom_format <- list("counts" = c("U","S","A"))

sce <- fishpond::loadFry("outputs/1_S16/af_quant/",
                         outputFormat = custom_format)

# Assuming 'sce' is your Single Cell Experiment object
# Replace Ensembl IDs in row names with gene symbols
ensembl_ids_sce <- rownames(sce)
gene_symbols_sce <- getGeneSymbols(ensembl_ids_sce)

# Replace Ensembl IDs with gene symbols in row names
# If no gene symbol is found, retain the original Ensembl ID
updated_rownames <- ifelse(is.na(gene_symbols_sce[ensembl_ids_sce]), ensembl_ids_sce, gene_symbols_sce[ensembl_ids_sce])
rownames(sce) <- updated_rownames

#create the seurat object by filtering for selected cells
S1<-CreateSeuratObject(counts(sce[,cells]),project = 'Organoid')
head(rownames(S1))
# An object of class Seurat 
# 58219 features across 8979 samples within 1 assay 
# Active assay: RNA (58219 features, 0 variable features)
# 1 layer present: counts

saveRDS(S1,file.path(out,'sample1.rds'))


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


#get the number of umis per cell
cells_count<-data.table(cell=rownames(mat),
                        count=rowSums(mat))
cells_count[,cell_rank:=rank(-count)]
#plot the 'knee' plot
ggplot(cells_count,aes(x=cell_rank,y=count))+geom_line()+
  scale_x_log10()+scale_y_log10()+theme_bw()+
  geom_hline(yintercept = min_genes) + geom_hline(yintercept = max_genes) + ggtitle('S15')
#here test different threshold to find the best one that is between the 2 knees, or at the middle the first knee
#we expecrt cutoff between 100 and 2000 UMIS (can be more if the sequencing depth or cell RNA amount is bigger)
#we expect to have between 2000 and 25k cells

sum(rowSums(mat)>cell_cutoff) #number of cells 

cells<-rownames(mat)[rowSums(mat)>cell_cutoff] 
#cellS5<-cells_count[count>cell_cutoff]$cell #same
matf<-mat[cells,]

S5<-CreateSeuratObject(t(as.matrix(matf)),project = 'Organoid')
# An object of class Seurat 
# 174657 features across 5460 samples within 1 assay 
# Active assay: RNA (174657 features, 0 variable features)
# 1 layer present: counts
S5
saveRDS(S5,file.path(out,'S5.rds'))

S5$sample<-'S5'

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
