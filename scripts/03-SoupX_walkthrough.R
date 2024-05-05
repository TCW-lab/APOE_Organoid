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

getGeneSymbols <- function(ensembl_ids) {
  require(org.Hs.eg.db)
  gene_symbols <- mapIds(org.Hs.eg.db,
                         keys = ensembl_ids,
                         column = "SYMBOL",
                         keytype = "ENSEMBL",
                         multiVals = "first")
  return(gene_symbols)
}
setwd('/projectnb/tcwlab/LabMember/akg/APOE_Jorganoid_project/')

out<-'outputs/simpleafSeurat'
source <- '/projectnb/tcwlab/LabMember/akg/APOE_Jorganoid_project/outputs/SoupX'


custom_format <- list("counts" = c("U","S","A"))

sce <- fishpond::loadFry("outputs/1_S16/af_quant/",
                         outputFormat = custom_format)

# Replace Ensembl IDs in row names with gene symbols
ensembl_ids_sce <- rownames(sce)
gene_symbols_sce <- getGeneSymbols(ensembl_ids_sce)

# Replace Ensembl IDs with gene symbols in row names
# If no gene symbol is found, keep original Ensembl ID
updated_rownames <- ifelse(is.na(gene_symbols_sce[ensembl_ids_sce]), 
                           ensembl_ids_sce, gene_symbols_sce[ensembl_ids_sce])
rownames(sce) <- updated_rownames

#create the seurat object by filtering for selected cells
S1Counts <- counts(sce)
S1Counts

library(SoupX)
## This is just for one object. I will have to split the Seurat object by sample ID
## From there, run each split object through the SoupX code (use lapply to do this)
## lapply is just a for loop that has input of a list, you will use splitObject of Seurat
## provide it with a list of Seurat objects. lapply wants a function
## return it with a 
# seurat_list <- lapply(seurat_object, function(){
#   return(seurat_object)
#}
#)

#raw background

library(Seurat)

Organoid_list<-SplitObject(Organoid,split.by = 'SampleID')

Organoid_list<-lapply(Organoid_list, function(seurat_obj){
  #SOUpX code
  
  
  seurat_obj[['SoupX_count']]<-soupxCount
  
  return(seurat_obj)
})


Organoid<-merge(Organoid_list[[1]],Organoid_list[1])



raw <- S1Counts

dim(S1Counts)
dim(S1)
S1 <- readRDS(file = "./outputs/simpleafSeurat/sample1.rds")
genes <- intersect(rownames(S1), rownames(raw))


setdiff(rownames(S1), genes)
#modify genes in raw


raw <- raw[genes,]
dim(raw)

#run SoupX algo (this package is what I ended up using instead of DecontX)
sc <- SoupChannel(raw, S1@assays$RNA@counts[genes,])

# normalize, findclusters, findneighbors, ETC.

sc <- setClusters(sc, S1$initial_seurat_clustering_0.5)
sc <- autoEstCont(sc, doPlot = FALSE)
soup_out <- adjustCounts(sc, roundToInt = TRUE) #adjusted counts matrix

TCW[["SoupX_counts"]] <- CreateAssayObject(counts = soup_out)
TCW@assays$RNA@counts <- soup_out #USED

#calculate percent of mRNA removed
percent_removed_soup <- (1 - (sum(TCW@assays$RNA@counts) / sum(TCW@assays$original_counts@counts)))*100












library(Seurat)

Organoid_list<-SplitObject(Organoid,split.by = 'SampleID')

Organoid_list<-lapply(Organoid_list, function(seurat_obj){
  #SOUpX code
  
  
  seurat_obj[['SoupX_count']]<-soupxCount
  
  return(seurat_obj)
})


Organoid<-merge(Organoid_list[[1]],Organoid_list[1])

