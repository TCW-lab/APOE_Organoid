### 10-Organoid-Generate-Pseudobulk (last edit: 04-12-2024)

library(Seurat)
options(Seurat.object.assay.version = 'v3')
library(data.table)
library(tibble)
library(Matrix)

setwd('/projectnb/tcwlab/LabMember/akg/APOE_Jorganoid_project')
dir <- '/projectnb/tcwlab/LabMember/akg/APOE_Jorganoid_project/outputs/01-Get-Pseudobulk'

fp<-function(...)file.path(...)

#we create a list containing the original cell type name 
cellTypes <- c('glutamatergic', 'GABAergic', 'Astrocyte', 'OPC', 
               'Epithelial', 'NPCs', 'VLMC', 'unknown')
cellType <- 'VLMC'
# dir w cellTypes .rds files
inputDir <- './outputs/01-Get-Pseudobulk/Seurat/'
outputDir <- './outputs/01-Get-Pseudobulk'

# seuratObject <- readRDS('./outputs/01-Get-Pseudobulk/Seurat/Astrocyte.rds')

# for each cell type
for(cellType in cellTypes) {
  # Seurat file for the current cell type
  seuratPath <- paste0(inputDir, cellType, '.rds')
  
  # read Seurat file
  seuratObject <- readRDS(seuratPath)
  
  genesDet <- rowSums(GetAssayData(seuratObject, assay = "SoupX_counts"))
  UMI_counts <- Matrix::colSums(GetAssayData(seuratObject, assay = "SoupX_counts"))
  
  
  seuratObject$median_UMIs <- UMI_counts
  seuratObject@meta.data$median_genes_det <- genesDet[rownames(seuratObject@meta.data)]
  
  
  pseudo_mat<-AggregateExpression(seuratObject, assays = 'SoupX_counts', 
                                  group.by = c('sample'), return.seurat = FALSE)
  
  #output aggregated expression into .csv
  output <- fp(dir, paste0(cellType, '_pseudobulk.csv.gz'))
  
  # get rid of cols with zero elements
  df <- as.data.frame(pseudo_mat$SoupX_counts)
  df <- df[, colSums(is.na(df)) < nrow(df)]
  
  fwrite(data.table(df, keep.rownames = 'gene_id'), file = output)
  

  
  
  # calc metrics for meta.data (mtd)
  mtd <- seuratObject@meta.data
  setDT(mtd, keep.rownames = TRUE)
  
  mtd[, n.cells := .N, by = .(sample)]
  #update number of cells per sample / genotype
  mtd[, n.cells.genotype:= .N, by=.(genotype)] #important
  
  mtd[, genotype_total := .N, by = genotype]
  # #calc proportion of cells of specific genotype / total cells per sample
  mtd[, prop.cells := n.cells / genotype_total, by = .(genotype)]
  # 
  # mtd[, median_UMIs := as.numeric(median_UMIs)]
  # mtd[, median_genes_det := as.numeric(median_genes_det)]
  # 
  # #calc median UMIs and genes detected per cell, per sample, per genotype
  # mtd[, med.umis.per.cell := median(median_UMIs, na.rm = TRUE), by = .(sample, genotype)]
  # mtd[, med.genes.per.cell := median(median_genes_det, na.rm = TRUE), by = .(sample, genotype)]
  
  mtd[,avg.pct.mt.per.cell:=mean(`percent.mt`),by=c('sample')] #important
  #flag donors with not enough cells
  # unique(mtd[,pass.threshold.n.cells:=n.cells>20])
  mtd[, pass.threshold.n.cells := n.cells > 20]
  mtd[, outlier.n.cells := !pass.threshold.n.cells, by = .(sample)]

  mtsc<-unique(mtd,by=c('sample', 'genotype'))
  # mtscf<-RemoveUselessColumns(mtsc,key_cols=c('sample','genotype'),pattern_to_exclude = 'ATAC|Multiome|Doublet|Number.of|Genes.') # not sure how to get this to work
  mtsc <- mtsc[,c("sample", "genotype", "median_UMIs", "n.cells", "pass.threshold.n.cells", "avg.pct.mt.per.cell", "prop.cells")] # subset to keep certain columns

  fwrite(mtsc, file.path(dir, paste0(cellType, '_sample_level_metadata.csv.gz')), row.names = TRUE)
  
}


