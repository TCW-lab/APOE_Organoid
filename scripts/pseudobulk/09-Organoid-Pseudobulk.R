# 04-09-2024 The following script is taking Alexandre Pelletier's SEAAD_pseudobulk.ipynb
# script and editing it to work with my Organoid data.

# "We create a pseudobulk matrix (aggregates cell counts by sample) for 
# each cell_type, and for each main_cell_type As a minimal working example (MWE), 
# it is perform like that for each cell type:" 
library(Seurat)
options(Seurat.object.assay.version = 'v3')
library(data.table)
library(tibble)
library(Matrix)

setwd('/projectnb/tcwlab/LabMember/akg/APOE_Jorganoid_project/')
dir <- '/projectnb/tcwlab/LabMember/akg/APOE_Jorganoid_project/outputs/01-Get-Pseudobulk/'

fp<-function(...)file.path(...)

?AggregateExpression()

#to match with the rds object name, we create a list containing the original cell type name 
cellTypes <- c('Glutamatergic', 'GABAergic', 'Glutamatergic', 'Astrocyte', 'OPC', 
               'Epithelial', 'NPCs', 'VLMC', 'unknown')

# dir w cellTypes .rds files
inputDir <- './outputs/01-Get-Pseudobulk/Seurat/'
outputDir <- './outputs/01-Get-Pseudobulk/'

seuratObject <- readRDS('./outputs/01-Get-Pseudobulk/Seurat/VLMC.rds')

# for each cell type
for(cellType in cellTypes) {
  # Seurat file for the current cell type
  seuratPath <- fp(dir, paste0(cellType, '.rds'))
  
  # read Seurat file
  seuratObject <- readRDS(seuratPath)
  
  #calculate the median UMIs and median genes detected
  UMI_counts <- Matrix::colSums(GetAssayData(seuratObject, assay = "SoupX_counts"))
  genesDet <- rowSums(GetAssayData(seuratObject, assay = 'SoupX_counts'))
  
  seuratObject$median_UMIs <- UMI_counts
  #seuratObject$median_genes_detected <- genesDet
  
  # aggregate expression per genotype
  pseudo_mat <- AggregateExpression(seuratObject, assays = 'SoupX_counts', 
                                       group.by = c('sample'),
                                       return.seurat = FALSE)
  
  cellType <- 'VLMC'
  output <- fp(dir, paste0(cellType, '_pseudobulk.csv.gz'))
  pseudo_mat$SoupX_counts
  fwrite(data.table(pseudo_mat$SoupX_counts, keep.rownames = "gene_id"), output)
  
  # calc metrics for meta.data (mtd)
  mtd <- seuratObject@meta.data
  
  mtd[,tot.cells.sample:=.N,by=.(`sample`)]
  mtd[,n.cells:=.N,by=.(`sample`,main_cell_type)]
  
  mtd[,prop.cells:=n.cells/tot.cells.donor,by=.(`Donor.ID`,main_cell_type)]
  
  mtd[,med.umis.per.cell:=median(Number.of.UMIs,na.rm = T),by=.(`Donor.ID`,main_cell_type)]
  mtd[,med.genes.per.cell:=median(Genes.detected,na.rm = T),by=.(`Donor.ID`,main_cell_type)]
  
  mtd[,avg.pct.mt.per.cell:=mean(`Fraction.mitochondrial.UMIs`),by=c('sample')]
  
  mtsc<-unique(mtd,by=c('sample'))
  
  
  # output to csv
  output <- fp(baseDir, paste0(cellType, '_pseudobulk.csv.gz'))
  fwrite(data.table(pseudo_matrix$RNA, keep.rownames = 'gene_id'), output)
}


VLMC <- readRDS('./outputs/01-Get-Pseudobulk/Seurat/VLMC.rds')

genesDet <- rowSums(GetAssayData(VLMC, assay = "SoupX_counts"))
UMI_counts <- Matrix::colSums(GetAssayData(VLMC, assay = "SoupX_counts"))

#replace with 'seuratObj'
VLMC$median_UMIs <- UMI_counts
VLMC@meta.data$median_genes_det <- genesDet[rownames(VLMC@meta.data)]


#replace 'VLMC' w seuratObj
pseudo_mat<-AggregateExpression(VLMC, assays = 'SoupX_counts', 
                                group.by = c('sample', 'genotype'), return.seurat = FALSE)

#output aggregated expression into .csv
# change 'VLMC' to 'cellType' for loop
output <- fp(dir, paste0('VLMC', '_pseudobulk.csv.gz'))
fwrite(data.table(pseudo_mat$SoupX_counts, keep.rownames = "gene_id"), output)



# calc metrics for meta.data (mtd)
mtd <- VLMC@meta.data
setDT(mtd)

mtd[, tot.cells.sample := .N, by = .(sample)]
#update number of cells per sample / genotype
mtd[,n.cells:=.N,by=.(sample, genotype)]

#calc proportion of cells of specific genotype / total cells per sample
mtd[,prop.cells:=n.cells/tot.cells.sample,by=.(sample, genotype)]

mtd[, median_UMIs := as.numeric(median_UMIs)]
mtd[, median_genes_det := as.numeric(median_genes_det)]

#calc median UMIs and genes detected per cell, per sample, per genotype
mtd[, med.umis.per.cell := median(median_UMIs, na.rm = TRUE), by = .(sample, genotype)]
mtd[, med.genes.per.cell := median(median_genes_det, na.rm = TRUE), by = .(sample, genotype)]

mtd[,avg.pct.mt.per.cell:=mean(`Fraction.mitochondrial.UMIs`),by=c('sample')]

mtsc<-unique(mtd,by=c('sample','genotype'))

#flag donors with not enough cells
mtd[,pass.threshold.n.cells:=n.cells>50,by=]
mtd[,outlier.n.cells:=!pass.threshold.n.cells,by=c('sample','genotype')]

# mtscf<-RemoveUselessColumns(mtsc,key_cols=c('sample','genotype'),pattern_to_exclude = 'ATAC|Multiome|Doublet|Number.of|Genes.')


fwrite(mtd,fp(dir,'all_final_RNAseq_nuclei_main_cell_type_level_metadata.csv.gz'))



fwrite(data.table(pseudo_mat$SoupX_counts,keep.rownames = 'gene_id'),
       'outputs/01-Get-Pseudobulk/VLMC_pseudobulk.csv.gz')

########### merge the split cell type into one pseudobulk data #################

# "for the cell type divided into different object/rds file(because too large data),
# we aggregate expression into one.
# we thus create a data.table pseudo_files_dt saving these information""

# I am going to run this, but across the 4 unique genotypes (APOE2, APOE3, APOE3Ch, APOE4)

genotype='APOE2'
out1<-fp(out,genotype)
mtd<-fread(fp(out1,'all_final_RNAseq_nuclei_metadata.csv.gz'))

#to match with the rds object name, we create a column containing the original cell type name 
mtd[,original_cell_type_name:=str_remove(str_replace_all(cell_type,' ','_'),'^Inh_|Exc_')]
mtd[,original_cell_type_name:=str_replace_all(original_cell_type_name,'/','-')]
fwrite(mtd,fp(out1,'all_final_RNAseq_nuclei_metadata.csv.gz')) 

mtsc<-unique(mtd,by=c('Donor.ID','cell_type'))

pseudo_files<-list.files(out1,pattern = '\\_pseudobulk\\.csv\\.gz',full.names = T)
pseudo_files_dt<-data.table(file=pseudo_files,
                            original_cell_type_name=str_remove(basename(pseudo_files),'_pseudobulk\\.csv\\.gz'))
pseudo_files_dt[,original_cell_type_name:=str_remove(original_cell_type_name,'\\_set[0-9]+')]


pseudo_files_dt<-merge(pseudo_files_dt,unique(mtsc[,.(original_cell_type_name,cell_type,main_cell_type)]))
pseudo_files_dt[,n_file:=.N,by='cell_type']

for(ct in unique(pseudo_files_dt[n_file>1]$original_cell_type_name)){
  pseudo_list<-lapply(pseudo_files_dt[original_cell_type_name==ct]$file, function(f)fread(f))
  
  #lacking samples column
  samples<-Reduce(union,lapply(pseudo_list,colnames))
  pseudo_list<-lapply(pseudo_list,function(x){
    samples_lacking<-setdiff(samples,colnames(x))
    x[,(samples_lacking):=0]
  })
  
  #transform as matrix
  pseudo_list<-lapply(pseudo_list, function(x)as.matrix(data.frame(x,row.names = 'gene_id')))
  
  #aggregate count by sample / featute
  features<-rownames(pseudo_list[[1]])
  samples<-colnames(pseudo_list[[1]])
  
  pseudo_merge<-Reduce(`+`,lapply(pseudo_list,function(x)x[features,samples]))
  print(head(pseudo_merge[,1:10]))
  fwrite(data.table(pseudo_merge,keep.rownames = 'gene_id'),fp(out1,paste0(ct,'pseudobulk.csv.gz')))
  
}


#we do not need anymore of the pseudobulk data from the divided cell type objects, so we remove them
system(paste('rm',paste(pseudo_files_dt[n_file>1]$file,collapse = ' ')))




