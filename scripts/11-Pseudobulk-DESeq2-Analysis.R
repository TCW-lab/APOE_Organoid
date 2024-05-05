### This script is taken from Alexandre Pelletier and changed to suit my #######
### Organoid data needs. 04-12-2024 ############################################
setwd('/projectnb/tcwlab/LabMember/akg/APOE_Jorganoid_project/')
out<-'outputs/11-Pseudobulk-DESeq2'
dir.create(out)
source('../utils/emmaplot.R') # Clone emmaplot.utils from GitHub and place it in your utils folder
source('../utils/pca_utils.R')
library(edgeR)
library(DESeq2)
library(sctransform)
library(data.table)
library(limma)
library(ggplot2)
fp<-function(...)file.path(...)

### in alphabetical order: Astrocyte, Epithelial, Astrocyte, glutamatergic
### NPCs, OPC, unknown, Epithelial
# load data 
pseudo_count<-fread('/projectnb/tcwlab/LabMember/akg/APOE_Jorganoid_project/outputs/01-Get-Pseudobulk/Epithelial_pseudobulk.csv.gz') # Change this
head(pseudo_count)

#transform to matrix
mat<-as.matrix(data.frame(pseudo_count,row.names = 'gene_id'))


# read in the metadata for Epithelials (e.g.)
mtd<-fread('/projectnb/tcwlab/LabMember/akg/APOE_Jorganoid_project/outputs/01-Get-Pseudobulk/Epithelial_sample_level_metadata.csv.gz')
#head(mtd)


# I) QC of samples and genes
# 
# We remove samples if flagged as outliers in previous QC / exploratory analysis (See SEAAD_QC.ipynb ).
# Here, all donors with clinical or cellular abnormalities are excluded, as well as pseudobulk samples derived from less than 50 cell (flag outlier.n.cells)
# 
#remove samples if flagged as outliers (See QC script) 
cat('outliers stats:')
mtdf<-mtd[(pass.threshold.n.cells)]
nrow(mtd)
nrow(mtdf)

to_keep <- mtdf$sample
cat(length(to_keep),' samples kept after QC\n')
mat <- mat[, to_keep]

#genes significantly express (thr= 1CPM) in less than 10% of sample are removed
isexpr <- rowSums(cpm(mat)>1) >= 0.1 * ncol(mat)
cat(sum(isexpr),' genes kept after QC\n')

matf <- mat[isexpr,]

# II) Influence of covariates
# Association with our factor of interest

## we have to factor correctly the categorical covariate (put as first position the 'reference')
# levels correctly categorical factor

mtdf[,genotype:=factor(genotype,levels = c('APOE33','APOE22', 'APOE33Ch', 'APOE44'))]
mtdf[,genotype:=factor(genotype,levels = c('APOE33', 'APOE33Ch', 'APOE44'))]


# Influence of the covariates on the gene expression
# 
# We then perform a PCA to assess influence of each covariates in the transcriptome variance. 
# These covariates contributing to the main sources of variance in our data have to be 
# integrated in our model to take into account these important factors of variability.
# 
# To study that, we first normalized and scale the data using a variance stabilizing 
# transformation (vst) designed for UMI count and single cell data. see this this article.
# We then perform PCA on this UMIs depth normalized matrix.

dim(matf)
covs_to_check<-list(categorical=c('genotype'),
                    nums=c('avg.pct.mt.per.cell', 'median_UMIs', 'n.cells'))






head(mtdf)
#normalized the data using a variance stabilizing transformation (vst)
matf_norm_scaled <- sctransform::vst(matf)$y #y is pearson residual, so already scaled data
# https://ucdavis-bioinformatics-training.github.io/2018-June-RNA-Seq-Workshop/thursday/DE.html
# Doing this for just 'unknown' cluster, as sctransform::vst is not working for this
# 
# library(edgeR)
# counts <- round(matf) 
# dge <- DGEList(counts=counts)
# 
# 
# keep <- filterByExpr(dge)
# dge <- dge[keep,, keep.lib.sizes=FALSE]
# 
# 
# # only genotype b/c lack of degrees of freedom
# design <- model.matrix(~ genotype +  avg.pct.mt.per.cell + median_UMIs + n.cells, data=mtdf)
# v <- voom(dge, design, plot=TRUE)
# 
# 
# fit <- lmFit(v, design)
# fit <- eBayes(fit)
# 
# res <- topTable(fit, adjust="BH", number=Inf)
names(res)[3] <- 'log2FC'
names(res)[4] <- 'stat'
# fwrite(res, "/projectnb/tcwlab/LabMember/akg/APOE_Jorganoid_project/outputs/11-Pseudobulk-DESeq2/res_pseudobulkDESeq2_APOE33Ch_44_vs_33_unknown_Cov_avg.pct.mt.per.cell_median.UMIs_n.cells.csv.gz")
# 
# fwrite(res,"/projectnb/tcwlab/LabMember/akg/APOE_Jorganoid_project/outputs/11-Pseudobulk-DESeq2/res_pseudobulkDESeq2_APOE22_vs_33_Epithelial_Cov_avg.pct.mt.per.cell_median.UMIs_n.cells.csv.gz"))



pca<-RunPca(matf_norm_scaled, scale = F,center = F)









#Let's check the influence of each covariates in the 2 first PCs
i <- 1
covs_to_plot<-head(unlist(covs_to_check)[i:length(unlist(covs_to_check))],4)
# covs_to_plot
print(PcaPlot(pca,mtdf,
              group.by = 'genotype',
              sample_col = 'sample',ncol=2))

for(i in seq.int(from = 1, to = length(unlist(covs_to_check)), by = 4)) {
  covs_to_plot <- head(unlist(covs_to_check)[i:length(unlist(covs_to_check))], 4)
  print(PcaPlot(pca, mtdf,
                group.by = covs_to_plot,
                sample_col = 'sample', ncol = 2),
                group.by = 'genotype')
  ggsave(fp(out, paste0('pca_Epithelial_pseudobulk_covs_', paste(covs_to_plot, collapse = '_'), '.png')), width = 7, height = 6)
  
}



# 
# res_cor_pcs<-CorrelCovarPCs(pca, mtdf,
#                             sample_col = 'sample',
#                             vars_num = covs_to_check$nums,
#                             vars_cat = covs_to_check$categorical)
# 
# plotPvalsHeatMap(res_cor_pcs,p.thr = 0.1,p_col='padj',
#                  labels_col =paste0(paste0('PC',1:10),"(",round(pctPC(pca,rngPCs = 1:10)*100,0),"%)") )



#ggsave(fp(out,ps('pca_VLMC_pseudobulk_covs_',paste(covs_to_plot,collapse = '_'),'.png')),width = 7,height = 6)



##################        III) run DESEQ2        ###############################

# Before performing the differential expression analysis, we scale the numerical 
# covariates that we decided to include to our model. It helps GLM convergence

#scale numerical covariates to improve GLM convergence
covs_to_scale<-c('median_UMIs','avg.pct.mt.per.cell','n.cells') #median_UMIs alrdy scaled, might not need to # CHECK PCA plot for these.
mtdf_scaled<-copy(mtdf)
mtdf_scaled[,(covs_to_scale):=lapply(.SD,scale),.SDcols=covs_to_scale]

## Create your linear reg. design based on w/e covariates are deemed important
design= ~ genotype + avg.pct.mt.per.cell + median_UMIs + n.cells


# We then create the DESEQ2 object and Run DESEQ2 model fitting.
# DESeq2 uses the raw count matrix to estimate the negative binomial model parameters. 
# it requires thus the raw/unnormalized (filtered) matrix.
# We also include the DESEq2 object the metadata and the chosen design of our model


#run DESEQ2
#running this on pseudobulk matrix (matf), design as design (genotype, potentially mt.pct and median umis), 
dds <- DESeqDataSetFromMatrix(matf, 
                              colData = data.frame(mtdf_scaled,row.names="sample")[colnames(matf),], 
                              design = design)
dds <- DESeq(dds)


# Next is to extract from this model the contrast of interest. A good way to do so 
# is to use the model.matrix() function. we can explicitly extract the case and control 
# condition to do the contrast. It is also very helpful when studying interaction 
# between 2 covariates See this tutorial for more information.

#extract contrast of interest
mod_mat <- model.matrix(design(dds), colData(dds))

APOE33 <- colMeans(mod_mat[dds$genotype== "APOE33", ])
APOE33Ch <- colMeans(mod_mat[dds$genotype == "APOE33Ch", ])
APOE22 <- colMeans(mod_mat[dds$genotype == "APOE22", ])
APOE44 <- colMeans(mod_mat[dds$genotype == "APOE44", ])

### Sequentially go through this one-by-one (APOE3 vs. APOE4, APOE3 vs. APOE3Ch, etc.) (go through all possible combinations)
res <- results(dds,contrast = APOE33Ch-APOE33,alpha = 0.05)

res<-data.table(as.data.frame(res),keep.rownames="gene")

# We can then check number of DEGs between APOE4 vs 3 after multiple test correction
res[padj<0.25][order(padj)]

# We save some useful information about this analysis directly in the differential
# Expression results dataframe, before save the results

#add some annotation to the results
res[,cell_type:="VLMC"]
res[,design:=paste0('~',as.character(design)[2])]
# fwrite for each comparison (obviously change the filename! with covariates you are using)
fwrite(res,fp(out,"res_pseudobulkDESeq2_APOE22_vs_33_VLMC_Cov_avg.pct.mt.per.cell_median.UMIs_n.cells.csv.gz"))





### Sequentially go through this one-by-one (APOE3 vs. APOE4, APOE3 vs. APOE3Ch, etc.) (go through all possible combinations)
res <- results(dds,contrast = APOE33Ch-APOE33,alpha = 0.05)

res<-data.table(as.data.frame(res),keep.rownames="gene")

# We can then check number of DEGs between APOE4 vs 3 after multiple test correction
res[padj<0.25][order(padj)]

# We save some useful information about this analysis directly in the differential
# Expression results dataframe, before save the results

#add some annotation to the results
res[,cell_type:="VLMC"]
res[,design:=paste0('~',as.character(design)[2])]
# fwrite for each comparison (obviously change the filename! with covariates you are using)
fwrite(res,fp(out,"res_pseudobulkDESeq2_APOE33Ch_vs_33_VLMC_Cov_avg.pct.mt.per.cell_median.UMIs_n.cells.csv.gz"))




### Sequentially go through this one-by-one (APOE3 vs. APOE4, APOE3 vs. APOE3Ch, etc.) (go through all possible combinations)
res <- results(dds,contrast = APOE44-APOE33,alpha = 0.05)

res<-data.table(as.data.frame(res),keep.rownames="gene")

# We can then check number of DEGs between APOE4 vs 3 after multiple test correction
res[padj<0.25][order(padj)]

# We save some useful information about this analysis directly in the differential
# Expression results dataframe, before save the results

#add some annotation to the results
res[,cell_type:="VLMC"]
res[,design:=paste0('~',as.character(design)[2])]
# fwrite for each comparison (obviously change the filename! with covariates you are using)
fwrite(res,fp(out,"res_pseudobulkDESeq2_APOE44_vs_33_VLMC_Cov_avg.pct.mt.per.cell_median.UMIs_n.cells.csv.gz"))


###### DEG Analysis ############################################################

# for each main cell type (excluding cluster 20 and others putative too small cluster),
# - barplot of the number of DEGs (padj<0.05 and abs(log2FC) > 0.25in each APOExx vs APOE33 comparison
# 
# For cell type and APOE comparison with more than 10 DEGs
# Emmaplot of the top 50 pathways enriched at padj <0.05
# Also the Barplot of the Proportion of each genotype in each cell type will be useful

######## Instructions for 04/18 ^ ##############################################



file_paths <- c('./outputs/11-Pseudobulk-DESeq2/res_pseudobulkDESeq2_APOE33Ch_vs_33_Astrocyte_Cov_avg.pct.mt.per.cell_median.UMIs_n.cells.csv.gz')
                './outputs/11-Pseudobulk-DESeq2/res_pseudobulkDESeq2_APOE22_44_vs_33_unknown_Cov_avg.pct.mt.per.cell_median.UMIs_n.cells.csv.gz')
                './outputs/11-Pseudobulk-DESeq2/res_pseudobulkDESeq2_APOE44_vs_33_unknown_Cov_avg.pct.mt.per.cell_median.UMIs_n.cells.csv.gz')


#Prep data
#duplicate gene names?
### Need to add a genotype column for each (shouldve done this prev step)
genotypes <- c('APOE33Ch')

res_de_list <- lapply(seq_along(file_paths), function(i) {
  df <- fread(file_paths[i])
  df[, genotype := genotypes[i]]
  return(df)
})

res_de <- rbindlist(res_de_list)

names(res_de)[3] <- "log2FC"

res_de[,is.dup:=duplicated(gene),by='genotype']

res_de[,nDEG_up:= sum(padj < 0.05 & log2FC > 0.25, na.rm = TRUE), by=.(genotype, cell_type)]
res_de[,nDEG_down:= sum(padj < 0.05 & log2FC < -0.25, na.rm = TRUE), by=.(genotype, cell_type)]
table(res_de$n.dup)



# Get unique rows based on genotype and cell_type
res_de_unique <- unique(res_de, by = c('genotype', 'cell_type'))

# melt the data to long format for ggplot
library(reshape2) 
res_de_long <- melt(res_de_unique, id.vars = c("genotype", "cell_type"), measure.vars = c("nDEG_up", "nDEG_down"), 
                    variable.name = "Regulation", value.name = "Count")

# barplot of each genotype vs. APOE33, count(DEGs)
ggplot(res_de_long, aes(x = genotype, y = Count, fill = Regulation)) +
  geom_col(position = position_dodge()) +
  facet_wrap(~cell_type) +
  labs(x = "Genotype", y = "Number of DEGs", title = "Number of Up- and Down-regulated DEGs by Genotype (vs. APOE33)") +
  theme_minimal() + 
  scale_fill_manual(values = c("nDEG_up" = "red4", "nDEG_down" = "navy"))

res_de_long
################# run FGSEA    #################################################
# 
# To run fgsea, you can use either the native way using the gmt file, or the pre
# build function that used the already preprocessed gmt files
# Native way
# 
# specify the path to MSIGDB gmt files and which pathway category to test

pathways_dir <- "/projectnb/tcwlab/MSigDB/"
pathways_to_test<-c("GO_all", 'CP_all') # choose pathway gene sets reference you want to test or let like that to test for all

pathways_info<-fread(file.path(pathways_dir,'all_CPandGOs_genesets_metadata.csv.gz'))
gmt_mtd<-fread(file.path(pathways_dir,'gmt_metadata.csv')) #contain the gmt files names and corresponding gene set category names

unique(res_de$genotype)

stat_column <- 'stat'
# res_gsea_all<-Reduce(rbind,lapply(unique(res_de$genotype), function(CG){
#   
#   message(paste('testing enrichement in',CG))
#   
#   res_def <- res_de[genotype == CG, ]
#   #calculate/extract the gene stats
#   gene_stats <- setNames(res_def[[stat_column]],res_def$gene)
#   
#   #run fgsea for every pathway source
#   res_gsea<-Reduce(rbind,lapply(pathways_to_test, function(p){
#     
#     message(paste('testing enrichement for',gmt_mtd[name==p]$desc))
#     
#     pathways<- gmtPathways(file.path(pathways_dir,gmt_mtd[name==p]$gmt))
#     
#     res<-fgsea(pathways,
#                stats=gene_stats,minSize=10,maxSize=2000,scoreType='std')
#     
#     return(res[,source:=p])
#     
#   }))
#   return(res_gsea[,genotype:=CG])
# }))





######            We add information about the pathways and save          ######


#res_gsea_all<-merge(res_gsea_all,pathways_info,by=c('pathway'))[order(brain_region,source,subcat,pval)]

# res_gsea_all <- fread(input = "outputs/12-fGSEA/res_fgsea_APOE2n3Chn4vs3_glutamatergics.csv.gz")

# Pre build function
# 
# This analysis can be run simply and more efficiently using the prebuild functions
# available in the alexandre-utils/deseq2_utils.R functions
# source('../utils/deseq2_utils.R')
res_gsea_all<-RunFgseaMsigdb(res_de = res_de, group.by = 'genotype',
                             score = 'stat')

head(res_gsea_all)

table(res_gsea_all[padj<0.05]$genotype)

fwrite(res_gsea_all,fp(out,"res_fgsea_APOE2n3Chn4vs3_unknown.csv.gz"))

##################### See the results using emmaplot ###########################

# Now we can check the number of pathways enriched (padj<0.05) and see the emmaplot

table(res_gsea_all[padj<0.05]$query)

# And see the emmaplot, for example the top100 of MTG brain region


# 
# res_gsea_sig<-res_gsea_all[padj<0.05]
# res_gsea_top50_APOE33Ch<-head(res_gsea_sig[genotype=='APOE33Ch'][order(padj)],50)
# res_gsea_top50_APOE22<-head(res_gsea_sig[genotype=='APOE22'][order(padj)],50)
# res_gsea_top50_APOE44<-head(res_gsea_sig[genotype=='APOE44'][order(padj)],50)

### Alternatively, we can select for the top-25 up and down- regulated:
library(dplyr)

res_gsea_sig <- res_gsea_all[padj < 0.05,]

# def a function to get the top pathways for both up- and down-regulated
get_top_pathways <- function(data, n=30) {
  # separate up- and down-regulated pathways
  up_regulated <- data[data$NES > 1][order(padj)][1:n, ]
  down_regulated <- data[data$NES < -1][order(padj)][1:n, ]
  
  # combine top up- and down-regulated pathways
  return(rbind(up_regulated, down_regulated))
}

# Apply function to each genotype
res_gsea_top50_APOE33Ch <- get_top_pathways(res_gsea_sig[query == 'APOE33Ch'], 30)
res_gsea_top50_APOE22 <- get_top_pathways(res_gsea_sig[query == 'APOE22'], 30)
res_gsea_top50_APOE44 <- get_top_pathways(res_gsea_sig[query == 'APOE44'], 30)

res_gsea_top50_APOE44 <- na.omit(res_gsea_top50_APOE44)




emmaplot(res_gsea_top50_APOE44)
emmaplot(res_gsea_top50_APOE33Ch)
emmaplot(res_gsea_top50_APOE22)




