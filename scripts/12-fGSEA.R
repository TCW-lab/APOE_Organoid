##### Script adapted from Alexandre's fGSEA script   04-16-2024  ###############
##### That script is available here: 
# https://github.com/TCW-lab/singlecell-brain/blob/main/notebooks/fgsea_analysis.ipynb 


###### Here we will describe how to perform fgsea analysis on MSIGDB data ######
# MSIGDB (for molecular signature database) gene sets (in gmt format) have been 
# downloaded from Msigdb website, and saved in projectnb/tcwlab/MSigDB.
# We have downloaded and analyzed only the CP: Canonical pathways, containing 
# notably the curated KEGG and REACTOME pathways, and C5: ontology gene sets, 
# containing all Gene ontology (GO) terms.


setwd('/projectnb/tcwlab/LabMember/akg/APOE_Jorganoid_project/')

out<-'outputs/11-fGSEA/'
dir.create(out)
source('../utils/pca_utils.R')
source('../utils/emmaplot.R')
source('../utils/r_utils.R')
library(fgsea)


##########  Prepare the differential expression results ########################
# we need to check that the genes are in gene symbol format (gene names), and that 
# there is no duplicated gene names.
# We also need here to define/calculate the gene statistics/score that we will used. 
# Here the score used is the statistic of the deseq2 results already save in the 
# 'stat' column, so we do not need to do supplementatl calculation. For others case 
# like differential expression results from FindMarkers() Seurat function, such 
# stat column is not present, we have to use another score derived from the differential 
# analysis results. Most of the time we are using sign(avg_log2FC)*-log10(p_value) as score.
# 
# NOTE: fgsea can also be used for any analysis generating a score/statistic by 
# gene Here we used the astrocyte differential expression generated previously



#set the stat column
stat_column='stat'

#Run next lines IF differential results from FindMarkers() results or DE results without 'stat'/'beta' column
#res_de[,diff_score:=sign(avg_log2FC)-log10(p_value)]
#stat_column='diff_score'



#load data
file_paths <- c('./outputs/11-Pseudobulk-DESeq2/res_pseudobulkDESeq2_APOE33Ch_vs_33_glutamatergic_Cov_avg.pct.mt.per.cell_median.UMIs_n.cells.csv.gz',
                             './outputs/11-Pseudobulk-DESeq2/res_pseudobulkDESeq2_APOE22_vs_33_glutamatergic_Cov_avg.pct.mt.per.cell_median.UMIs_n.cells.csv.gz',
                             './outputs/11-Pseudobulk-DESeq2/res_pseudobulkDESeq2_APOE44_vs_33_glutamatergic_Cov_avg.pct.mt.per.cell_median.UMIs_n.cells.csv.gz')
#Prep data
#duplicate gene names?
### Need to add a genotype column for each (should've done this prev step)
genotypes <- c('APOE33Ch', 'APOE22', 'APOE44')

res_de_list <- lapply(seq_along(file_paths), function(i) {
  df <- fread(file_paths[i])
  df[, genotype := genotypes[i]]
  return(df)
})

res_de <- rbindlist(res_de_list)

res_de[,is.dup:=duplicated(gene),by='genotype']
res_de[,n.dup:=sum(is.dup),by=.(genotype,gene)]
table(res_de$n.dup)

# if duplicates, remove the duplicate with the smallest score/statistic
#res_de[order(brain_region,-abs(stat)),to_rm:=duplicated(gene),by='brain_region']
#res_de<-res_de[!(to_rm)]



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

res_gsea_all<-Reduce(rbind,lapply(unique(res_de$genotype), function(CG){
  
  message(paste('testing enrichement in',CG))
  
  res_def <- res_de[genotype == CG, ]
  #calculate/extract the gene stats
  gene_stats <- setNames(res_def[[stat_column]],res_def$gene)
  
  #run fgsea for every pathway source
  res_gsea<-Reduce(rbind,lapply(pathways_to_test, function(p){
    
    message(paste('testing enrichement for',gmt_mtd[name==p]$desc))
    
    pathways<- gmtPathways(file.path(pathways_dir,gmt_mtd[name==p]$gmt))
    
    res<-fgsea(pathways,
               stats=gene_stats,minSize=10,maxSize=2000,scoreType='std')
    
    return(res[,source:=p])
    
  }))
  return(res_gsea[,genotype:=CG])
}))





######            We add information about the pathways and save          ######


#res_gsea_all<-merge(res_gsea_all,pathways_info,by=c('pathway'))[order(brain_region,source,subcat,pval)]
table(res_gsea_all[padj<0.05]$genotype)

fwrite(res_gsea_all,fp(out,"res_fgsea_APOE2n3Chn4vs3_glutamatergic.csv.gz"))


# Pre build function
# 
# This analysis can be run simply and more efficiently using the prebuild functions
# available in the alexandre-utils/deseq2_utils.R functions
source('../utils/deseq2_utils.R')
res_gsea_all<-RunFgseaMsigdb(res_de = res_de,group.by = 'genotype',
                             score = stat_column)

head(res_gsea_all)



##################### See the results using emmaplot ###########################

# Now we can check the number of pathways enriched (padj<0.05) and see the emmaplot

table(res_gsea_all[padj<0.05]$query)

# And see the emmaplot, for example the top100 of MTG brain region

source('../utils/emmaplot.R')

res_gsea_sig<-res_gsea_all[padj<0.05]
res_gsea_top100_APOE33Ch<-head(res_gsea_sig[query=='APOE33Ch'][order(padj)],50)

emmaplot(res_gsea_top100_APOE33Ch)

