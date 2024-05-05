library(Seurat)
library(data.table)
library(dplyr)
setwd('/projectnb/tcwlab/LabMember/akg/APOE_Jorganoid_project/')

organoid <- readRDS('outputs/03-CellTypeAnnotation/organoid_post_harmony.rds')

DimPlot(organoid, reduction = "humap", group.by = "ident", label = TRUE)

# organoid is obj, removing ident = 24 (deemed a doublet)
cells_to_remove <- WhichCells(organoid, idents = "24")

# create a new Seurat object w/o those cells
organoid <- subset(organoid, cells = setdiff(Cells(organoid), cells_to_remove))



#Idents:
c("GABAergic", "GABAergic", "Glutamatergic", "Astrocyte", "Astrocyte", "Glutamatergic", "OPC/Oligo", "GABAergic", "GABAergic", "Glutamatergic",
  "NPCs - cycling", "NPCs - cycling", "Astrocyte", "Astrocyte", "Astrocyte", "Astrocyte", "Glutamatergic", "Astrocyte", "NPCs - cycling", "VLMC",
  20, "NPCs - cycling", "Glutamatergic", "Pigmented Epithelial", 24)
organoid <- RenameIdents(organoid, "19" = "VLMC")



##### Cluster 3 ################################################################

markers3 <- read.csv('outputs/05-FindMarkers/Astrocyte_Markers_Cluster_3.csv')
markers3 <- rename(markers3, gene = X)
markers3 <- markers3[markers3$avg_log2FC > 0.25 & markers3$p_val < 0.001, ]


Astrocyte1_markers = c("VAMP5", "ISG15", "PTN", "LGALS1", "PTX3", "IFI44L",
                       "FOS", "CAV1", "SMOC1", "IFIT1", "COL4A1", "IFI6", "EN2")
ast_markers1 <- markers3[markers3$X %in% Astrocyte1_markers ,]

Astrocyte2_markers = c("CCND1", "C1QL1", "SPP1", "FXYD5", "RPS4Y1", "FXYD7",
                       "COTL1", "CTGF", "CYP26B1", "PURPL", "CDKN1A", "LITAF", "OCIAD2")
ast_markers2 <- markers3[markers3$X %in% Astrocyte2_markers ,]

Astrocyte3_markers = c("SPP1", "VGF", "TMEM158", "KCNE4", "ITGA7", "EDNRB", "TPM2",
                       "TFPI2", "SCG2", "TMSB4X", "LGALS3", "RDH10", "RCAN1", "GATM")
ast_markers3 <- markers3[markers3$X %in% Astrocyte3_markers ,]

Astrocyte4_markers = c("PMP2", "NKX6-2", "DBI", "BCAN", "AC004540.2", "PTPRZ1",
                       "BBOX1", "C1orf61", "PTGDS", "S100B", "TIMP3", "GFAP")
ast_markers4 <- markers3[markers3$gene %in% Astrocyte4_markers ,]

Astrocyte5_markers = c("DLK1", "LY6H", "TESC", "COL1A2", "APOE", "FGFBP2", "RBP4",
                       "FBLN5", "RPRML", "BCKDHB", "SV2C", "PEG10", "IFITM2", "PLS3")
ast_markers5 <- markers3[markers3$gene %in% Astrocyte5_markers ,]

Astrocyte6_markers = c("TFF3", "CARD16", "SPON1", "SULF1", "PI16", "S100A10", "DCN",
                       "CCL2", "GOLIM4", "ITIH5", "NTN1", "SLIT2", "IGFBP4", "SULN2")
ast_markers6 <- markers3[markers3$gene %in% Astrocyte6_markers ,]

# organoid <- RenameIdents(organoid, 'A4 (PMP2/PTPRZ1)' = 'A4 \n (PMP2/PTPRZ1)')
FeaturePlot(organoid, features = c('PMP2', "PTPRZ1"), reduction = "humap", max.cutoff = "q95", label = TRUE, label.size = 3 )
#Idents(organoid)<-'seurat_clusters'

###### Cluster 4 ###############################################################
markers4 <- read.csv('outputs/05-FindMarkers/Astrocyte_Markers_Cluster_4.csv')
markers4 <- rename(markers4, gene = X)
markers4 <- markers3[markers4$avg_log2FC > 0.25 & markers4$p_val < 0.001, ]

Astrocyte1_markers = c("VAMP5", "ISG15", "PTN", "LGALS1", "PTX3", "IFI44L",
                       "FOS", "CAV1", "SMOC1", "IFIT1", "COL4A1", "IFI6", "EN2")
ast_markers1 <- markers4[markers4$X %in% Astrocyte1_markers ,]

Astrocyte2_markers = c("CCND1", "C1QL1", "SPP1", "FXYD5", "RPS4Y1", "FXYD7",
                       "COTL1", "CTGF", "CYP26B1", "PURPL", "CDKN1A", "LITAF", "OCIAD2")
ast_markers2 <- markers4[markers4$X %in% Astrocyte2_markers ,]

Astrocyte3_markers = c("SPP1", "VGF", "TMEM158", "KCNE4", "ITGA7", "EDNRB", "TPM2",
                       "TFPI2", "SCG2", "TMSB4X", "LGALS3", "RDH10", "RCAN1", "GATM")
ast_markers4 <- markers4[markers4$X %in% Astrocyte3_markers ,]

Astrocyte4_markers = c("PMP2", "NKX6-2", "DBI", "BCAN", "AC004540.2", "PTPRZ1",
                       "BBOX1", "C1orf61", "PTGDS", "S100B", "TIMP3", "GFAP")
ast_markers4 <- markers4[markers4$gene %in% Astrocyte4_markers ,]

Astrocyte5_markers = c("DLK1", "LY6H", "TESC", "COL1A2", "APOE", "FGFBP2", "RBP4",
                       "FBLN5", "RPRML", "BCKDHB", "SV2C", "PEG10", "IFITM2", "PLS3")
ast_markers5 <- markers4[markers4$gene %in% Astrocyte5_markers ,]

Astrocyte6_markers = c("TFF3", "CARD16", "SPON1", "SULF1", "PI16", "S100A10", "DCN",
                       "CCL2", "GOLIM4", "ITIH5", "NTN1", "SLIT2", "IGFBP4", "SULN2")
ast_markers6 <- markers4[markers4$gene %in% Astrocyte6_markers ,]

organoid <- RenameIdents(organoid, '4' = 'A6 \n (CCL2)')
FeaturePlot(organoid, features = c('CCL2'), reduction = "humap", max.cutoff = "q95", label = TRUE, label.size = 3 )

####### Cluster 12 #############################################################
markers12 <- read.csv('outputs/05-FindMarkers/Astrocyte_Markers_Cluster_12.csv')
markers12 <- rename(markers12, gene = X)

markers12 <- markers3[markers12$avg_log2FC > 0.25 & markers12$p_val < 0.001, ]

Astrocyte1_markers = c("VAMP5", "ISG15", "PTN", "LGALS1", "PTX3", "IFI44L",
                       "FOS", "CAV1", "SMOC1", "IFIT1", "COL4A1", "IFI6", "EN2")
ast_markers_1 <- markers12[markers12$X %in% Astrocyte1_markers ,]

Astrocyte2_markers = c("CCND1", "C1QL1", "SPP1", "FXYD5", "RPS4Y1", "FXYD7",
                       "COTL1", "CTGF", "CYP26B1", "PURPL", "CDKN1A", "LITAF", "OCIAD2")
ast_markers_2 <- markers12[markers12$X %in% Astrocyte2_markers ,]

Astrocyte3_markers = c("SPP1", "VGF", "TMEM158", "KCNE4", "ITGA7", "EDNRB", "TPM2",
                       "TFPI2", "SCG2", "TMSB4X", "LGALS3", "RDH10", "RCAN1", "GATM")
ast_markers_3 <- markers12[markers12$X %in% Astrocyte3_markers ,]

Astrocyte4_markers = c("PMP2", "NKX6-2", "DBI", "BCAN", "AC004540.2", "PTPRZ1",
                       "BBOX1", "C1orf61", "PTGDS", "S100B", "TIMP3", "GFAP")
ast_markers_4 <- markers12[markers12$gene %in% Astrocyte4_markers ,]

Astrocyte5_markers = c("DLK1", "LY6H", "TESC", "COL1A2", "APOE", "FGFBP2", "RBP4",
                       "FBLN5", "RPRML", "BCKDHB", "SV2C", "PEG10", "IFITM2", "PLS3")
ast_markers_5 <- markers12[markers12$gene %in% Astrocyte5_markers ,]

Astrocyte6_markers = c("TFF3", "CARD16", "SPON1", "SULF1", "PI16", "S100A10", "DCN",
                       "CCL2", "GOLIM4", "ITIH5", "NTN1", "SLIT2", "IGFBP4", "SULN2")
ast_markers_6 <- markers12[markers12$gene %in% Astrocyte6_markers ,]

# organoid <- RenameIdents(organoid, '12' = 'A \n ()') #no matches
##### Cluster 13 ###############################################################

markers13 <- read.csv('outputs/05-FindMarkers/Astrocyte_Markers_Cluster_13.csv')
markers13 <- rename(markers13, gene = X)
markers13 <- markers13[markers13$avg_log2FC > 0.25 & markers13$p_val < 0.001, ]

Astrocyte1_markers = c("VAMP5", "ISG15", "PTN", "LGALS1", "PTX3", "IFI44L",
                       "FOS", "CAV1", "SMOC1", "IFIT1", "COL4A1", "IFI6", "EN2")
ast_markers_1 <- markers13[markers13$X %in% Astrocyte1_markers ,]

Astrocyte2_markers = c("CCND1", "C1QL1", "SPP1", "FXYD5", "RPS4Y1", "FXYD7",
                       "COTL1", "CTGF", "CYP26B1", "PURPL", "CDKN1A", "LITAF", "OCIAD2")
ast_markers_2 <- markers13[markers13$X %in% Astrocyte2_markers ,]

Astrocyte3_markers = c("SPP1", "VGF", "TMEM158", "KCNE4", "ITGA7", "EDNRB", "TPM2",
                       "TFPI2", "SCG2", "TMSB4X", "LGALS3", "RDH10", "RCAN1", "GATM")
ast_markers_3 <- markers13[markers13$X %in% Astrocyte3_markers ,]

Astrocyte4_markers = c("PMP2", "NKX6-2", "DBI", "BCAN", "AC004540.2", "PTPRZ1",
                       "BBOX1", "C1orf61", "PTGDS", "S100B", "TIMP3", "GFAP")
ast_markers_4 <- markers13[markers13$gene %in% Astrocyte4_markers ,]

Astrocyte5_markers = c("DLK1", "LY6H", "TESC", "COL1A2", "APOE", "FGFBP2", "RBP4",
                       "FBLN5", "RPRML", "BCKDHB", "SV2C", "PEG10", "IFITM2", "PLS3")
ast_markers_5 <- markers13[markers13$gene %in% Astrocyte5_markers ,]

Astrocyte6_markers = c("TFF3", "CARD16", "SPON1", "SULF1", "PI16", "S100A10", "DCN",
                       "CCL2", "GOLIM4", "ITIH5", "NTN1", "SLIT2", "IGFBP4", "SULN2")
ast_markers_6 <- markers13[markers13$gene %in% Astrocyte6_markers ,]

organoid <- RenameIdents(organoid, '13' = 'A5 \n (PEG10)') #1match


####### Cluster 14 #############################################################
markers14 <- read.csv('outputs/05-FindMarkers/Astrocyte_Markers_Cluster_14.csv')
markers14 <- rename(markers14, gene = X)
markers14 <- markers14[markers14$avg_log2FC > 0.25 & markers14$p_val <= 0.001, ]


Astrocyte1_markers = c("VAMP5", "ISG15", "PTN", "LGALS1", "PTX3", "IFI44L",
                       "FOS", "CAV1", "SMOC1", "IFIT1", "COL4A1", "IFI6", "EN2")
ast_markers_1 <- markers14[markers14$X %in% Astrocyte1_markers ,]

Astrocyte2_markers = c("CCND1", "C1QL1", "SPP1", "FXYD5", "RPS4Y1", "FXYD7",
                       "COTL1", "CTGF", "CYP26B1", "PURPL", "CDKN1A", "LITAF", "OCIAD2")
ast_markers_2 <- markers14[markers14$X %in% Astrocyte2_markers ,]

Astrocyte3_markers = c("SPP1", "VGF", "TMEM158", "KCNE4", "ITGA7", "EDNRB", "TPM2",
                       "TFPI2", "SCG2", "TMSB4X", "LGALS3", "RDH10", "RCAN1", "GATM")
ast_markers_3 <- markers14[markers14$X %in% Astrocyte3_markers ,]

Astrocyte4_markers = c("PMP2", "NKX6-2", "DBI", "BCAN", "AC004540.2", "PTPRZ1",
                       "BBOX1", "C1orf61", "PTGDS", "S100B", "TIMP3", "GFAP")
ast_markers_4 <- markers14[markers14$gene %in% Astrocyte4_markers ,]

Astrocyte5_markers = c("DLK1", "LY6H", "TESC", "COL1A2", "APOE", "FGFBP2", "RBP4",
                       "FBLN5", "RPRML", "BCKDHB", "SV2C", "PEG10", "IFITM2", "PLS3")
ast_markers_5 <- markers14[markers14$gene %in% Astrocyte5_markers ,]

Astrocyte6_markers = c("TFF3", "CARD16", "SPON1", "SULF1", "PI16", "S100A10", "DCN",
                       "CCL2", "GOLIM4", "ITIH5", "NTN1", "SLIT2", "IGFBP4", "SULN2")
ast_markers_6 <- markers14[markers14$gene %in% Astrocyte6_markers ,]

organoid <- RenameIdents(organoid, '14' = 'A5 \n (PEG10)') #1match


######### Cluster 17 ###########################################################

## Definitely Astrocyte (some specific labels like S100B and AQP4)
markers17 <- read.csv('outputs/05-FindMarkers/Astrocyte_Markers_Cluster_17.csv')
markers17 <- rename(markers17, gene = X)
markers17 <- markers17[markers17$avg_log2FC > 0.25 & markers17$p_val_adj <= 0.001, ]

Astrocyte1_markers = c("VAMP5", "ISG15", "PTN", "LGALS1", "PTX3", "IFI44L",
                       "FOS", "CAV1", "SMOC1", "IFIT1", "COL4A1", "IFI6", "EN2")
ast_markers_1 <- markers17[markers17$X %in% Astrocyte1_markers ,]

Astrocyte2_markers = c("CCND1", "C1QL1", "SPP1", "FXYD5", "RPS4Y1", "FXYD7",
                       "COTL1", "CTGF", "CYP26B1", "PURPL", "CDKN1A", "LITAF", "OCIAD2")
ast_markers_2 <- markers17[markers17$X %in% Astrocyte2_markers ,]

Astrocyte3_markers = c("SPP1", "VGF", "TMEM158", "KCNE4", "ITGA7", "EDNRB", "TPM2",
                       "TFPI2", "SCG2", "TMSB4X", "LGALS3", "RDH10", "RCAN1", "GATM")
ast_markers_3 <- markers17[markers17$X %in% Astrocyte3_markers ,]

Astrocyte4_markers = c("PMP2", "NKX6-2", "DBI", "BCAN", "AC004540.2", "PTPRZ1",
                       "BBOX1", "C1orf61", "PTGDS", "S100B", "TIMP3", "GFAP")
ast_markers_4 <- markers17[markers17$gene %in% Astrocyte4_markers ,]

Astrocyte5_markers = c("DLK1", "LY6H", "TESC", "COL1A2", "APOE", "FGFBP2", "RBP4",
                       "FBLN5", "RPRML", "BCKDHB", "SV2C", "PEG10", "IFITM2", "PLS3")
ast_markers_5 <- markers17[markers17$gene %in% Astrocyte5_markers ,]

Astrocyte6_markers = c("TFF3", "CARD16", "SPON1", "SULF1", "PI16", "S100A10", "DCN",
                       "CCL2", "GOLIM4", "ITIH5", "NTN1", "SLIT2", "IGFBP4", "SULN2")
ast_markers_6 <- markers17[markers17$gene %in% Astrocyte6_markers ,]

organoid <- RenameIdents(organoid, = ) # no match

########### Cluster 20 #########################################################
# Apoe33Ch 20 vs. other
markers20 <- read.csv('outputs/05-FindMarkers/APOE33Ch_Markers_Cluster_20.csv')
markers20 <- rename(markers20, gene = X)
markers20 <- markers20[markers20$avg_log2FC > 0.25 & markers20$p_val <= 0.001, ]




############# Cluster 22 #######################################################
markers22 <- read.csv('outputs/05-FindMarkers/Markers_Cluster_22.csv')
markers22 <- rename(markers22, gene = X)
markers22 <- markers23[markers22$avg_log2FC > 0 & markers22$p_val <= 0.05, ]




############## Cluster 23 ######################################################
# epithelial (check slack)  pigmented epithelium  
# 
markers23 <- read.csv('outputs/05-FindMarkers/Markers_Cluster_23.csv')
markers23 <- rename(markers23, gene = X)
markers23 <- markers23[markers23$avg_log2FC > 0.25 & markers23$p_val <= 0.05, ]

indiv_markers <- c("ITGA8", "CLIC6")

microglia_markers23 <- markers23[markers23$gene %in% indiv_markers ,]

microglia_cell_markers <- c("HEXB", "SALL1", "GPR34", "TMEM119", "P2RY12",
                            "OLFML3", "TGFBR1", "MERTK", "CD68", "CD206", "CD45")
microglia_markers23 <- markers23[markers23$gene %in% microglia_cell_markers ,]

oligo_markers <- c("MBP", "MAG")
oligo_markers23 <- markers23[markers23$gene %in% oligo_markers ,]


# Exc Neuron - Serotonin Receptor
## General Neuron Markers

Neuronal_markers <- c("DCX", "INSM1", "ST18", "NHLH1", "NEUROD1", "NEUROD6",
                      "PRDM8", "MYT1L", "MAPT", "SNAP25", "SYP", "DLX6")
### Cluster 20: ANK3. 






################ Cluster 24 ####################################################
# not sure if doublet or not 
markers24 <- read.csv('outputs/05-FindMarkers/Markers_Cluster_24.csv')
markers24 <- rename(markers24, gene = X)
markers24 <- markers24[markers24$avg_log2FC > 0.25 & markers24$p_val <= 0.05, ]

dendritic_cell_markers <- c("LYZ", "RGS1", "CD74", "DQA1", "CD11C", "LIN", "CD1A")
dentritic_cells24 <- markers24[markers24$gene %in% dendritic_cell_markers, ] 
# ^ No matches.
### determined to be doublet btwn Gluta and Oligo. Remove this.


################## Examining SingleR output ####################################

## What I need to do (03/15 - 03/17):
# Run FindMarkers() on cluster 20, but only with APOE33Ch data
# Examine the data using Deepti's newer markers
# 
# We'll look at individual unique markers for each of the Astrocyte clusters (for example)
# and from there annotate them as 'Astrocyte - <unique_marker>' <- must be unique to this
# subtype of astrocyte. Use FindMarker() specifically to only compare Astrocyte clusters
# to each other. Once you do this, you can analyze whether this is indeed unique
# to this cluster.

## Examine code to make sure you used DoubletFinder.
## Probably not a doublet.

# organoid$nFeature_RNA
# VlnPlot(organoid, group.by = "seurat_clusters", features = "nCount_RNA", pt.size = 0)


#### 03/22 ####
# 1. look to see if your ast markers (compared to each other) is similar to what Deepti has in her MCC data.
# 2. is the cluster we find in our organoids similar to the one from the paper's organoid data?
# -> check if for each astrocyte if we can also identify this astro subtype in Apoe33Ch paper data.
# Is APOE33Ch-specific cluster also present in the paper. 
# 2-D vs. 3-D models, so we expect a difference. Annotate based on 


### For cluster 20, check if there is an overlap between SEAAD neurons of inhibitory
# (specifically, look into motor neurons, if possible, likely not) 
# Cerebellum scRNA-Seq data. Look to see if you can find a match of these gene markers
# from cluster 20 to one of the clusters from Cerebellum data.

## Before this, you'll look through the APOE3Ch paper, 


## Try adjusting the n.neighbors and other parameters of the UMAP to see if you'll
# find an interaction between two clusters. This could be revealing of some of the 
# nature how these are related to each other.


