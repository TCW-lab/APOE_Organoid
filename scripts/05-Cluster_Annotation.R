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
# for miQC / flexmix:
library(flexmix)
library(miQC)
library(singleCellTK)
library(fgsea)
library(DropletUtils)
library(SeuratData)
library(tidyr)
library(SeuratWrappers)
library(data.table)
library(SoupX)
library(fishpond)
library(remotes)
library(ggplot2)
library(presto)
setwd("/projectnb/tcwlab/LabMember/akg/projects/APOE_Organoid")

dir <- "/projectnb/tcwlab/LabMember/akg/projects/APOE_Organoid"
out <- '/outputs/05-Cluster_Annotation'


### Step 1. Read in RDS object and plot FeaturePlots for various markers   #####

organoid <- readRDS(paste0(dir, '/outputs/05-Cluster_Annotation', '/organoid.rds'))

DimPlot(organoid, reduction = 'umap', group.by = 'seurat_clusters', label = TRUE)


organoid$harmony_clusters <- Idents(organoid)


Final_Markers_Progenitors = c("TOP2A", "MKI67", "FOXM1", "CENPF", "PAX6", "PCNA", "E2F2", "HOPX", "LHX2", "OTX2", "GLI3", "HES2", "HES6") 
Final_Markers_General_Neuron = c("DCX", "INSM1", "ST18", "NHLH1", "NEUROD1", "NEUROD6", "PRDM8", "MYT1L", "MAPT", "SNAP25", "SYP", "DLX6")
Final_Markers_Glut_Neuron = c("BCL11B", "DLG4", "STMN2", "SLC17A6", "GRIN2B", "RELN", "GRIA2")
Final_Markers_Gaba_Neuron = c("GAD1", "GAD2", "DLX1", "DLX2", "SST", "NRXN1", "ANK3") #rm "DLX6-AS1"
Final_Markers_Astrocyte = c("APOE", "GFAP", "PEA15", "S100B", "ALDH1L1", "FGFR3", "AGT", "AQP4")
Final_Markers_OPC = c("OLIG1", "OLIG2", "PCDH15", "LHFPL3", "MBP") #added MPC, only expr in oligos 
Final_Markers_VLMC = c("LUM", "PDGFRA", "DCN", "POSTN", "OGN", "APOD", "RSPO3") 
Final_Markers_Pericytes = c("CSPG4", "PDGFRB", "VTN", "ACTA2", "KCNJ8", "ABCC9", "ACE2", "ART3", "ATP13A5") #rm CD146
Final_Markers_PigEp = c("CLIC6") # https://www.nature.com/articles/s41467-020-15326-5
Final_Markers_Chloroid = c("ITGA8", "TTR", "FOLR1", "PRLR")
Final_Markers_Radial_Glia = c('PAX6', 'GLAST', 'SOX2', 'MEIS2', 'FOXP2', "TLE4", "VIM")
## Visualize Canonical Cell Type Marker Expression
#Feature Plots

Final_Marker_List = list(Final_Markers_Progenitors, Final_Markers_General_Neuron, Final_Markers_Glut_Neuron, Final_Markers_Gaba_Neuron,
                         Final_Markers_Astrocyte, Final_Markers_OPC, Final_Markers_VLMC, Final_Markers_Pericytes,
                         Final_Markers_PigEp, Final_Markers_Chloroid, Final_Markers_Radial_Glia)
Final_Marker_List_Names = c("Progenitor", "Neuron", "Glutamatergic", "GABA",
                            "Astrocyte", "OPC", "VLMC", "Pericyte",
                            "PigEp", "Chloroid", "Radial_Glia")

count = 1
for (x in Final_Marker_List) {
  for (y in x) {
    
    name = Final_Marker_List_Names[count]
    png(file = paste(dir, "/outputs/images/05-Cluster_Annotation/Feature_Plots/", name, "_", y, "_feature_plot.png", sep = ""), width = 15, height = 15, units = "in", res = 300)
    
    plot <- FeaturePlot(organoid, reduction = "humap", features = y, label = TRUE, label.size = 6) + 
      theme(
        axis.text.x = element_text(size = 24),
        axis.text.y = element_text(size = 24),
        axis.title.x = element_text(size = 24),
        axis.title.y = element_text(size = 24),
        plot.title = element_text(size = 24),
        legend.text = element_text(size = 24),  # Adjust the size of legend text
        legend.title = element_text(size = 24)  # Adjust the size of legend title
      )
    
    print(plot)
    dev.off()
  }
  count = count + 1
}

names(organoid@reductions)
FeaturePlot(organoid, reduction = "humap", label = TRUE, label.size = 6)
# Feature plots should be similar to below:
#Dot Plots 



feats <- list(Final_Markers_Progenitors, Final_Markers_General_Neuron, Final_Markers_Glut_Neuron, Final_Markers_Gaba_Neuron,
              Final_Markers_Astrocyte, Final_Markers_OPC, Final_Markers_VLMC, Final_Markers_Pericytes,
              Final_Markers_PigEp, Final_Markers_Chloroid, Final_Markers_Radial_Glia)
marker_names <- list("Progenitors", "General_Neuron", "Glut", "GABA", "Astrocyte", 
                     "OPC", "VLMC", "Pericyte", "PigEp", "Chloroid", "Radial_Glia")

DotPlot(
  organoid,
  features = c("OLIG1", "OLIG2", "PCDH15", "LHFPL3", "MBP"),
  cols = c("light grey", "blue"),
  col.min = -2.5,
  col.max = 2.5,
  dot.min = 0,
  dot.scale = 10,
  idents = NULL,
  group.by = NULL,
  split.by = NULL,
  cluster.idents = FALSE,
  scale = TRUE,
  scale.by = "radius"
  # scale.min = NA,
  # scale.max = NA
)


count = 0
for (y in feats) {
  
  count = count + 1
  name = marker_names[count]
  
  #dotplot
  png(paste(dir, "/outputs/images/05-Cluster_Annotation/Dot_Plots/", name, "_dot_plot.png", sep = ""), width = 25, height = 15, units = "in", res = 300)
  
  plot <- DotPlot(
    organoid,
    features = y,
    cols = c("light grey", "blue"),
    col.min = -2.5,
    col.max = 2.5,
    dot.min = 0,
    dot.scale = 10,
    idents = NULL,
    group.by = NULL,
    split.by = NULL,
    cluster.idents = FALSE,
    scale = TRUE,
    scale.by = "radius"
    # scale.min = NA,
    # scale.max = NA
  ) + theme(
    axis.text.x = element_text(size = 30),
    axis.text.y = element_text(size = 30),
    axis.title.x = element_text(size = 30),
    axis.title.y = element_text(size = 30),
    legend.text = element_text(size = 30),  # Adjust the size of legend text
    legend.title = element_text(size = 30)
  )
  
  print(plot)
  
  dev.off()
}

## Step 2. Actual Annotation                                                ####

## Astrocyte:
DimPlot(organoid, reduction = "humap", label = TRUE, label.size = 6)

#Idents (For Reference):
c("GABAergic", "Glutamatergic", "NPCs", "GABAergic", "Astrocyte", "Astrocyte", "GABAergic", 7, "Glutamatergic", "Immature Astrocyte",
  "Glutamatergic", "GABAergic", "NPCs (cycling)", 13, 14, 15, "NPCs (cycling)", 17, "Astrocyte", "NPCs (cycling)",
  "VLMC", "NPCs", "Pigmented Epithelial", 23)

### Top level of annotation: Astrocyte
### 4,5, 18, 9, 
### 14 (This is a different sub-level annotation)
### 
organoid <- RenameIdents(organoid, "4" = "Astrocyte")
organoid <- RenameIdents(organoid, "5" = "Astrocyte")
organoid <- RenameIdents(organoid, "18" = "Astrocyte")
organoid <- RenameIdents(organoid, "14" = "Astrocyte") # temporary annotation.
# DimPlot(organoid, reduction = "humap", label = TRUE, label.size = 5)
# 9 is immature Astrocyte
organoid <- RenameIdents(organoid, "9" = "Immature Astrocyte")

## GABAergic:
# 0, 3, 6, 11 (somewhat questionable)
organoid <- RenameIdents(organoid, "0" = "GABAergic")
organoid <- RenameIdents(organoid, "3" = "GABAergic")
organoid <- RenameIdents(organoid, "6" = "GABAergic")
organoid <- RenameIdents(organoid, "11" = "GABAergic")


# DimPlot(organoid, reduction = "humap", label = TRUE, label.size = 5)

## Glutamatergic:
# 1, 8, 10
organoid <- RenameIdents(organoid, "1" = "Glutamatergic")
organoid <- RenameIdents(organoid, "8" = "Glutamatergic")
organoid <- RenameIdents(organoid, "10" = "Glutamatergic")
# DimPlot(organoid, reduction = "humap", label = TRUE, label.size = 5)


## NPCs:
# NPCs (cycling) 12, 16, 19, 13 is dubious
organoid <- RenameIdents(organoid, "12" = "NPCs (cycling)")
organoid <- RenameIdents(organoid, "13" = "NPCs (cycling)")
organoid <- RenameIdents(organoid, "16" = "NPCs (cycling)")
organoid <- RenameIdents(organoid, "19" = "NPCs (cycling)")
# DimPlot(organoid, reduction = 'humap', label = TRUE, label.size = 5)

# NPCs (non-cycling): 
organoid <- RenameIdents(organoid, "2" = "NPCs")
organoid <- RenameIdents(organoid, "21" = "NPCs")
## OPCs:
# 7
organoid <- RenameIdents(organoid, "7" = "OPC")
# DimPlot(organoid, reduction = 'humap', label = TRUE, label.size = 5)

## Pericyte:
# 22
organoid <- RenameIdents(organoid, "22" = "Pericyte")
# DimPlot(organoid, reduction = 'humap', label = TRUE, label.size = 5)
# Pericyte Progenitor
organoid <- RenameIdents(organoid, '15' = 'Pericyte Progenitor Cell')

## VLMC:
# 20
organoid <- RenameIdents(organoid, "20" = "VLMC")
# DimPlot(organoid, reduction = 'humap', label = TRUE, label.size = 4.5)

DimPlot(organoid, reduction = 'humap', label = TRUE)

## Pigmented Epithelial
# 22 for certain.
# organoid <- RenameIdents(organoid, "22" = "Pigmented Epithelial")
## Chloroid Plexus
# Could be 20, but 20 seems more likely to be VLMC

## Unknown cluster (17)
organoid <- RenameIdents(organoid, "17" = "unknown")


?FindMarkers
markers <- FindMarkers(organoid, ident.1 = "15")
#### Add labels that do not overlap with each other:
# Install and load required packages

library(Seurat)
library(ggplot2)
library(ggrepel)

# Create the DimPlot without labels
p <- DimPlot(organoid, reduction = 'humap', label = FALSE)

# Add labels using geom_text_repel to avoid overlapping
p <- p + geom_text_repel(aes(label = organoid@meta.data$cluster), size = 5)

# Print the plot
print(p)




##### 14 will currently be listed as Astrocyte (just for Today / tomorrow)
##### 15 will be determined when Alexandre views the expr markers
##### 23 does not need to be annotated for right now due to it being 
## fGSEA -> emmaplot of each genotype. Comparing APOE3CH vs. APOE33, APOE4 vs. APOE33, APOE22 vs. APOE33
## Proportion of each cell type in each genotype (looking for differences) All clusters. Cluster 17.
## Compute percentage within each sample, percentage of each cell type. Boxplot by genotype.
## 3 samples by genotype. Aggregate these 3 points (each line of DF is a sample / celltype, cell column which is proportion,)
## (sample, cell_type/Cluster, percentage) (S1, Astr, 0.5%) (etc., etc. etc.)
## Give this to ggplot to create boxplots. Have a ggplot for Astrocyte, Gluta, etc.
## FindMarkers, ident.1 = 17, (you dont need to provide ident.2). Heatmap function from Seurat, put top 20 markers.
## Provide Alexandre with .csv files of 15 (just create it from the dataframe output from FindMarkers).
## 

library(pheatmap)

# Assume 'markers' is your dataframe from FindMarkers
markers <- FindMarkers(organoid, ident.1 = "unknown")

# Filter the markers based on the given criteria
filtered_markers <- markers %>%
  filter(abs(avg_log2FC) > 1.0, pct.1 > 0.25, p_val_adj < 0.001)

# Separate up-regulated and down-regulated genes
upregulated <- filtered_markers %>%
  filter(avg_log2FC > 0) %>%
  arrange(desc(abs(avg_log2FC))) %>%
  head(20)

downregulated <- filtered_markers %>%
  filter(avg_log2FC < 0) %>%
  arrange(desc(abs(avg_log2FC))) %>%
  head(20)

# Combine the top 20 up- and down-regulated genes into a single dataframe
top_genes <- bind_rows(upregulated, downregulated)

# Extract the avg_log2FC values
avg_log2FC_values <- top_genes$avg_log2FC
names(avg_log2FC_values) <- rownames(top_genes)

# Create a dataframe for the heatmap
heatmap_data <- data.frame(Gene = names(avg_log2FC_values), avg_log2FC = avg_log2FC_values)
rownames(heatmap_data) <- heatmap_data$Gene
heatmap_data <- heatmap_data[, -1, drop = FALSE]  # Remove the Gene column

pheatmap(heatmap_data,
         cluster_rows = TRUE,
         cluster_cols = FALSE,
         show_rownames = TRUE,
         show_colnames = FALSE,
         color = colorRampPalette(c("blue", "white", "red"))(100),
         main = "Top 20 Up- and Down-Regulated Genes Heatmap",
         cellwidth = 40,  # Adjust the cell width to make the heatmap more compact
         fontsize = 15,   # Adjust the font size to fit better
         display_numbers = TRUE,  # Display the avg_log2FC values
         number_format = "%.2f",  # Format numbers to two decimal places
         number_color = "white"   # Set the text color to white
)

saveRDS(organoid, file = paste0(dir, '/outputs/05-Cluster_Annotation/organoid1.rds'))
 
### Step 3. Subset Seurat Objects for Pseudobulk Prep                       ####


cellTypes <- c('Glutamatergic', 'GABAergic', 'Astrocyte', 'OPC', 
               'NPCs', 'VLMC', 'Pericyte Progenitor Cell', 'Pericyte', 'unknown', '23')

# Gluta
Glutamatergic <- subset(organoid, idents = "Glutamatergic")
saveRDS(Glutamatergic, file = paste0(dir, out, '/Glutamatergic.rds'))
rm(Glutamatergic)

# GABA
GABAergic <- subset(organoid, idents = "GABAergic")
saveRDS(GABAergic, file = paste0(dir, out, '/GABAergic.rds'))
rm(GABAergic)

# NPCs
NPCs <- subset(organoid, idents = "NPCs")
saveRDS(NPCs, file = paste0(dir, out, '/NPCs.rds'))
rm(NPCs)

# Imm. Astrocyte
ImmatureAstrocyte <- subset(organoid, idents = "Immature Astrocyte")
saveRDS(ImmatureAstrocyte, file = paste0(dir, out, '/ImmatureAstrocyte.rds'))
rm(ImmatureAstrocyte)

# Astrocyte
Astrocyte <- subset(organoid, idents = "Astrocyte")
saveRDS(Astrocyte, file = paste0(dir, out, '/Astrocyte.rds'))
rm(Astrocyte)

# OPC
OPC <- subset(organoid, idents = "OPC")
saveRDS(OPC, file = paste0(dir, out, '/OPC.rds'))
rm(OPC)

# unknown
unknown <- subset(organoid, idents = 'unknown')
saveRDS(unknown, file = paste0(dir, out, '/unknown.rds'))
rm(unknown)

# Pericyte
Pericyte <- subset(organoid, idents = "Pericyte")
saveRDS(Pericyte, file = paste0(dir, out, '/Pericyte.rds'))
rm(Pericyte)

# PericyteProgenitor
PericyteProgenitor <- subset(organoid, idents = "Pericyte Progenitor Cell")
saveRDS(PericyteProgenitor, file = paste0(dir, out, '/PericyteProgenitor.rds'))
rm(PericyteProgenitor)

# Gluta
cluster23 <- subset(organoid, idents = "23")
saveRDS(cluster23, file = paste0(dir, out, '/23.rds'))
rm(cluster23)

# VLMC
VLMC <- subset(organoid, idents = 'VLMC')
saveRDS(VLMC, file = paste0(dir, out, '/VLMC.rds'))