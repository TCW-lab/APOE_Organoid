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

organoid <- readRDS(paste0(dir, '/outputs/05-Cluster_Annotation', '/organoid_unannotated.rds'))
# ^ organoid is the original, unannotated ^

unique(organoid$seurat_clusters)

DimPlot(organoid, reduction = 'humap', group.by = 'seurat_clusters', label = TRUE,
        label.size = 4)

DimPlot(organoid, reduction = 'humap', group.by = 'seurat_clusters', label = TRUE)

Idents(organoid) <- organoid$seurat_clusters
organoid$seurat_clusters <- Idents(organoid)
unique(Idents(organoid))
unique(organoid$seurat_clusters)

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
Final_Markers_Radial_Glia = c('PAX6', 'SOX2', 'MEIS2', 'FOXP2', "TLE4", "VIM")
Final_Markers_Neural_Crest = c('FOXD3', 'SNAI1', 'SNAI2', 'SOX8', 'SOX9', 'SOX10')
## Visualize Canonical Cell Type Marker Expression
#Feature Plots

Final_Marker_List = list(Final_Markers_Progenitors, Final_Markers_General_Neuron, Final_Markers_Glut_Neuron, Final_Markers_Gaba_Neuron,
                         Final_Markers_Astrocyte, Final_Markers_OPC, Final_Markers_VLMC, Final_Markers_Pericytes,
                         Final_Markers_PigEp, Final_Markers_Chloroid, Final_Markers_Radial_Glia, Final_Markers_Neural_Crest)
Final_Marker_List_Names = c("Progenitor", "Neuron", "Glutamatergic", "GABA",
                            "Astrocyte", "OPC", "VLMC", "Pericyte",
                            "PigEp", "Chloroid", "Radial_Glia", "Neural_Crest")

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
              Final_Markers_PigEp, Final_Markers_Chloroid, Final_Markers_Radial_Glia, Final_Markers_Neural_Crest)
marker_names <- list("Progenitors", "General_Neuron", "Glut", "GABA", "Astrocyte", 
                     "OPC", "VLMC", "Pericyte", "PigEp", "Chloroid", "Radial_Glia", "Neural_Crest")


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
DimPlot(organoid, reduction = "humap", label = TRUE, label.size = 4)


#Idents (For Reference):
c("GABAergic", "GABAergic", 'Glutamatergic', "NPCs", "Astrocyte", 
  "Astrocyte", "OPC/Oligo", "GABAergic", 'Glutamatergic', "Astroglia",
  'Glutamatergic', "NPCs (cycling)", "NPCs (cycling)", 'Pericyte Progenitor', 'Astrocyte', 
  "NPCs (cycling)", 'unknown1', "Astrocyte", "NPCs (cycling)", "VLMC",
  "NPCs", "Choroid Plexus", 'unknown2')

### Top level of annotation: Astrocyte (4.33x1.5)
### 18, 3, 4
### 
organoid <- RenameIdents(organoid, "3" = "Astrocyte")
organoid <- RenameIdents(organoid, "4" = "Astrocyte")
organoid <- RenameIdents(organoid, '15' = 'Astrocyte')
organoid <- RenameIdents(organoid, "18" = "Astrocyte")
DimPlot(organoid, reduction = "humap", label = TRUE, label.size = 4)
# 7 is Astroglia
organoid <- RenameIdents(organoid, "7" = "Astroglia")


## GABAergic:
# 0,2,5
organoid <- RenameIdents(organoid, "0" = "GABAergic")
organoid <- RenameIdents(organoid, "2" = "GABAergic")
organoid <- RenameIdents(organoid, "5" = "GABAergic")

DimPlot(organoid, reduction = "humap", label = TRUE, label.size = 4)

## Glutamatergic:
# 1, 8, 9
organoid <- RenameIdents(organoid, "1" = "Glutamatergic")
organoid <- RenameIdents(organoid, "9" = "Glutamatergic")
organoid <- RenameIdents(organoid, "8" = "Glutamatergic")
# DimPlot(organoid, reduction = "humap", label = TRUE, label.size = 5)

## Choroid Plexus
# 21
# organoid <- RenameIdents(organoid, "21" = "Choroid Plexus") # renaming this Pericyte


## NPCs:
# NPCs (cycling) 13, 16, 19
organoid <- RenameIdents(organoid, "11" = "NPCs (cycling)")
organoid <- RenameIdents(organoid, "13" = "NPCs (cycling)")
organoid <- RenameIdents(organoid, "16" = "NPCs (cycling)")
DimPlot(organoid, reduction = 'humap', label = TRUE, label.size = 4)
organoid <- RenameIdents(organoid, "12" = "NPCs")
# NPCs (non-cycling): 
# organoid <- RenameIdents(organoid, "9" = "NPCs")
# organoid <- RenameIdents(organoid, "11" = "NPCs")

## OPCs / Oligo:
# 5
organoid <- RenameIdents(organoid, "6" = "OPC / Oligo")
organoid <- RenameIdents(organoid, "10" = "OPC / Oligo")
DimPlot(organoid, reduction = 'humap', label = TRUE, label.size = 4)

## Pericyte (progenitor cells):
# 20
organoid <- RenameIdents(organoid, '20' = 'Pericyte') 
# Pericyte Progenitor 14
organoid <- RenameIdents(organoid, '14' = 'Pericyte Progenitor') 

DimPlot(organoid, reduction = 'humap', label = TRUE, label.size = 4)
## VLMC:
# 19
organoid <- RenameIdents(organoid, "19" = "VLMC")
# DimPlot(organoid, reduction = 'humap', label = TRUE, label.size = 4.5)

## Pigmented Epithelial
# 22 for certain.
# organoid <- RenameIdents(organoid, "22" = "Pigmented Epithelial")
## Chloroid Plexus
# Could be 20, but 20 seems more likely to be VLMC

## Unknown cluster (17)
organoid <- RenameIdents(organoid, "17" = "unknown")


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
# DimPlot(organoid, reduction = 'humap', label = TRUE, size = 4)

saveRDS(organoid, file = paste0(dir, out, '/organoid_annotated.rds'))

### Step 3. Visualization for certain metrics                               ####

organoid <- readRDS(paste0(dir, out, '/organoid_annotated.rds'))

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

### Visualize the expression across each of the 4 genotypes:
library(ggplot2)
library(gridExtra)
# Ensure 'genotype' is a factor with correct levels
organoid$genotype <- factor(organoid$genotype, levels = c("APOE22", "APOE33", "APOE33Ch", "APOE44"))

# Create a base plot with all cells in light gray
base_plot <- DimPlot(organoid, reduction = "humap", cells = colnames(organoid), 
                     cols = "lightgray") + 
  theme_minimal()

# Function to create individual plots for each genotype
create_genotype_plot <- function(genotype) {
  DimPlot(organoid, reduction = "humap", group.by = "genotype", cells.highlight = WhichCells(organoid, expression = genotype == !!genotype), 
          cols.highlight = "blue", cols = "lightgray") + 
    ggtitle(genotype) + 
    theme_minimal() +
    NoLegend()
}

# Create individual plots for each genotype
plot_apoe22 <- create_genotype_plot("APOE22")
plot_apoe33 <- create_genotype_plot("APOE33")
plot_apoe33ch <- create_genotype_plot("APOE33Ch")
plot_apoe44 <- create_genotype_plot("APOE44")

# Combine the plots into a 4-panel plot
combined_plot <- grid.arrange(plot_apoe22, plot_apoe33, plot_apoe33ch, plot_apoe44, ncol = 2)

# Display the combined plot
print(combined_plot)






library(Seurat)
library(ggplot2)
library(gridExtra)

# Check the available reductions
print(Reductions(organoid))

# Create individual plots for each genotype
plot_apoe2 <- FeaturePlot(organoid, reduction = "humap", features = "genotype", 
                          cells = WhichCells(organoid, expression = genotype == "APOE22")) + 
  ggtitle("APOE22")

plot_apoe3 <- FeaturePlot(organoid, reduction = "humap", features = "genotype", 
                          cells = WhichCells(organoid, expression = genotype == "APOE33")) + 
  ggtitle("APOE33")

plot_apoe3ch <- FeaturePlot(organoid, reduction = "humap", features = "genotype", 
                            cells = WhichCells(organoid, expression = genotype == "APOE33Ch")) + 
  ggtitle("APOE33Ch")

plot_apoe4 <- FeaturePlot(organoid, reduction = "humap", features = "genotype", 
                          cells = WhichCells(organoid, expression = genotype == "APOE44")) + 
  ggtitle("APOE44")

# Combine the plots into a 4-panel plot
grid.arrange(plot_apoe2, plot_apoe3, plot_apoe3ch, plot_apoe4, ncol = 2)




# unique(organoid$seurat_clusters)



organoid <- readRDS(paste0(dir, '/outputs/05-Cluster_Annotation/organoid1.rds'))
# saveRDS(organoid, file = paste0(dir, '/outputs/05-Cluster_Annotation/organoid1.rds'))

### Step 4. Subset Seurat Objects for Pseudobulk Prep                       ####


cellTypes <- c('Glutamatergic', 'GABAergic', 'Astrocyte', 'OPC', 
               'NPCs', 'VLMC', 'Pericyte Progenitor Cell', 'Pericyte', 'unknown')

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

# NPCs (cycling)
NPCsCycling <- subset(organoid, idents = "NPCs (cycling)")
saveRDS(NPCsCycling, file = paste0(dir, out, '/NPCsCycling.rds'))
rm(NPCsCycling)




# Imm. Astrocyte
Astroglia <- subset(organoid, idents = "Astroglia")
saveRDS(Astroglia, file = paste0(dir, out, '/Astroglia.rds'))
rm(Astroglia)

# Astrocyte
Astrocyte <- subset(organoid, idents = "Astrocyte")
saveRDS(Astrocyte, file = paste0(dir, out, '/Astrocyte.rds'))
rm(Astrocyte)

# OPC
OPC_Oligo <- subset(organoid, idents = "OPC / Oligo")
saveRDS(OPC_Oligo, file = paste0(dir, out, '/OPC_Oligo.rds'))
rm(OPC_Oligo)


# Pericyte
Pericyte <- subset(organoid, idents = "Pericyte")
saveRDS(Pericyte, file = paste0(dir, out, '/Pericyte.rds'))
rm(Pericyte)

# PericyteProgenitor
PericyteProgenitor <- subset(organoid, idents = "Pericyte Progenitor")
saveRDS(PericyteProgenitor, file = paste0(dir, out, '/PericyteProgenitor.rds'))
rm(PericyteProgenitor)

# unknown1
unknown1 <- subset(organoid, idents = "Unknown")
saveRDS(unknown1, file = paste0(dir, out, '/Unknown.rds'))
rm(Unknown)

# VLMC
VLMC <- subset(organoid, idents = 'VLMC')
saveRDS(VLMC, file = paste0(dir, out, '/VLMC.rds'))
rm(VLMC)
#### Plot by genotype:


DimPlot(organoid, group.by = 'genotype', reduction = 'humap')

organoid$seurat_clusters <- Idents(organoid)
# DimPlot(organoid, reduction = 'humap', label = TRUE)
library(patchwork)
# Calculate the total number of cells per cluster
cluster_counts <- table(organoid@meta.data$seurat_clusters)

# Order the clusters by the number of cells in descending order
ordered_clusters <- names(sort(cluster_counts, decreasing = TRUE))
custom_theme <- theme(
  legend.position = "none",        # Remove legend
  panel.grid.major = element_line(color = "gray", size = 0.25),  # Major grid lines
  panel.grid.minor = element_line(color = "gray", size = 0.25), # Minor grid lines
  panel.background = element_blank(),  # Remove panel background
  plot.title = element_text(hjust = 0.5) # Center the plot titles
)


# Plot nCount_RNA with a title and no legend, ordered by the number of cells per cluster
plot1 <- VlnPlot(organoid, features = "nCount_RNA", group.by = "seurat_clusters", pt.size = 0) +
  ggtitle("nCount_RNA Distribution Across Seurat Clusters") +
  custom_theme +
  scale_x_discrete(limits = ordered_clusters) +
  scale_y_log10()

# Plot percent.mt with a title and no legend, ordered by the number of cells per cluster
plot2 <- VlnPlot(organoid, features = "percent.mt", group.by = "seurat_clusters", pt.size = 0) +
  ggtitle("percent.mt Distribution Across Seurat Clusters") +
  custom_theme +
  scale_x_discrete(limits = ordered_clusters) +
  scale_y_log10()

# Plot nFeature_RNA with a title and no legend, ordered by the number of cells per cluster
# Plot nFeature_RNA with a title, no legend, background grid, and custom log scale
plot3 <- VlnPlot(organoid, features = "nFeature_RNA", group.by = "seurat_clusters", pt.size = 0) +
  ggtitle("nFeature_RNA Distribution Across Seurat Clusters") +
  custom_theme +
  scale_x_discrete(limits = ordered_clusters) +
  scale_y_log10(breaks = c(100, 1000, 10000), labels = c("1e+02", "1e+03", "1e+04"))


print(plot3)

# Combine plots
combined_plot <- plot1 + plot2 + plot3
print(combined_plot)











