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



BiocManager::install(version = "3.18")

setwd("/projectnb/tcwlab/LabMember/akg/projects/APOE_Organoid")

dir <- "/projectnb/tcwlab/LabMember/akg/projects/APOE_Organoid"
out <- '/outputs/06-Cluster_Annotation'


### Step 1. Read in RDS object and plot FeaturePlots for various markers   #####

organoid <- readRDS(paste0(dir, '/outputs/05-Doublet_Finder', '/organoid.rds'))

organoid <- Seurat::SCTransform(organoid, vars.to.regress = 'percent.mt', vst.flavor="v2",
                                useNames = FALSE)

unique(organoid$RNA_snn_res.0.6)
saveRDS(organoid, file = paste0(dir, '/outputs/05-Doublet_Finder/organoid.rds'))







### First inspect the DimPlot of the humap to see what we're working with:
DimPlot(organoid, reduction = "umap_rotated", group.by = 'seurat_clusters', label = TRUE, label.size = 6) + 
  labs(color = 'Seurat Clusters') + 
  ggtitle('UMAP DimPlot')


library(Seurat)

# Extract UMAP coordinates
umap_coords <- Embeddings(organoid, "umap")

# Apply 90 degrees clockwise rotation
rotation_matrix <- matrix(c(0, 1, -1, 0), nrow = 2)
rotated_coords <- umap_coords %*% rotation_matrix

# Update the Seurat object with the rotated coordinates
organoid[['umap_rotated']] <- CreateDimReducObject(embeddings = rotated_coords, key = "UMAP_", assay = DefaultAssay(organoid))

# Plot the rotated UMAP
DimPlot(organoid, reduction = "umap_rotated", group.by = 'seurat_clusters', label = TRUE, label.size = 6)






organoid$harmony_clusters <- Idents(organoid)


Final_Markers_Progenitors = c("TOP2A", "MKI67", "FOXM1", "CENPF", "PAX6", "PCNA", "E2F2", "HOPX", "LHX2", "OTX2", "GLI3", "HES2", "HES6") 
Final_Markers_General_Neuron = c("DCX", "INSM1", "ST18", "NHLH1", "NEUROD1", "NEUROD6", "PRDM8", "MYT1L", "MAPT", "SNAP25", "SYP", "DLX6")
Final_Markers_Glut_Neuron = c("BCL11B", "DLG4", "STMN2", "SLC17A6", "GRIN2B", "RELN", "GRIA2")
Final_Markers_Gaba_Neuron = c("GAD1", "GAD2", "DLX1", "DLX2", "SST", "NRXN1", "ANK3") #rm "DLX6-AS1"
Final_Markers_Astrocyte = c("APOE", "GFAP", "PEA15", "S100B", "ALDH1L1", "FGFR3", "AGT", "AQP4")
Final_Markers_OPC = c("OLIG1", "OLIG2", "PCDH15", "LHFPL3", "MBP") #added MPC, only expr in oligos 
Final_Markers_VLMC = c("LUM", "PDGFRA", "DCN", "POSTN", "OGN", "APOD", "RSPO3") 
Final_Markers_Pericytes = c("CSPG4", "PDGFRB", "VTN", "ACTA2", "KCNJ8", "ABCC9", "ACE2", "ART3", "ATP13A5") #rm CD146

## Visualize Canonical Cell Type Marker Expression
#Feature Plots

Final_Marker_List = list(Final_Markers_Progenitors, Final_Markers_General_Neuron, Final_Markers_Glut_Neuron, Final_Markers_Gaba_Neuron,
                         Final_Markers_Astrocyte, Final_Markers_OPC, Final_Markers_VLMC, Final_Markers_Pericytes)
Final_Marker_List_Names = c("Progenitor", "Neuron", "Glutamatergic", "GABA",
                            "Astrocyte", "OPC", "VLMC", "Pericyte")

count = 1
for (x in Final_Marker_List) {
  for (y in x) {
    
    name = Final_Marker_List_Names[count]
    pdf(file = paste(dir, "/outputs/images/05-Cluster_Annotation/Feature_Plots/", name, "_", y, "_feature_plot.pdf", sep = ""), width = 15, height = 15)
    print(FeaturePlot(organoid, reduction = "humap", features = y, label = TRUE, label.size = 6))
    dev.off()
  }
  count = count + 1
}

names(organoid@reductions)
FeaturePlot(organoid, reduction = "humap", label = TRUE, label.size = 6)
# Feature plots should be similar to below:
#Dot Plots 

Final_Markers_Neuron[4:24]

feats <- list(Final_Markers_Progenitors, Final_Markers_General_Neuron, Final_Markers_Glut_Neuron, Final_Markers_Gaba_Neuron,
              Final_Markers_Astrocyte, Final_Markers_OPC, Final_Markers_VLMC, Final_Markers_Pericytes)
marker_names <- list("Progenitors", "General_Neuron", "Glut", "GABA", "Astrocyte", "OPC", "VLMC", "Pericyte")

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
  pdf(paste(dir, "/outputs/images/05-Cluster_Annotation/Dot_Plots/", name, "_dot_plot.pdf", sep = ""), width = 25, height = 15)
  
  print(DotPlot(
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
  )
  )
  
  dev.off()
}
