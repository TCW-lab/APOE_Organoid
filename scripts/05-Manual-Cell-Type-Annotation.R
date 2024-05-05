################################################################################
##################### Cell Type Annotation #####################################
################################################################################
## Read in organoid4

### Harmony Clusters
###
setwd('/projectnb/tcwlab/LabMember/akg/APOE_Jorganoid_project/')
# install.packages("harmony")
library(harmony)
library(Seurat)
organoid <- readRDS('outputs/03-CellTypeAnnotation/organoid_post_harmony.rds')

organoid <- organoid
rm(organoid)
organoid <- FindVariableFeatures(organoid)
organoid

organoid <- ScaleData(organoid)
organoid <- RunPCA(organoid)
## Run Harmony for Batch Correction:
organoid <- RunHarmony(organoid, group.by.vars = "sample")

# After Harmony integration, typically proceed with re-clustering
organoid <- FindNeighbors(organoid, reduction = "harmony")
organoid <- FindClusters(organoid)

# to visualize the results, you can run UMAP or another dimensionality reduction technique
organoid <- RunUMAP(organoid, dims = 1:50, reduction = "harmony", verbose = FALSE)

DimPlot(organoid, reduction = "umap", group.by = "genotype")


##### 02-02-2024 do the individual plots for all of these markers ##############
# organoid$seurat_clusters <- Idents(organoid)
organoid$harmony_clusters <- Idents(organoid)


Final_Markers_Progenitors = c("TOP2A", "MKI67", "FOXM1", "CENPF", "PAX6", "PCNA", "E2F2", "HOPX", "LHX2", "OTX2", "GLI3", "HES2", "HES6") 
Final_Markers_General_Neuron = c("DCX", "INSM1", "ST18", "NHLH1", "NEUROD1", "NEUROD6", "PRDM8", "MYT1L", "MAPT", "SNAP25", "SYP", "DLX6")
Final_Markers_Glut_Neuron = c("BCL11B", "DLG4", "STMN2", "SLC17A6", "GRIN2B", "RELN", "GRIA2")
Final_Markers_Gaba_Neuron = c("GAD1", "GAD2", "DLX1", "DLX2", "SST", "NRXN1", "ANK3") #rm "DLX6-AS1"
Final_Markers_Astrocyte = c("APOE", "GFAP", "PEA15", "S100B", "ALDH1L1", "FGFR3", "AGT", "AQP4")
Final_Markers_OPC = c("OLIG1", "OLIG2", "PCDH15", "LHFPL3", "MBP") #added MPC, only expr in oligos 
Final_Markers_VLMC = c("LUM", "PDGFRA", "DCN", "POSTN", "OGN", "APOD", "RSPO3") 
Final_Markers_Pericytes = c("CSPG4", "PDGFRB", "VTN", "CSPG4", "ACTA2", "KCNJ8", "ABCC9", "ACE2", "ART3", "ATP13A5") #rm CD146

## Visualize Canonical Cell Type Marker Expression
#Feature Plots

Final_Marker_List = list(Final_Markers_Progenitors, Final_Markers_General_Neuron, Final_Markers_Glut_Neuron, Final_Markers_Gaba_Neuron,
                         Final_Markers_Astrocyte, Final_Markers_OPC, Final_Markers_VLMC, Final_Markers_Pericytes)
Final_Marker_List_Names = c("Progenitor", "Neuron", "Glutamatergic", "GABA",
                            "Astrocyte", "OPC", "VLMC", "Pericyte")

dir <- "/projectnb/organoidlab/LabMember/akg/APOE_Jorganoid_project/outputs"
count = 1
for (x in Final_Marker_List) {
  for (y in x) {
    
    name = Final_Marker_List_Names[count]
    pdf(file = paste(dir, "/Feature_Plots/", name, "_", y, "_feature_plot.pdf", sep = ""), width = 15, height = 15)
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
  pdf(paste(dir, "/Dot_Plots/", name, "_dot_plot.pdf", sep = ""), width = 25, height = 15)
  
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


## Barcode rank plot

# Retrieve the total counts per cell
cell_counts <- organoid@meta.data$nCount_RNA

# Order cells by total counts
sorted_counts <- sort(cell_counts, decreasing = TRUE)

plot(log10(sorted_counts), pch = 19, cex = 0.5, main = "Barcode Rank Plot", 
     xlab = "Cell Rank", ylab = "Log10 Total Counts")

