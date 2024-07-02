(1.) First, perform alignment using Simpleaf quant in order to generate mat (matrix) files for each sample, which contain read alignment data. 
Convert these matrices to 24 SingleCellExperiment Objects, then perform nCount_RNA cutoff to remove cells with fewer than ~1000 mRNAs detected. We then 
(2.) take the 24 individual Seurat objects, perform miQC (in order to QC for percent.mt and nCount_RNA (upper bounds)) and perform 
DoubletFinder (in order to remove possible doublets) (total of about 20% of cells removed from combined steps). The next step is 
(3.) using SoupX to correct for ambientRNA molecules (~3.5% of molecules removed). We then combine the 24 samples Seurat files into one, then 
(4.) perform Normalization, FindVariableFeatures, ScaleData, PCA, UMAP, FindClusters, Harmony batch-correction, then we run FindNeighbors, FindClusters, and UMAP (for Harmony reduction). 
(5.) We then annotate the clusters, based on cell marker expression profiles. 
(6.) In an attempt to more robustly assess effects of APOE genotype on expression profiles, we perform Pseudobulk / DESeq2 analysis on the data from each cell type, this involves subsetting the organoid object for each specific cell type, then performing DESeq2 analysis. Finally, we 
(7.) Perform fGSEA for identifying enriched pathways, validation, and hypothesis generation.
