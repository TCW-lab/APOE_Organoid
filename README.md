## Analysis Workflow

1. **Alignment and Matrix Generation**
   - Perform alignment using Simpleaf quant to generate `mat` (matrix) files for each sample, containing read alignment data.
   - Convert these matrices into 24 SingleCellExperiment objects.
   - Apply `nCount_RNA` cutoff to remove cells with fewer than ~1000 mRNAs detected.

2. **Seurat Object Creation and Quality Control**
   - Convert the 24 SingleCellExperiment objects into individual Seurat objects.
   - Perform `miQC` to QC for `percent.mt` and `nCount_RNA` (upper bounds).
   - Apply `DoubletFinder` to remove possible doublets (total of about 20% of cells removed from combined steps).

3. **Ambient RNA Correction**
   - Use `SoupX` to correct for ambient RNA molecules (~3.5% of molecules removed).
   - Combine the 24 Seurat objects into one Seurat object.

4. **Normalization and Clustering**
   - Perform `Normalization`, `FindVariableFeatures`, `ScaleData`, `PCA`, `UMAP`, and `FindClusters`.
   - Apply `Harmony` for batch correction.
   - Run `FindNeighbors`, `FindClusters`, and `UMAP` (for Harmony reduction).

5. **Cluster Annotation**
   - Annotate clusters based on cell marker expression profiles.

6. **Pseudobulk / DESeq2 Analysis**
   - To assess the effects of APOE genotype on expression profiles, perform Pseudobulk / DESeq2 analysis for each cell type.
   - Subset the organoid object for each specific cell type and perform DESeq2 analysis.

7. **Pathway Enrichment Analysis**
   - Perform `fGSEA` for identifying enriched pathways, validation, and hypothesis generation.
