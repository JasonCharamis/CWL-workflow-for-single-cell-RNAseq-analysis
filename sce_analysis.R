
## Script for basic scRNAseq data analysis. Includes:

# 1. Quality Control (QC) of cells
# 2. Normalization of gene expression
# 3. Feature selection
# 4. PCA
# 5. Clustering and tSNE, UMAP visualizations
# 6. Marker gene detection

## Load libraries
suppressPackageStartupMessages(library(seurat))
suppressPackageStartupMessages(library(ggplot2))


## Get input file from standard input ------------------------------------------------------------------
args <- commandArgs(trailingOnly = TRUE)
stopifnot(length(args) > 0, file.exists(args))

## Read file as SingleCellExperimentClass
sce <- Read10X(args[1])

# Quality control based on MT genes percent
mito_genes <- grep("^MT-", rownames(sce), value = TRUE)
sce <- AddMetaData(sce, metadata = perCellQCMetrics(sce, subsets = list(Mito = mito_genes)))

# Normalization
sce <- NormalizeData(sce)

# Feature selection
sce <- FindVariableFeatures(sce, selection.method = "vst", nfeatures = 2000)

# PCA
sce <- RunPCA(sce, features = VariableFeatures(object = sce), npcs = 25)

# Clustering
sce <- FindNeighbors(sce, dims = 1:25)
sce <- FindClusters(sce, resolution = 0.5)

# tSNE and UMAP
tSNE <- RunTSNE(sce, dims = 1:25)
UMAP <- RunUMAP(sce, dims = 1:25)

ggsave(plot=tSNE,"tSNE.svg", dpi = 1000)
ggsave(plot=UMAP,"UMAP.svg", dpi = 1000)

# Marker gene detection
markers <- FindAllMarkers(sce, only.pos = TRUE, min.pct = 0.25)
write.table(markers, file = "Marker_genes_detected.txt", quote = FALSE, sep = "\t")
