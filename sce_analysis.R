
## Script for basic scRNAseq data analysis. Includes:

# 1. Quality Control (QC) of cells
# 2. Normalization of gene expression
# 3. Feature selection
# 4. PCA
# 5. Clustering and tSNE, UMAP visualizations
# 6. Marker gene detection

## Load libraries
suppressPackageStartupMessages(library(scRNAseq))
suppressPackageStartupMessages(library(scater))
suppressPackageStartupMessages(library(scran))
suppressPackageStartupMessages(library(bluster))
suppressPackageStartupMessages(library(ggplot2))
suppressPackageStartupMessages(library(DropletUtils))


## Get input file from standard input ------------------------------------------------------------------
args <- commandArgs(trailingOnly = TRUE)
stopifnot(length(args) > 0, file.exists(args))

## Read file as SingleCellExperimentClass
sce <- read10xCounts(args[1])

# Quality control based on MT genes percent
is.mito <- grepl("^MT-", rownames(sce))
qcstats <- perCellQCMetrics(sce, subsets=list(Mito=is.mito))
filtered <- quickPerCellQC(qcstats, percent_subsets="subsets_Mito_percent")

# Normalization
set.seed(1234)
clusters <- quickCluster(sce)
sce <- computeSumFactors(sce, cluster=clusters) 
sce <- logNormCounts(sce)

# Feature selection
dec <- modelGeneVar(sce)
hvg <- getTopHVGs(dec, prop=0.1)

# PCA with 25 components (default)
set.seed(1234)
sce <- runPCA(sce, ncomponents=25, subset_row=hvg)
ggsave(plot=pca,"PCA.svg", dpi = 1000)

# Clustering
snn.gr <- buildSNNGraph(sce, use.dimred="PCA")
colLabels(sce) <- factor(igraph::cluster_walktrap(snn.gr)$membership)

tSNE <- plotTSNE(sce, colour_by="label")
UMAP <- plotUMAP(sce, colour_by="label")

ggsave(plot=tSNE,"tSNE.svg", dpi = 1000)
ggsave(plot=UMAP,"UMAP.svg", dpi = 1000)

## Marker gene detection
markers <- findMarkers(sce, direction="up")
write.table("Marker_genes_detected.txt", markers, quote = F)
