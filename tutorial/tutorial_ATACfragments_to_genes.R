# ==============================================================================
# Tutorial: Generating a Gene Activity Matrix from scATAC-seq fragments using ArchR
#
# Goal: Convert ATAC fragment files from multiple samples into a single
#       "cells x genes" gene activity matrix (GeneScoreMatrix in ArchR),
#       for downstream integration with scRNA-seq (e.g., scPositioner).
#
# Core idea: When creating Arrow files, setting addGeneScoreMat = TRUE makes
#            ArchR automatically compute, for each gene, a weighted fragment
#            count over promoter + gene body regions, approximating expression.
#
# Requirements: R >= 4.0, with ArchR and Matrix installed.
# ==============================================================================

# ------------------------------------------------------------------------------
# Step 0: Environment setup and working directory
# ------------------------------------------------------------------------------
setwd('your/path/to/ArchR/')

library(ArchR)
library(Matrix)

# Set multithreading (adjust based on available cores; 16~32 recommended for large datasets)
addArchRThreads(threads = 32)


# ------------------------------------------------------------------------------
# Step 1: Prepare example data and reference genome
# We use two small samples (B001-A-006 and B001-A-302) to demonstrate the full pipeline. 
# Demo data can be downloaded via https://drive.google.com/drive/folders/1DntqOZcHvWIsGM990ZPKxUtLhh-p4rw4?usp=drive_link
# ------------------------------------------------------------------------------

# ---- B001 batch (using 2 samples as demo) ----
b01006 <- "../B001-A-006/B001-A-006_atac_fragments.tsv.gz"
b01302 <- "../B001-A-302/B001-A-302_atac_fragments.tsv.gz"

# Combine into a named vector (order must match names below)
paths <- c(b01006, b01302)
names(paths) <- c("b01006", "b01302")

cat(sprintf("Defined %d sample path(s).\n", length(paths)))

addArchRGenome("hg38")  # Built-in human genome annotation in ArchR


# ------------------------------------------------------------------------------
# Step 2: Create Arrow files from fragments (core step)
# ------------------------------------------------------------------------------
ArrowFiles <- ArchR::createArrowFiles(
  inputFiles = paths,
  sampleNames = names(paths),
  filterTSS = 4,          # Min TSS enrichment score per cell (QC)
  filterFrags = 1000,     # Drop low-quality cells
  addTileMat = TRUE,      # Tile matrix for downstream dim reduction (optional)
  minFrags = 500,         # Min fragments per cell
  maxFrags = 1e+05,       # Max fragments per cell (remove doublets/large nuclei)
  addGeneScoreMat = TRUE, # ★ Automatically computes GeneScoreMatrix
  nChunk = 1,
  threads = 32
)


# ------------------------------------------------------------------------------
# Step 3: Create an ArchRProject object
# ------------------------------------------------------------------------------
proj <- ArchRProject(
  ArrowFiles = ArrowFiles,
  outputDirectory = "Project_output",
  copyArrows = TRUE  # Keep an unmodified copy of Arrow files for reuse
)

# Check available matrices (should include "GeneScoreMatrix")
print(getAvailableMatrices(proj))

use_proj <- proj


# ------------------------------------------------------------------------------
# Step 4: Extract the GeneScoreMatrix
# ------------------------------------------------------------------------------
genescore <- getMatrixFromProject(
  ArchRProj = use_proj,
  useMatrix = "GeneScoreMatrix",
  useSeqnames = NULL,
  verbose = TRUE,
  binarize = FALSE,       # FALSE = continuous scores (recommended)
  threads = getArchRThreads(),
  logFile = createLogFile("getMatrixFromProject")
)


# ------------------------------------------------------------------------------
# Step 5: Format conversion and export
# ------------------------------------------------------------------------------
# Transpose: genes x cells -> cells x genes, then convert to dgCMatrix sparse format
genescore_matrix <- t(genescore@assays@data$GeneScoreMatrix)
genescore_matrix <- as(genescore_matrix, "dgCMatrix")

# Get gene names (features)
genescore_genes <- getFeatures(use_proj, useMatrix = "GeneScoreMatrix")

# Get cell metadata (e.g., Sample, cellNames), aligned to matrix row names
genescore_meta <- data.frame(use_proj@cellColData)
genescore_meta <- genescore_meta[rownames(genescore_matrix), ]

# Create output dir and export
dir.create("../processed", showWarnings = FALSE)
writeMM(genescore_matrix, "../processed/genescore.txt")   # cells x genes matrix
write.csv(genescore_meta, "../processed/meta.csv")        # cell metadata
write.csv(genescore_genes, "../processed/genenames.csv")  # gene list

cat("Gene activity matrix exported to ../processed/\n")
cat(sprintf("Matrix dimensions: %d cells x %d genes\n", nrow(genescore_matrix), ncol(genescore_matrix)))


# ------------------------------------------------------------------------------
# Step 6: Load in Python and convert to h5ad (AnnData)
# Run the following code in Python (requires scanpy / anndata / scipy)
# ------------------------------------------------------------------------------
# import scipy.io
# import pandas as pd
# import anndata as ad
#
# # Load the sparse matrix exported from R (cells x genes, MatrixMarket format)
# # NOTE: R's writeMM does NOT store row/column names, so we assign them manually.
# mat = scipy.io.mmread("../processed/genescore.txt").tocsr()
#
# # Load cell metadata and gene list
# meta = pd.read_csv("../processed/meta.csv", index_col=0)
# genes = pd.read_csv("../processed/genenames.csv", index_col=0)
#
# # Sanity checks: ensure dimensionality matches
# assert mat.shape[0] == meta.shape[0], "Row count mismatch between matrix and metadata"
# assert mat.shape[1] == genes.shape[0], "Column count mismatch between matrix and gene list"
#
# # Build AnnData object
# # genescore_genes is a single-column CSV; take its first column as var_names
# adata = ad.AnnData(X=mat)
# adata.obs = meta        # meta.index (cell names) becomes obs_names automatically
# adata.var_names = genes.iloc[:, 0].values
#
# # Save as h5ad for downstream use (e.g., scPositioner, scanpy)
# adata.write("../processed/genescore.h5ad")
#
# print(f"AnnData created: {adata.shape[0]} cells x {adata.shape[1]} genes")
# print("Saved to ../processed/genescore.h5ad")
