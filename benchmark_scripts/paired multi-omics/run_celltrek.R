setwd("C:/Users/11437/Desktop/school/scpositioner/paired_omics")
library(Seurat)
library(Matrix)
library(SeuratObject)
library(CellTrek)

rna_expr <- read.csv("human_rna_expr.csv", row.names = 1, check.names = FALSE)
rna_meta <- read.csv("human_rna_meta.csv", row.names = 1)
atac_expr <- read.csv("human_atac_expr.csv", row.names = 1, check.names = FALSE)
atac_meta <- read.csv("human_atac_meta.csv", row.names = 1)

rownames(atac_expr) <- gsub("^atac-", "", rownames(atac_expr))
rownames(rna_expr) <- as.character(rownames(rna_expr))
common_features <- intersect(rownames(rna_expr), rownames(atac_expr))

message("shared features: ", length(common_features))

rna_expr <- rna_expr[common_features, , drop = FALSE]
atac_expr <- atac_expr[common_features, , drop = FALSE]

rna_mat <- Matrix(as.matrix(rna_expr), sparse = TRUE)

sc_adata <- CreateSeuratObject(
  counts = rna_mat,
  assay = "RNA"
)

sc_adata <- AddMetaData(sc_adata, rna_meta)
sc_adata$celltype <- "celltype"
sc_adata$orig.ident <- "scRNA"
atac_mat <- Matrix(as.matrix(atac_expr), sparse = TRUE)

st_adata <- CreateSeuratObject(
  counts = atac_mat,
  assay = "RNA"
)

st_adata <- AddMetaData(st_adata, atac_meta)

st_adata$orig.ident <- "ATAC_spatial"
st_adata$celltype <- "celltype"

coords <- data.frame(
  tissue = 1,
  row = atac_meta$x,
  col = atac_meta$y,
  imagerow = atac_meta$x,
  imagecol = atac_meta$y
)

rownames(coords) <- rownames(atac_meta)


colnames(st_adata) <- paste0(colnames(st_adata), "-1")
rownames(atac_meta) <- colnames(st_adata)
rownames(coords) <- colnames(st_adata)

fake_image <- array(0, dim = c(10, 10, 3))

sf <- scalefactors(
  spot = 1,
  fiducial = 1,
  hires = 1,
  lowres = 1
)

st_adata@images$slice1 <- new(
  Class = "VisiumV1",
  image = fake_image,
  coordinates = coords,
  scale.factors = sf,
  assay = "RNA",
  key = "slice1_"
)

sc_adata <- NormalizeData(sc_adata)
sc_adata <- FindVariableFeatures(sc_adata)
sc_adata <- ScaleData(sc_adata)
sc_adata <- RunPCA(sc_adata, npcs = 30)

st_adata <- NormalizeData(st_adata)
st_adata <- FindVariableFeatures(st_adata)
st_adata <- ScaleData(st_adata)
st_adata <- RunPCA(st_adata, npcs = 30)

traint <- CellTrek::traint(
  st_data = st_adata,
  sc_data = sc_adata,
  st_assay = "RNA",
  sc_assay = "RNA",
  cell_names = "celltype"
)

celltrek_obj <- CellTrek::celltrek(
  st_sc_int = traint,
  int_assay = "traint",
  sc_data = sc_adata,
  sc_assay = "RNA",
  reduction = "pca",
  intp = TRUE,
  intp_pnt = 5000,
  intp_lin = FALSE,
  nPCs = 30,
  ntree = 1000,
  dist_thresh = 0.55,
  top_spot = 5,
  spot_n = 50,
  repel_r = 0.0001,
  repel_iter = 20,
  keep_model = TRUE
)$celltrek

write.csv(celltrek_obj@meta.data, file = "celltrek_human.csv",row.names = TRUE)


#mouse p22
rna_expr <- read.csv("mousep22_rna_expr.csv", row.names = 1, check.names = FALSE)
rna_meta <- read.csv("mousep22_rna_meta.csv", row.names = 1)
atac_expr <- read.csv("mousep22_atac_expr.csv", row.names = 1, check.names = FALSE)
atac_meta <- read.csv("mousep22_atac_meta.csv", row.names = 1)

rownames(atac_expr) <- gsub("^atac-", "", rownames(atac_expr))
rownames(rna_expr) <- as.character(rownames(rna_expr))
common_features <- intersect(rownames(rna_expr), rownames(atac_expr))

message("shared features: ", length(common_features))

rna_expr <- rna_expr[common_features, , drop = FALSE]
atac_expr <- atac_expr[common_features, , drop = FALSE]

rna_mat <- Matrix(as.matrix(rna_expr), sparse = TRUE)

sc_adata <- CreateSeuratObject(
  counts = rna_mat,
  assay = "RNA"
)

sc_adata <- AddMetaData(sc_adata, rna_meta)
sc_adata$celltype <- "celltype"
sc_adata$orig.ident <- "scRNA"
atac_mat <- Matrix(as.matrix(atac_expr), sparse = TRUE)

st_adata <- CreateSeuratObject(
  counts = atac_mat,
  assay = "RNA"
)

st_adata <- AddMetaData(st_adata, atac_meta)

st_adata$orig.ident <- "ATAC_spatial"
st_adata$celltype <- "celltype"

coords <- data.frame(
  tissue = 1,
  row = atac_meta$x,
  col = atac_meta$y,
  imagerow = atac_meta$x,
  imagecol = atac_meta$y
)

rownames(coords) <- rownames(atac_meta)


colnames(st_adata) <- paste0(colnames(st_adata), "-1")
rownames(atac_meta) <- colnames(st_adata)
rownames(coords) <- colnames(st_adata)

fake_image <- array(0, dim = c(10, 10, 3))

sf <- scalefactors(
  spot = 1,
  fiducial = 1,
  hires = 1,
  lowres = 1
)

st_adata@images$slice1 <- new(
  Class = "VisiumV1",
  image = fake_image,
  coordinates = coords,
  scale.factors = sf,
  assay = "RNA",
  key = "slice1_"
)

sc_adata <- NormalizeData(sc_adata)
sc_adata <- FindVariableFeatures(sc_adata)
sc_adata <- ScaleData(sc_adata)
sc_adata <- RunPCA(sc_adata, npcs = 30)

st_adata <- NormalizeData(st_adata)
st_adata <- FindVariableFeatures(st_adata)
st_adata <- ScaleData(st_adata)
st_adata <- RunPCA(st_adata, npcs = 30)

traint <- CellTrek::traint(
  st_data = st_adata,
  sc_data = sc_adata,
  st_assay = "RNA",
  sc_assay = "RNA",
  cell_names = "celltype"
)

celltrek_obj <- CellTrek::celltrek(
  st_sc_int = traint,
  int_assay = "traint",
  sc_data = sc_adata,
  sc_assay = "RNA",
  reduction = "pca",
  intp = TRUE,
  intp_pnt = 5000,
  intp_lin = FALSE,
  nPCs = 30,
  ntree = 1000,
  dist_thresh = 0.55,
  top_spot = 5,
  spot_n = 50,
  repel_r = 0.0001,
  repel_iter = 20,
  keep_model = TRUE
)$celltrek

write.csv(celltrek_obj@meta.data, file = "celltrek_mousep22.csv",row.names = TRUE)


## verse
rna_expr <- read.csv("human_rna_expr.csv", row.names = 1, check.names = FALSE)
rna_meta <- read.csv("human_rna_meta.csv", row.names = 1)
atac_expr <- read.csv("human_atac_expr.csv", row.names = 1, check.names = FALSE)
atac_meta <- read.csv("human_atac_meta.csv", row.names = 1)

rownames(atac_expr) <- gsub("^atac-", "", rownames(atac_expr))
rownames(rna_expr) <- as.character(rownames(rna_expr))
common_features <- intersect(rownames(rna_expr), rownames(atac_expr))

message("shared features: ", length(common_features))

rna_expr <- rna_expr[common_features, , drop = FALSE]
atac_expr <- atac_expr[common_features, , drop = FALSE]

rna_mat <- Matrix(as.matrix(rna_expr), sparse = TRUE)
atac_mat <- Matrix(as.matrix(atac_expr), sparse = TRUE)

sc_adata <- CreateSeuratObject(
  counts = atac_mat,
  assay = "RNA"
)

sc_adata <- AddMetaData(sc_adata, atac_meta)
sc_adata$celltype <- "celltype"
sc_adata$orig.ident <- "atac"

st_adata <- CreateSeuratObject(
  counts = rna_mat,
  assay = "RNA"
)

st_adata <- AddMetaData(st_adata, atac_meta)

st_adata$orig.ident <- "RNA_spatial"
st_adata$celltype <- "celltype"

coords <- data.frame(
  tissue = 1,
  row = rna_meta$x,
  col = rna_meta$y,
  imagerow = rna_meta$x,
  imagecol = rna_meta$y
)

rownames(coords) <- rownames(rna_meta)


colnames(st_adata) <- paste0(colnames(st_adata), "-1")
rownames(rna_meta) <- colnames(st_adata)
rownames(coords) <- colnames(st_adata)

fake_image <- array(0, dim = c(10, 10, 3))

sf <- scalefactors(
  spot = 1,
  fiducial = 1,
  hires = 1,
  lowres = 1
)

st_adata@images$slice1 <- new(
  Class = "VisiumV1",
  image = fake_image,
  coordinates = coords,
  scale.factors = sf,
  assay = "RNA",
  key = "slice1_"
)

sc_adata <- NormalizeData(sc_adata)
sc_adata <- FindVariableFeatures(sc_adata)
sc_adata <- ScaleData(sc_adata)
sc_adata <- RunPCA(sc_adata, npcs = 30)

st_adata <- NormalizeData(st_adata)
st_adata <- FindVariableFeatures(st_adata)
st_adata <- ScaleData(st_adata)
st_adata <- RunPCA(st_adata, npcs = 30)

traint <- CellTrek::traint(
  st_data = st_adata,
  sc_data = sc_adata,
  st_assay = "RNA",
  sc_assay = "RNA",
  cell_names = "celltype"
)

celltrek_obj <- CellTrek::celltrek(
  st_sc_int = traint,
  int_assay = "traint",
  sc_data = sc_adata,
  sc_assay = "RNA",
  reduction = "pca",
  intp = TRUE,
  intp_pnt = 5000,
  intp_lin = FALSE,
  nPCs = 30,
  ntree = 1000,
  dist_thresh = 0.55,
  top_spot = 5,
  spot_n = 50,
  repel_r = 0.0001,
  repel_iter = 20,
  keep_model = TRUE
)$celltrek

write.csv(celltrek_obj@meta.data, file = "celltrek_human_verse.csv",row.names = TRUE)


#mouse p22
rna_expr <- read.csv("mousep22_rna_expr.csv", row.names = 1, check.names = FALSE)
rna_meta <- read.csv("mousep22_rna_meta.csv", row.names = 1)
atac_expr <- read.csv("mousep22_atac_expr.csv", row.names = 1, check.names = FALSE)
atac_meta <- read.csv("mousep22_atac_meta.csv", row.names = 1)

rownames(atac_expr) <- gsub("^atac-", "", rownames(atac_expr))
rownames(rna_expr) <- as.character(rownames(rna_expr))
common_features <- intersect(rownames(rna_expr), rownames(atac_expr))

message("shared features: ", length(common_features))

rna_expr <- rna_expr[common_features, , drop = FALSE]
atac_expr <- atac_expr[common_features, , drop = FALSE]

rna_mat <- Matrix(as.matrix(rna_expr), sparse = TRUE)
atac_mat <- Matrix(as.matrix(atac_expr), sparse = TRUE)

sc_adata <- CreateSeuratObject(
  counts = atac_mat,
  assay = "RNA"
)

sc_adata <- AddMetaData(sc_adata, atac_meta)
sc_adata$celltype <- "celltype"
sc_adata$orig.ident <- "atac"

st_adata <- CreateSeuratObject(
  counts = rna_mat,
  assay = "RNA"
)

st_adata <- AddMetaData(st_adata, atac_meta)

st_adata$orig.ident <- "RNA_spatial"
st_adata$celltype <- "celltype"

coords <- data.frame(
  tissue = 1,
  row = rna_meta$x,
  col = rna_meta$y,
  imagerow = rna_meta$x,
  imagecol = rna_meta$y
)

rownames(coords) <- rownames(rna_meta)


colnames(st_adata) <- paste0(colnames(st_adata), "-1")
rownames(rna_meta) <- colnames(st_adata)
rownames(coords) <- colnames(st_adata)

fake_image <- array(0, dim = c(10, 10, 3))

sf <- scalefactors(
  spot = 1,
  fiducial = 1,
  hires = 1,
  lowres = 1
)

st_adata@images$slice1 <- new(
  Class = "VisiumV1",
  image = fake_image,
  coordinates = coords,
  scale.factors = sf,
  assay = "RNA",
  key = "slice1_"
)

sc_adata <- NormalizeData(sc_adata)
sc_adata <- FindVariableFeatures(sc_adata)
sc_adata <- ScaleData(sc_adata)
sc_adata <- RunPCA(sc_adata, npcs = 30)

st_adata <- NormalizeData(st_adata)
st_adata <- FindVariableFeatures(st_adata)
st_adata <- ScaleData(st_adata)
st_adata <- RunPCA(st_adata, npcs = 30)

traint <- CellTrek::traint(
  st_data = st_adata,
  sc_data = sc_adata,
  st_assay = "RNA",
  sc_assay = "RNA",
  cell_names = "celltype"
)

celltrek_obj <- CellTrek::celltrek(
  st_sc_int = traint,
  int_assay = "traint",
  sc_data = sc_adata,
  sc_assay = "RNA",
  reduction = "pca",
  intp = TRUE,
  intp_pnt = 5000,
  intp_lin = FALSE,
  nPCs = 30,
  ntree = 1000,
  dist_thresh = 0.55,
  top_spot = 5,
  spot_n = 50,
  repel_r = 0.0001,
  repel_iter = 20,
  keep_model = TRUE
)$celltrek

write.csv(celltrek_obj@meta.data, file = "celltrek_mousep22_verse.csv",row.names = TRUE)


#mouse p22-h3k27ac
rna_expr <- read.csv("mousep22h3k27ac_rna_expr.csv", row.names = 1, check.names = FALSE)
rna_meta <- read.csv("mousep22h3k27ac_rna_meta.csv", row.names = 1)
atac_expr <- read.csv("mousep22h3k27ac_expr.csv", row.names = 1, check.names = FALSE)
atac_meta <- read.csv("mousep22h3k27ac_meta.csv", row.names = 1)

rownames(atac_expr) <- gsub("^atac-", "", rownames(atac_expr))
rownames(rna_expr) <- as.character(rownames(rna_expr))
common_features <- intersect(rownames(rna_expr), rownames(atac_expr))

message("shared features: ", length(common_features))

rna_expr <- rna_expr[common_features, , drop = FALSE]
atac_expr <- atac_expr[common_features, , drop = FALSE]

rna_mat <- Matrix(as.matrix(rna_expr), sparse = TRUE)

sc_adata <- CreateSeuratObject(
  counts = rna_mat,
  assay = "RNA"
)

sc_adata <- AddMetaData(sc_adata, rna_meta)
sc_adata$celltype <- "celltype"
sc_adata$orig.ident <- "scRNA"
atac_mat <- Matrix(as.matrix(atac_expr), sparse = TRUE)

st_adata <- CreateSeuratObject(
  counts = atac_mat,
  assay = "RNA"
)

st_adata <- AddMetaData(st_adata, atac_meta)

st_adata$orig.ident <- "ATAC_spatial"
st_adata$celltype <- "celltype"

coords <- data.frame(
  tissue = 1,
  row = atac_meta$x,
  col = atac_meta$y,
  imagerow = atac_meta$x,
  imagecol = atac_meta$y
)

rownames(coords) <- rownames(atac_meta)


colnames(st_adata) <- paste0(colnames(st_adata), "-1")
rownames(atac_meta) <- colnames(st_adata)
rownames(coords) <- colnames(st_adata)

fake_image <- array(0, dim = c(10, 10, 3))

sf <- scalefactors(
  spot = 1,
  fiducial = 1,
  hires = 1,
  lowres = 1
)

st_adata@images$slice1 <- new(
  Class = "VisiumV1",
  image = fake_image,
  coordinates = coords,
  scale.factors = sf,
  assay = "RNA",
  key = "slice1_"
)

sc_adata <- NormalizeData(sc_adata)
sc_adata <- FindVariableFeatures(sc_adata)
sc_adata <- ScaleData(sc_adata)
sc_adata <- RunPCA(sc_adata, npcs = 30)

st_adata <- NormalizeData(st_adata)
st_adata <- FindVariableFeatures(st_adata)
st_adata <- ScaleData(st_adata)
st_adata <- RunPCA(st_adata, npcs = 30)

traint <- CellTrek::traint(
  st_data = st_adata,
  sc_data = sc_adata,
  st_assay = "RNA",
  sc_assay = "RNA",
  cell_names = "celltype"
)

celltrek_obj <- CellTrek::celltrek(
  st_sc_int = traint,
  int_assay = "traint",
  sc_data = sc_adata,
  sc_assay = "RNA",
  reduction = "pca",
  intp = TRUE,
  intp_pnt = 5000,
  intp_lin = FALSE,
  nPCs = 30,
  ntree = 1000,
  dist_thresh = 0.55,
  top_spot = 5,
  spot_n = 50,
  repel_r = 0.0001,
  repel_iter = 20,
  keep_model = TRUE
)$celltrek

write.csv(celltrek_obj@meta.data, file = "celltrek_mousep22h3k27ac.csv",row.names = TRUE)

#mouse p22-h3k27ac
rna_expr <- read.csv("mousep22h3k27ac_rna_expr.csv", row.names = 1, check.names = FALSE)
rna_meta <- read.csv("mousep22h3k27ac_rna_meta.csv", row.names = 1)
atac_expr <- read.csv("mousep22h3k27ac_expr.csv", row.names = 1, check.names = FALSE)
atac_meta <- read.csv("mousep22h3k27ac_meta.csv", row.names = 1)

rownames(atac_expr) <- gsub("^atac-", "", rownames(atac_expr))
rownames(rna_expr) <- as.character(rownames(rna_expr))
common_features <- intersect(rownames(rna_expr), rownames(atac_expr))

message("shared features: ", length(common_features))

rna_expr <- rna_expr[common_features, , drop = FALSE]
atac_expr <- atac_expr[common_features, , drop = FALSE]

rna_mat <- Matrix(as.matrix(rna_expr), sparse = TRUE)
atac_mat <- Matrix(as.matrix(atac_expr), sparse = TRUE)

sc_adata <- CreateSeuratObject(
  counts = atac_mat,
  assay = "RNA"
)

sc_adata <- AddMetaData(sc_adata, atac_meta)
sc_adata$celltype <- "celltype"
sc_adata$orig.ident <- "atac"

st_adata <- CreateSeuratObject(
  counts = rna_mat,
  assay = "RNA"
)

st_adata <- AddMetaData(st_adata, atac_meta)

st_adata$orig.ident <- "RNA_spatial"
st_adata$celltype <- "celltype"

coords <- data.frame(
  tissue = 1,
  row = rna_meta$x,
  col = rna_meta$y,
  imagerow = rna_meta$x,
  imagecol = rna_meta$y
)

rownames(coords) <- rownames(rna_meta)


colnames(st_adata) <- paste0(colnames(st_adata), "-1")
rownames(rna_meta) <- colnames(st_adata)
rownames(coords) <- colnames(st_adata)

fake_image <- array(0, dim = c(10, 10, 3))

sf <- scalefactors(
  spot = 1,
  fiducial = 1,
  hires = 1,
  lowres = 1
)

st_adata@images$slice1 <- new(
  Class = "VisiumV1",
  image = fake_image,
  coordinates = coords,
  scale.factors = sf,
  assay = "RNA",
  key = "slice1_"
)

sc_adata <- NormalizeData(sc_adata)
sc_adata <- FindVariableFeatures(sc_adata)
sc_adata <- ScaleData(sc_adata)
sc_adata <- RunPCA(sc_adata, npcs = 30)

st_adata <- NormalizeData(st_adata)
st_adata <- FindVariableFeatures(st_adata)
st_adata <- ScaleData(st_adata)
st_adata <- RunPCA(st_adata, npcs = 30)

traint <- CellTrek::traint(
  st_data = st_adata,
  sc_data = sc_adata,
  st_assay = "RNA",
  sc_assay = "RNA",
  cell_names = "celltype"
)

celltrek_obj <- CellTrek::celltrek(
  st_sc_int = traint,
  int_assay = "traint",
  sc_data = sc_adata,
  sc_assay = "RNA",
  reduction = "pca",
  intp = TRUE,
  intp_pnt = 5000,
  intp_lin = FALSE,
  nPCs = 30,
  ntree = 1000,
  dist_thresh = 0.55,
  top_spot = 5,
  spot_n = 50,
  repel_r = 0.0001,
  repel_iter = 20,
  keep_model = TRUE
)$celltrek

write.csv(celltrek_obj@meta.data, file = "celltrek_mousep22h3k27ac_verse.csv",row.names = TRUE)

