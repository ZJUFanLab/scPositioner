options(stringsAsFactors = F)
library("CellTrek")
library("dplyr")
library("Seurat")
library("viridis")
library("ConsensusClusterPlus")
library(SeuratDisk)


######################## preprocess ########################
Convert(paste0('workspace/scPositioner/tissue_structure/Kidney/processed/sc_adata.h5ad'), "h5seurat",overwrite = TRUE,assay = "RNA")
sc_adata <- LoadH5Seurat('workspace/scPositioner/tissue_structure/Kidney/processed/sc_adata.h5seurat', meta.data = FALSE, misc = FALSE)
sc_obs <- read.csv('workspace/scPositioner/tissue_structure/Kidney/processed/sc_obs.csv', row.names = 1)
sc_adata <- AddMetaData(sc_adata, sc_obs)

st_adata <- Load10X_Spatial('workspace/scPositioner/tissue_structure/Kidney/ST/')
colnames(st_adata@assays$Spatial@counts) <- gsub('-', '\\.', colnames(st_adata@assays$Spatial@counts))
colnames(st_adata@assays$Spatial@data) <- gsub('-', '\\.', colnames(st_adata@assays$Spatial@data))
rownames(st_adata@meta.data) <- gsub('-', '\\.', rownames(st_adata@meta.data) )
rownames(st_adata@images[["slice1"]]@coordinates) <- gsub('-', '\\.', rownames(st_adata@images[["slice1"]]@coordinates)  )
names(st_adata@active.ident) <- gsub('-', '\\.', names(st_adata@active.ident))

common_genes <- intersect(rownames(sc_adata), rownames(st_adata))
sc_adata <- sc_adata[common_genes, ]
st_adata <- st_adata[common_genes, ]

SpatialDimPlot(st_adata)
######################## run ########################
t1 <- Sys.time()
traint <- CellTrek::traint(st_data=st_adata, sc_data=sc_adata, sc_assay='RNA', cell_names='Cell_type')
DimPlot(traint, group.by = "type") 
spot_n = 5
repel_r = 0.0001

celltrek_obj <- CellTrek::celltrek(
  st_sc_int=traint, 
  int_assay='traint',
  sc_data=sc_adata, 
  sc_assay = 'RNA', 
  reduction='pca', 
  intp=T, 
  intp_pnt=5000, 
  intp_lin=F, 
  nPCs=30, 
  ntree=1000, 
  dist_thresh=0.55, 
  top_spot=5, 
  spot_n=spot_n, 
  repel_r=repel_r, 
  repel_iter=20, 
  keep_model=T)$celltrek

t2 <- Sys.time()
t2 - t1

saveRDS(celltrek_obj, file = 'workspace/scPositioner/tissue_structure/Kidney/results/CellTrek_result.rds')

SaveH5Seurat(celltrek_obj, filename = "workspace/scPositioner/tissue_structure/Kidney/results/CellTrek_result.h5Seurat")
Convert("workspace/scPositioner/tissue_structure/Kidney/results/CellTrek_result.h5Seurat", dest = "h5ad")
