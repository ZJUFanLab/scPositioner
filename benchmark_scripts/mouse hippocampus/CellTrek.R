options(stringsAsFactors = F)
library("CellTrek")
library("dplyr")
library("Seurat")
library("viridis")
library("ConsensusClusterPlus")
library(SeuratDisk)


######################## Cerebellum ########################
######################## preprocess ########################
noise_list <- c('0', '05', '10', '20', '40')
n_list <- c(5, 10, 15)

for (noise in noise_list) {
  for (n in n_list) {
    # sc
    Convert(
      paste0('workspace/scPositioner/benchmark/Cerebellum-Puck_180430_6/output/Cerebellum_sc_noise', noise, '.h5ad'), 
      "h5seurat",overwrite = TRUE,assay = "RNA"
      )
    sc_adata <- LoadH5Seurat(
      paste0('workspace/scPositioner/benchmark/Cerebellum-Puck_180430_6/output/Cerebellum_sc_noise', noise, '.h5seurat'), 
      meta.data = FALSE, misc = FALSE
      )
    
    sc_obs <- read.csv('workspace/scPositioner/benchmark/Cerebellum-Puck_180430_6/output/sc_obs.csv', row.names = 1)
    sc_adata <- AddMetaData(sc_adata, sc_obs)
    colnames(sc_adata@assays$RNA@counts) <- gsub('-', '\\.', colnames(sc_adata@assays$RNA@counts))
    colnames(sc_adata@assays$RNA@data) <- gsub('-', '\\.', colnames(sc_adata@assays$RNA@data))
    colnames(sc_adata@assays$RNA@scale.data) <- gsub('-', '\\.', colnames(sc_adata@assays$RNA@scale.data))
    names(sc_adata@active.ident) <- gsub('-', '\\.', names(sc_adata@active.ident))
    sc_adata$orig.ident <- Idents(sc_adata)
    
    # st
    Convert(
      paste0('workspace/scPositioner/benchmark/Cerebellum-Puck_180430_6/output/Cerebellum_st_n', n, '.h5ad'), 
      "h5seurat",overwrite = TRUE,assay = "RNA"
      )
    
    st_adata <- LoadH5Seurat(
      paste0('workspace/scPositioner/benchmark/Cerebellum-Puck_180430_6/output/Cerebellum_st_n', n, '.h5seurat'), 
      meta.data = FALSE, misc = FALSE
      )
    
    sfs <- scalefactors(spot = 1, fiducial = 1, hires = 1, lowres = 1)
    images <- Read10X_Image('workspace/scPositioner/tumor2/HGSOC/ST-GSE211956/spatial/')
    coords <- data.frame(
      tissue = 1,
      row = st_adata@reductions$spatial@cell.embeddings[, 1],
      col = st_adata@reductions$spatial@cell.embeddings[, 2],
      imagerow = st_adata@reductions$spatial@cell.embeddings[, 1],
      imagecol = st_adata@reductions$spatial@cell.embeddings[, 2]
    )
    rownames(coords) <- rownames(st_adata@reductions$spatial@cell.embeddings)
    
    st_adata@images$slice1 =  new(
      image = images@image,
      Class = 'VisiumV1',
      assay = "Spatial",
      key = "slice1_",
      coordinates = coords,
      scale.factors = sfs
    )
    st_adata$orig.ident <- Idents(st_adata)
    
    # run 
    traint <- CellTrek::traint(st_data=st_adata, sc_data=sc_adata, st_assay="RNA", sc_assay='RNA', cell_names='CellType')
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
      spot_n=n, 
      repel_r=repel_r, 
      repel_iter=20, 
      keep_model=T)$celltrek
    
    save_name <- paste0('Cerebellum_CellTrek_noise', noise, '_n', n)
    
    saveRDS(celltrek_obj, file = paste0('workspace/scPositioner/benchmark/Cerebellum-Puck_180430_6/results/', save_name, '.rds'))
    SaveH5Seurat(celltrek_obj, filename = paste0('workspace/scPositioner/benchmark/Cerebellum-Puck_180430_6/results/', save_name, '.h5Seurat'))
    Convert(paste0('workspace/scPositioner/benchmark/Cerebellum-Puck_180430_6/results/', save_name, '.h5Seurat'), dest = "h5ad")
    
  }
}


######################## Hippocampus ########################
######################## preprocess ########################
noise_list <- c('0', '05', '10', '20', '40')
n_list <- c(5, 10, 15)

for (noise in noise_list) {
  for (n in n_list) {
    # sc
    Convert(
      paste0('workspace/scPositioner/benchmark/Hippocampus-Puck_180611_2/output/Hippocampus_sc_noise', noise, '.h5ad'), 
      "h5seurat",overwrite = TRUE,assay = "RNA"
    )
    sc_adata <- LoadH5Seurat(
      paste0('workspace/scPositioner/benchmark/Hippocampus-Puck_180611_2/output/Hippocampus_sc_noise', noise, '.h5seurat'), 
      meta.data = FALSE, misc = FALSE
    )
    
    sc_obs <- read.csv('workspace/scPositioner/benchmark/Hippocampus-Puck_180611_2/output/sc_obs.csv', row.names = 1)
    sc_adata <- AddMetaData(sc_adata, sc_obs)
    colnames(sc_adata@assays$RNA@counts) <- gsub('-', '\\.', colnames(sc_adata@assays$RNA@counts))
    colnames(sc_adata@assays$RNA@data) <- gsub('-', '\\.', colnames(sc_adata@assays$RNA@data))
    colnames(sc_adata@assays$RNA@scale.data) <- gsub('-', '\\.', colnames(sc_adata@assays$RNA@scale.data))
    names(sc_adata@active.ident) <- gsub('-', '\\.', names(sc_adata@active.ident))
    sc_adata$orig.ident <- Idents(sc_adata)
    
    # st
    Convert(
      paste0('workspace/scPositioner/benchmark/Hippocampus-Puck_180611_2/output/Hippocampus_st_n', n, '.h5ad'), 
      "h5seurat",overwrite = TRUE,assay = "RNA"
    )
    
    st_adata <- LoadH5Seurat(
      paste0('workspace/scPositioner/benchmark/Hippocampus-Puck_180611_2/output/Hippocampus_st_n', n, '.h5seurat'), 
      meta.data = FALSE, misc = FALSE
    )
    
    sfs <- scalefactors(spot = 1, fiducial = 1, hires = 1, lowres = 1)
    images <- Read10X_Image('workspace/scPositioner/tumor2/HGSOC/ST-GSE211956/spatial/')
    coords <- data.frame(
      tissue = 1,
      row = st_adata@reductions$spatial@cell.embeddings[, 1],
      col = st_adata@reductions$spatial@cell.embeddings[, 2],
      imagerow = st_adata@reductions$spatial@cell.embeddings[, 1],
      imagecol = st_adata@reductions$spatial@cell.embeddings[, 2]
    )
    rownames(coords) <- rownames(st_adata@reductions$spatial@cell.embeddings)
    
    st_adata@images$slice1 =  new(
      image = images@image,
      Class = 'VisiumV1',
      assay = "Spatial",
      key = "slice1_",
      coordinates = coords,
      scale.factors = sfs
    )
    st_adata$orig.ident <- Idents(st_adata)
    
    # run 
    traint <- CellTrek::traint(st_data=st_adata, sc_data=sc_adata, st_assay="RNA", sc_assay='RNA', cell_names='CellType')
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
      spot_n=n, 
      repel_r=repel_r, 
      repel_iter=20, 
      keep_model=T)$celltrek
    
    save_name <- paste0('Hippocampus_CellTrek_noise', noise, '_n', n)
    
    saveRDS(celltrek_obj, file = paste0('workspace/scPositioner/benchmark/Hippocampus-Puck_180611_2/results/', save_name, '.rds'))
    SaveH5Seurat(celltrek_obj, filename = paste0('workspace/scPositioner/benchmark/Hippocampus-Puck_180611_2/results/', save_name, '.h5Seurat'))
    Convert(paste0('workspace/scPositioner/benchmark/Hippocampus-Puck_180611_2/results/', save_name, '.h5Seurat'), dest = "h5ad")
    
  }
}
