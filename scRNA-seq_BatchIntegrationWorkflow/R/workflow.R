# R/workflow.R
make_seurat_one_batch <- function(counts, batch_name) {
  obj <- Seurat::CreateSeuratObject(counts = counts, project = batch_name)
  obj$batch <- batch_name
  obj[["percent.mt"]] <- Seurat::PercentageFeatureSet(obj, pattern = "^MT-")
  obj
}

apply_qc <- function(obj, min_genes, max_genes, max_mt) {
  subset(
    obj,
    subset = nFeature_RNA >= min_genes &
      nFeature_RNA <= max_genes &
      percent.mt <= max_mt
  )
}


preprocess_seurat <- function(obj, nfeatures, npcs) {
  obj <- Seurat::NormalizeData(obj)
  obj <- Seurat::FindVariableFeatures(obj, selection.method = "vst", nfeatures = nfeatures)
  obj <- Seurat::ScaleData(obj, verbose = FALSE)
  obj <- Seurat::RunPCA(obj, npcs = npcs, verbose = FALSE)
  obj <- Seurat::FindNeighbors(obj, dims = 1:min(npcs, 30))
  obj <- Seurat::FindClusters(obj, resolution = 0.5)
  obj <- Seurat::RunUMAP(obj, dims = 1:min(npcs, 30))
  obj
}
