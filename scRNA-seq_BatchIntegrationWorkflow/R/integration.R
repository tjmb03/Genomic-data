# R/integration.R
has_harmony <- function() requireNamespace("harmony", quietly = TRUE)
has_batchelor <- function() requireNamespace("batchelor", quietly = TRUE) &&
  requireNamespace("SingleCellExperiment", quietly = TRUE) &&
  requireNamespace("scater", quietly = TRUE)

integrate_seurat_anchors <- function(objs, nfeatures = 3000, npcs = 30) {
  features <- Seurat::SelectIntegrationFeatures(object.list = objs, nfeatures = nfeatures)

  objs2 <- lapply(objs, function(o) {
    o <- Seurat::ScaleData(o, features = features, verbose = FALSE)
    o <- Seurat::RunPCA(o, features = features, npcs = npcs, verbose = FALSE)
    o
  })

  anchors <- Seurat::FindIntegrationAnchors(
    object.list = objs2,
    anchor.features = features,
    reduction = "rpca"
  )

  obj.int <- Seurat::IntegrateData(anchorset = anchors)
  Seurat::DefaultAssay(obj.int) <- "integrated"

  obj.int <- Seurat::ScaleData(obj.int, verbose = FALSE)
  obj.int <- Seurat::RunPCA(obj.int, npcs = npcs, verbose = FALSE)
  obj.int <- Seurat::FindNeighbors(obj.int, dims = 1:min(npcs, 30))
  obj.int <- Seurat::FindClusters(obj.int, resolution = 0.5)
  obj.int <- Seurat::RunUMAP(obj.int, dims = 1:min(npcs, 30))
  obj.int
}

integrate_harmony <- function(obj.merge, group_var = "batch", npcs = 30) {
  if (!has_harmony()) stop("Harmony not installed. Run: install.packages('harmony')")
  obj.merge <- Seurat::ScaleData(obj.merge, verbose = FALSE)
  obj.merge <- Seurat::RunPCA(obj.merge, npcs = npcs, verbose = FALSE)
  obj.merge <- harmony::RunHarmony(obj.merge, group.by.vars = group_var)
  # Keep both harmony and UMAP for visualization
  obj.merge <- Seurat::RunUMAP(obj.merge, reduction = "harmony", dims = 1:min(npcs, 30))
  obj.merge <- Seurat::FindNeighbors(obj.merge, reduction = "harmony", dims = 1:min(npcs, 30))
  obj.merge <- Seurat::FindClusters(obj.merge, resolution = 0.5)
  obj.merge
}

integrate_fastmnn <- function(objs, npcs = 30, k = 20) {
  if (!has_batchelor()) {
    stop("fastMNN requires Bioconductor packages: batchelor, SingleCellExperiment, scater. See README.")
  }

  # Convert each batch to SingleCellExperiment with log-normalized counts
  sce_list <- lapply(names(objs), function(nm) {
    o <- objs[[nm]]
    sce <- Seurat::as.SingleCellExperiment(o)
    sce$batch <- nm

    # BUG FIX: 10x gene symbol tables frequently contain duplicate names (alias
    # genes, pseudogene loci).  fastMNN propagates them into "reconstructed", and
    # CreateSeuratObject -> LogMap then rejects the matrix with
    # "Duplicate rownames not allowed".  Deduplicate here, before fastMNN runs,
    # so the correction itself also sees clean gene names.
    if (anyDuplicated(rownames(sce))) {
      rownames(sce) <- make.unique(rownames(sce))
    }
    sce
  })
  names(sce_list) <- names(objs)

  sce_list <- lapply(sce_list, function(sce) {
    sce <- scater::logNormCounts(sce)
    sce
  })

  # fastMNN batch correction; returns a corrected SCE with reducedDim "corrected"
  # BUG FIX: fastMNN takes batches as individual ... args, not a list.
  # Passing sce_list directly treated all cells as a single batch (no correction).
  corrected <- do.call(
    batchelor::fastMNN,
    c(sce_list, list(k = k, d = npcs, BSPARAM = BiocSingular::ExactParam()))
  )

  # Build a Seurat object for convenience: use first batch counts (not used for corrected embedding)
  # We'll store the corrected embedding as a reduction "mnn"
  # Also attach metadata (including batch)
  # BUG FIX: batchelor::fastMNN output SCE contains assay "reconstructed"
  # (batch-corrected log-expression), NOT "counts". Using "counts" throws an error.
  # We use "reconstructed" as a proxy expression matrix for the Seurat shell object;
  # all downstream analysis uses the corrected "mnn" embedding, not these values.
  mat_counts <- as.matrix(SummarizedExperiment::assay(corrected, "reconstructed"))

  # Deduplicate gene names (rownames) — aliases / pseudogene loci.
  if (anyDuplicated(rownames(mat_counts))) {
    rownames(mat_counts) <- make.unique(rownames(mat_counts))
  }

  # BUG FIX: the actual crash — LogMap(y = cells.all) fails with
  # "Duplicate rownames not allowed" because its "rownames" ARE the cell barcodes.
  # fastMNN concatenates SCEs from all batches without adding any batch prefix,
  # so identical barcodes from different batches (e.g. "AAACCTGA-1" in both
  # BatchA and BatchB) collide in colnames(mat_counts).
  # Fix: stamp each cell name with its batch label before building the Seurat object,
  # mirroring what Seurat::merge(..., add.cell.ids = ...) does for merged objects.
  batch_labels <- corrected$batch
  if (anyDuplicated(colnames(mat_counts))) {
    colnames(mat_counts) <- paste(batch_labels, colnames(mat_counts), sep = "_")
  }

  obj <- Seurat::CreateSeuratObject(mat_counts, project = "fastMNN")
  obj$batch <- batch_labels

  # Add corrected embedding
  emb <- SingleCellExperiment::reducedDim(corrected, "corrected")
  rownames(emb) <- colnames(obj)
  # BUG FIX: DefaultAssay() was unqualified; use Seurat:: to be safe in all call contexts.
  obj[["mnn"]] <- Seurat::CreateDimReducObject(embeddings = emb, key = "MNN_", assay = Seurat::DefaultAssay(obj))

  # Graph + clusters + UMAP using corrected embedding
  obj <- Seurat::FindNeighbors(obj, reduction = "mnn", dims = 1:min(npcs, ncol(emb)))
  obj <- Seurat::FindClusters(obj, resolution = 0.5)
  obj <- Seurat::RunUMAP(obj, reduction = "mnn", dims = 1:min(npcs, ncol(emb)))
  obj
}
