# R/diagnostics.R
# Metrics:
# - Silhouette mean (batch, clusters)
# - Variance explained (R^2) of batch vs clusters on embedding (ANOVA-style)
#
# Computed on a subsample for speed.

.subsample_cells <- function(obj, max_cells = 20000) {
  n <- ncol(obj)
  if (n <= max_cells) return(colnames(obj))
  sample(colnames(obj), size = max_cells, replace = FALSE)
}

.choose_reduction <- function(obj, method_name) {
  # Prefer method-specific latent space where applicable
  if (method_name == "harmony" && "harmony" %in% names(obj@reductions)) return("harmony")
  if (method_name == "fastmnn" && "mnn" %in% names(obj@reductions)) return("mnn")
  # seurat integrated + merged typically have PCA
  if ("pca" %in% names(obj@reductions)) return("pca")
  if ("umap" %in% names(obj@reductions)) return("umap")
  return(names(obj@reductions)[1])
}

.get_embedding <- function(obj, reduction) {
  Seurat::Embeddings(obj, reduction)
}

compute_silhouette_mean <- function(emb, labels) {
  if (length(unique(labels)) < 2) return(NA_real_)
  if (nrow(emb) < 3) return(NA_real_)
  if (!requireNamespace("cluster", quietly = TRUE)) stop("Install 'cluster' for silhouette metrics.")
  labs <- as.integer(as.factor(labels))
  d <- stats::dist(emb)
  sil <- cluster::silhouette(labs, d)
  mean(sil[, "sil_width"], na.rm = TRUE)
}

r2_oneway <- function(y, group) {
  # One-way ANOVA R^2: 1 - SSE/SST
  g <- as.factor(group)
  if (nlevels(g) < 2) return(NA_real_)
  y <- as.numeric(y)
  mu <- mean(y, na.rm = TRUE)
  sst <- sum((y - mu)^2, na.rm = TRUE)
  if (sst == 0) return(0)
  # within-group SSE
  sse <- 0
  for (lv in levels(g)) {
    idx <- which(g == lv)
    yi <- y[idx]
    mi <- mean(yi, na.rm = TRUE)
    sse <- sse + sum((yi - mi)^2, na.rm = TRUE)
  }
  1 - (sse / sst)
}

compute_r2_mean <- function(emb, group, weights = NULL) {
  # Compute mean R^2 across embedding dimensions, optionally weighted
  k <- ncol(emb)
  r2s <- vapply(seq_len(k), function(j) r2_oneway(emb[, j], group), numeric(1))
  if (is.null(weights)) {
    return(mean(r2s, na.rm = TRUE))
  } else {
    w <- weights[seq_len(min(length(weights), k))]
    r2s <- r2s[seq_along(w)]
    if (all(is.na(r2s))) return(NA_real_)
    sum(w * r2s, na.rm = TRUE) / sum(w[!is.na(r2s)])
  }
}

compute_method_metrics <- function(obj, method_name, max_cells = 20000, dims_use = 30) {
  if (is.null(obj)) return(NULL)

  keep <- .subsample_cells(obj, max_cells = max_cells)
  meta <- obj@meta.data[keep, , drop = FALSE]

  if (!("batch" %in% colnames(meta))) meta$batch <- "batch1"
  if (!("seurat_clusters" %in% colnames(meta))) meta$seurat_clusters <- "0"

  red <- .choose_reduction(obj, method_name)
  emb <- .get_embedding(obj, red)
  emb <- emb[keep, , drop = FALSE]
  emb <- emb[, 1:min(dims_use, ncol(emb)), drop = FALSE]

  # weights: for PCA use stdev^2 if available
  weights <- NULL
  if (red == "pca" && !is.null(obj@reductions$pca@stdev)) {
    v <- obj@reductions$pca@stdev^2
    weights <- v / sum(v)
  }

  data.frame(
    method = method_name,
    reduction_used = red,
    n_cells_total = ncol(obj),
    n_cells_metrics = length(keep),
    dims_used = ncol(emb),
    batch_sil_mean = compute_silhouette_mean(emb, meta$batch),
    cluster_sil_mean = compute_silhouette_mean(emb, meta$seurat_clusters),
    batch_r2_mean = compute_r2_mean(emb, meta$batch, weights = weights),
    cluster_r2_mean = compute_r2_mean(emb, meta$seurat_clusters, weights = weights),
    stringsAsFactors = FALSE
  )
}

compute_all_metrics <- function(candidates, max_cells = 20000, dims_use = 30) {
  rows <- list()
  for (nm in names(candidates)) {
    obj <- candidates[[nm]]
    if (is.null(obj)) next
    rows[[nm]] <- compute_method_metrics(obj, nm, max_cells = max_cells, dims_use = dims_use)
  }
  # BUG FIX: variable was named 'summary', shadowing base::summary(). Renamed to metrics_df.
  metrics_df <- do.call(rbind, rows)
  list(summary = metrics_df)
}

decide_best_method <- function(metrics) {
  df <- metrics$summary
  if (nrow(df) == 0) return(data.frame(note = "No metrics available."))

  norm01 <- function(x, invert = FALSE) {
    if (all(is.na(x))) return(rep(NA_real_, length(x)))
    rng <- range(x, na.rm = TRUE)
    if (diff(rng) == 0) return(rep(0.5, length(x)))
    y <- (x - rng[1]) / diff(rng)
    if (invert) y <- 1 - y
    y
  }

  # Higher score = better
  s_batch_sil <- norm01(df$batch_sil_mean, invert = TRUE)      # lower better
  s_batch_r2  <- norm01(df$batch_r2_mean, invert = TRUE)       # lower better
  s_clu_sil   <- norm01(df$cluster_sil_mean, invert = FALSE)   # higher better
  s_clu_r2    <- norm01(df$cluster_r2_mean, invert = FALSE)    # higher better

  # Weights: mixing 50%, conservation 50%, each split across silhouette/R2
  score <- 0.25*s_batch_sil + 0.25*s_batch_r2 + 0.25*s_clu_sil + 0.25*s_clu_r2

  out <- df
  out$score <- score
  out$recommended <- ""
  best_idx <- which.max(out$score)
  if (length(best_idx) == 1 && !is.na(best_idx)) out$recommended[best_idx] <- "✅ RECOMMENDED"

  # Overcorrection heuristic vs merged
  out$overcorrection_flag <- ""
  ref <- out[out$method == "merged", , drop = FALSE]
  if (nrow(ref) == 1) {
    for (i in seq_len(nrow(out))) {
      if (out$method[i] == "merged") next
      if (is.na(out$batch_sil_mean[i]) || is.na(out$cluster_sil_mean[i])) next
      batch_improve <- ref$batch_sil_mean - out$batch_sil_mean[i]
      cluster_drop <- ref$cluster_sil_mean - out$cluster_sil_mean[i]
      if (batch_improve > 0.05 && cluster_drop > 0.05) out$overcorrection_flag[i] <- "⚠ possible overcorrection"
    }
  }

  out[order(-out$score), ]
}

plot_metric_bars <- function(summary_df, metric_col, title) {
  df <- summary_df
  df$method <- factor(df$method, levels = df$method)
  ggplot(df, aes(x = method, y = .data[[metric_col]])) +
    geom_col() +
    labs(title = title, x = "", y = metric_col) +
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 20, hjust = 1))
}
