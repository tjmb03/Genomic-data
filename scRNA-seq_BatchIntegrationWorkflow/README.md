# scRNA-seq Batch Integration Workflow


A reproducible, publication-ready R pipeline for benchmarking and selecting the best batch correction method for 10x scRNA-seq data. The workflow runs three integration methods head-to-head, scores them on quantitative metrics, and auto-recommends the winner — all inside a single rendered HTML report.

🔗 Live App:
 https://tjmb03.github.io/Genomic-data/

---

## Overview

```
10x input(s)  ──►  QC & preprocessing  ──►  3 integration methods
  .h5 / MTX /                                 │
  archive                                      ├── Seurat Anchors (rPCA)
                                               ├── Harmony
                                               └── fastMNN (batchelor)
                                                        │
                                               Benchmarking metrics
                                               (silhouette + R²)
                                                        │
                                               Decision table + Figures
                                                        │
                                               outputs/*.rds / *.csv
```

---

## Features

- **Flexible input** — accepts `.h5`, MTX directories, or archives (`.tar.gz`, `.tgz`, `.tar`, `.zip`) containing a valid 10x directory; auto-extracted at runtime
- **Three correction methods** — Seurat rPCA anchors, Harmony, fastMNN — each individually toggleable via report params
- **Dual-axis benchmarking** — every method is scored on both batch mixing (lower is better) and biology conservation (higher is better) using both silhouette width and ANOVA-style R², computed on a fast subsample
- **Overcorrection detection** — flags methods where batch mixing improved but cluster structure degraded relative to the unintegrated baseline
- **Publication figures** — three patchwork figure panels (QC, UMAP comparison, metric bar charts) ready for manuscripts
- **Saved artifacts** — `.rds` objects and `.csv` tables written to `outputs/` for downstream use

---

## Project Structure

```
scRNA-integration/
├── integration_report.Rmd   # Entry point — render this
├── R/
│   ├── io_10x.R             # Input loading (.h5, MTX dir, archives)
│   ├── workflow.R           # Per-batch QC, normalisation, PCA
│   ├── integration.R        # Seurat anchors, Harmony, fastMNN
│   ├── diagnostics.R        # Silhouette, R², scoring, best-method logic
│   └── reporting_tables.R   # Benchmark table with grade labels
├── data/
│   ├── batchA/              # .h5, MTX directory, or archive
│   └── batchB/
└── outputs/                 # Generated on render
    ├── integration_figure1_qc_baseline.png
    ├── integration_figure2_umap_comparison.png
    ├── integration_figure3_benchmarking.png
    ├── integration_metrics.csv
    └── integration_decision.csv
```

---

## Installation

**CRAN packages:**
```r
install.packages(c("Seurat", "Matrix", "ggplot2", "dplyr",
                   "cluster", "rmarkdown", "patchwork", "harmony"))
```

**Bioconductor packages** (required for fastMNN):
```r
if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")
BiocManager::install(c("batchelor", "SingleCellExperiment",
                       "scater", "BiocSingular", "SummarizedExperiment"))
```

---

## Usage

Place your 10x data under `data/` (any supported format), then render from the project root:

```r
rmarkdown::render(
  "integration_report.Rmd",
  params = list(
    tenx_h5      = c("data/batchA/filtered_feature_bc_matrix.h5",
                     "data/batchB/filtered_feature_bc_matrix.h5"),
    batch_names  = c("BatchA", "BatchB"),
    max_cells_metrics = 20000
  )
)
```

The report is **location-independent** — paths are resolved relative to the Rmd file, so it renders correctly regardless of your R session's working directory.

---

## Parameters

| Parameter | Default | Description |
|---|---|---|
| `tenx_h5` | — | Vector of input paths (`.h5`, MTX dir, or archive) |
| `batch_names` | — | Labels for each batch; auto-derived from paths if omitted |
| `min_genes` | `300` | Min features per cell (QC filter) |
| `max_genes` | `6000` | Max features per cell (QC filter) |
| `max_percent_mt` | `20` | Max mitochondrial % (QC filter) |
| `nfeatures` | `3000` | Variable features for HVG selection |
| `npcs` | `30` | Number of PCs / latent dimensions |
| `max_cells_metrics` | `20000` | Subsample size for benchmarking |
| `run_seurat_integration` | `true` | Toggle Seurat anchors method |
| `run_harmony` | `true` | Toggle Harmony method |
| `run_fastmnn` | `true` | Toggle fastMNN method |
| `fastmnn_k` | `20` | k-nearest neighbours for fastMNN |
| `save_png` | `true` | Save figures 1–3 as PNG files to `outputs/` |
| `png_width` | `14` | PNG width in inches |
| `png_height` | `10` | PNG height in inches |
| `png_dpi` | `150` | PNG resolution (use 300 for print) |
| `out_prefix` | `outputs/integration` | Prefix for all output files |

---

## Benchmarking Metrics

Each method is evaluated on four scores, then ranked by a composite:

| Metric | Axis | Direction |
|---|---|---|
| Batch silhouette mean | Batch mixing | Lower → better mixing |
| Batch R² (ANOVA-style) | Batch mixing | Lower → less batch variance |
| Cluster silhouette mean | Biology conservation | Higher → tighter clusters |
| Cluster R² (ANOVA-style) | Biology conservation | Higher → more cluster variance |

Scores are min-max normalised across methods and averaged with equal weight (25% each). An overcorrection flag is raised when a method strongly reduces batch silhouette but also degrades cluster silhouette relative to the unintegrated baseline.

---

## Output Report Sections

- **Figure 1** — pre-filter QC distributions (violin + scatter) and baseline UMAP with no correction
- **Figure 2** — side-by-side UMAP panels for every method, coloured by batch and by cluster
- **Figure 3** — metric bar charts for all four scores across methods
- **Decision table** — full numeric scores, grade labels (Excellent / Good / Moderate / Poor), runtimes, recommended method, and overcorrection flags

---

## Supported Input Formats

| Format | Example |
|---|---|
| 10x HDF5 | `filtered_feature_bc_matrix.h5` |
| MTX directory | folder with `matrix.mtx.gz`, `barcodes.tsv.gz`, `features.tsv.gz` |
| Gzipped tar | `pbmc_filtered.tar.gz` containing an MTX directory |
| Zip archive | `pbmc_filtered.zip` containing an MTX directory |

---

## Dependencies

| Package | Source | Role |
|---|---|---|
| Seurat ≥ 5 | CRAN | Core scRNA-seq object + Seurat anchors |
| harmony | CRAN | Harmony integration |
| batchelor | Bioconductor | fastMNN integration |
| SingleCellExperiment | Bioconductor | SCE conversion for fastMNN |
| scater | Bioconductor | Log-normalisation for SCE |
| BiocSingular | Bioconductor | Exact SVD backend for fastMNN |
| SummarizedExperiment | Bioconductor | Assay access on fastMNN output |
| patchwork | CRAN | Figure layout |
| cluster | CRAN | Silhouette computation |
| ggplot2, dplyr | CRAN | Plotting and data wrangling |


