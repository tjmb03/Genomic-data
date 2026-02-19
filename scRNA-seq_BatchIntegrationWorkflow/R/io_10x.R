# R/io_10x.R
# Robust 10x input loader for local (non-Shiny) workflows.
# Supports:
#   - .h5 (10x feature-barcode matrix) via Seurat::Read10X_h5()
#   - directory containing matrix.mtx(.gz) + barcodes.tsv(.gz) + features.tsv(.gz)/genes.tsv(.gz) via Seurat::Read10X()
#   - archives: .tar.gz/.tgz/.tar/.zip that contain a 10x directory (same trio), auto-extracted to a temp dir

read_10x_h5 <- function(path) {
  stopifnot(file.exists(path))
  counts <- Seurat::Read10X_h5(path)
  if (is.list(counts)) {
    if ("Gene Expression" %in% names(counts)) return(counts[["Gene Expression"]])
    return(counts[[1]])
  }
  counts
}

.find_10x_dir <- function(root_dir) {
  all_dirs <- unique(c(root_dir, list.dirs(root_dir, recursive = TRUE, full.names = TRUE)))
  is_10x_dir <- function(d) {
    files <- list.files(d, full.names = FALSE)
    has_mtx <- any(grepl("^matrix\\.mtx(\\.gz)?$", files, ignore.case = TRUE))
    has_bc  <- any(grepl("^barcodes\\.tsv(\\.gz)?$", files, ignore.case = TRUE))
    has_feat <- any(grepl("^(features|genes)\\.tsv(\\.gz)?$", files, ignore.case = TRUE))
    has_mtx && has_bc && has_feat
  }
  hit <- all_dirs[vapply(all_dirs, is_10x_dir, logical(1))]
  if (length(hit) == 0) return(NULL)
  hit[1]
}

read_10x_dir <- function(dir_path) {
  stopifnot(dir.exists(dir_path))
  counts <- Seurat::Read10X(data.dir = dir_path)
  if (is.list(counts)) {
    if ("Gene Expression" %in% names(counts)) return(counts[["Gene Expression"]])
    return(counts[[1]])
  }
  counts
}

read_10x_input <- function(path) {
  if (dir.exists(path)) return(read_10x_dir(path))
  if (!file.exists(path)) stop("Input not found: ", path)

  ext <- tolower(tools::file_ext(path))

  if (ext == "h5") return(read_10x_h5(path))

  if (ext %in% c("zip", "tar", "gz", "tgz")) {
    tmp <- tempfile("tenx_unpack_")
    dir.create(tmp, recursive = TRUE, showWarnings = FALSE)

    if (ext == "zip") {
      utils::unzip(path, exdir = tmp)
    } else {
      utils::untar(path, exdir = tmp)
    }

    tenx_dir <- .find_10x_dir(tmp)
    if (is.null(tenx_dir)) {
      stop(
        "Archive extracted but no valid 10x MTX directory was found.\n",
        "Expected matrix.mtx(.gz), barcodes.tsv(.gz), and features.tsv(.gz) or genes.tsv(.gz) in the same folder."
      )
    }
    return(read_10x_dir(tenx_dir))
  }

  stop("Unsupported input type: ", path, "\nProvide .h5, a 10x MTX directory, or an archive (.tar.gz/.tgz/.tar/.zip).")
}

derive_batch_names <- function(paths) {
  stems <- tools::file_path_sans_ext(basename(paths))
  parents <- basename(dirname(paths))
  make.unique(paste0(parents, "_", stems))
}
