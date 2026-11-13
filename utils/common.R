#' =============================================================================
#' Seurat CLI - Common Utility Functions
#' =============================================================================
#'
#' Shared utility functions used across all Seurat CLI scripts.
#'
#' @author Generated for Seurat CLI Project
#' =============================================================================

# -----------------------------------------------------------------------------
# Package Management
# -----------------------------------------------------------------------------

#' Install and load required packages
#'
#' @param packages Character vector of package names
#' @param bioc_packages Character vector of Bioconductor package names
#' @param quiet Logical, suppress messages
install_and_load <- function(packages = NULL, bioc_packages = NULL, quiet = FALSE) {

  # Install CRAN packages if not present
  if (!is.null(packages)) {
    for (pkg in packages) {
      if (!requireNamespace(pkg, quietly = TRUE)) {
        if (!quiet) message(sprintf("Installing package: %s", pkg))
        install.packages(pkg, repos = "https://cloud.r-project.org", quiet = quiet)
      }
      suppressPackageStartupMessages(library(pkg, character.only = TRUE))
    }
  }

  # Install Bioconductor packages if not present
  if (!is.null(bioc_packages)) {
    if (!requireNamespace("BiocManager", quietly = TRUE)) {
      install.packages("BiocManager", repos = "https://cloud.r-project.org", quiet = quiet)
    }
    for (pkg in bioc_packages) {
      if (!requireNamespace(pkg, quietly = TRUE)) {
        if (!quiet) message(sprintf("Installing Bioconductor package: %s", pkg))
        BiocManager::install(pkg, ask = FALSE, update = FALSE, quiet = quiet)
      }
      suppressPackageStartupMessages(library(pkg, character.only = TRUE))
    }
  }
}

#' Oldest Seurat release the scripts support
#'
#' install_seurat.R and the Docker image install at least this version.
seurat_min_version <- "6.0.0"

#' Stop with install instructions if the installed Seurat is too old
check_seurat_version <- function() {
  installed <- utils::packageVersion("Seurat")
  if (installed < seurat_min_version) {
    stop(sprintf(paste0(
      "Seurat >= %s is required, but %s is installed.\n",
      "Upgrade it with: Rscript install_seurat.R (from the seurat-cli directory).\n",
      "conda-forge may not have Seurat %s yet; the script installs it from CRAN."),
      seurat_min_version, installed, seurat_min_version), call. = FALSE)
  }
  invisible(installed)
}

#' Use n threads for Seurat's multithreaded steps
#'
#' Seurat (>= 5.6) runs NormalizeData, FindVariableFeatures, ScaleData,
#' SCTransform, RunPCA, FindNeighbors, FindClusters, RunUMAP and CCA/RPCA
#' integration in multithreaded C++ code, set with setThreads(). The threads
#' share one copy of the data, unlike future's multisession workers, which each
#' hold a copy. Call after load_seurat_deps().
#'
#' @param n Requested number of threads
#' @param log_env Logging environment from init_logging()
#' @return Number of threads used (invisibly)
setup_threads <- function(n, log_env = NULL) {
  n <- max(1L, as.integer(n))
  available <- parallel::detectCores()
  if (!is.na(available) && n > available) {
    log_message(sprintf("Requested %d threads but only %d cores are available; using %d",
                        n, available, available), log_env)
    n <- as.integer(available)
  }
  Seurat::setThreads(n, verbose = FALSE)
  if (n > 1) log_message(sprintf("Using %d threads", n), log_env)
  invisible(n)
}

#' Load core Seurat dependencies
#'
#' SeuratData is only needed for demo mode, so it is loaded by
#' load_demo_dataset() rather than here.
#'
#' @param quiet Logical, suppress messages
load_seurat_deps <- function(quiet = FALSE) {
  if (requireNamespace("Seurat", quietly = TRUE)) check_seurat_version()
  core_packages <- c("Seurat", "dplyr", "ggplot2", "patchwork")
  install_and_load(packages = core_packages, quiet = quiet)
}

#' Load a SeuratData demo dataset, updated to the current Seurat object format
#'
#' Installs the dataset on first use. SeuratData itself is on GitHub, not
#' CRAN, so it must already be installed (setup.R does this).
#'
#' @param dataset Dataset name, e.g. "pbmc3k" or "ifnb"
#' @return Seurat object
load_demo_dataset <- function(dataset) {
  if (!requireNamespace("SeuratData", quietly = TRUE)) {
    stop("Demo mode needs the SeuratData package (GitHub only). ",
         "Install it with: Rscript setup.R")
  }
  suppressPackageStartupMessages(library(SeuratData))

  if (!requireNamespace(paste0(dataset, ".SeuratData"), quietly = TRUE)) {
    message(sprintf("Installing demo dataset: %s", dataset))
    SeuratData::InstallData(dataset)
  }

  data_env <- new.env()
  utils::data(list = dataset, envir = data_env)
  UpdateSeuratObject(data_env[[dataset]])
}

#' Directory for downloaded demo files
#'
#' Defaults to the user cache directory for seurat-cli; set SEURAT_CLI_CACHE
#' to use a different one.
cache_dir <- function() {
  dir <- Sys.getenv("SEURAT_CLI_CACHE", tools::R_user_dir("seurat-cli", which = "cache"))
  dir.create(dir, recursive = TRUE, showWarnings = FALSE)
  dir
}

#' Download a file once, then reuse the cached copy
#'
#' @param url File URL
#' @param filename Name for the cached file
#' @return Path to the cached file
cached_download <- function(url, filename = basename(url)) {
  path <- file.path(cache_dir(), filename)
  if (!file.exists(path)) {
    message(sprintf("Downloading %s", url))
    partial <- paste0(path, ".part")
    status <- tryCatch(utils::download.file(url, partial, mode = "wb", quiet = TRUE),
                       error = function(e) 1L)
    if (!identical(as.integer(status), 0L) || !file.exists(partial)) {
      unlink(partial)
      stop(sprintf("Could not download %s", url), call. = FALSE)
    }
    file.rename(partial, path)
  }
  path
}

#' Add donor IDs to the ifnb demo dataset
#'
#' ifnb pools cells from 8 donors per condition. Donor identities come from
#' the demuxlet calls published with Kang et al. 2018, as in Seurat's
#' differential expression tutorial: only confident singlet calls are kept,
#' and cells without a donor ID are dropped. The files (~5 MB) are downloaded
#' once and cached.
#'
#' @param ifnb ifnb Seurat object from load_demo_dataset("ifnb")
#' @return ifnb with a donor_id column, subset to cells with a donor ID
add_ifnb_donor_ids <- function(ifnb) {
  base_url <- "https://raw.githubusercontent.com/yelabucsf/demuxlet_paper_code/master/fig3/"
  read_calls <- function(file) {
    utils::read.table(cached_download(paste0(base_url, file)), header = TRUE,
                      stringsAsFactors = FALSE)
  }
  info <- rbind(read_calls("ye1.ctrl.8.10.sm.best"), read_calls("ye2.stim.8.10.sm.best"))

  # ifnb cell names use "." where the demuxlet barcodes use "-"
  info$BARCODE <- gsub("-", ".", info$BARCODE, fixed = TRUE)
  info <- info[grep("SNG", info$BEST), ]
  # Barcodes seen in both conditions can't be assigned unambiguously
  info <- info[!duplicated(info$BARCODE) & !duplicated(info$BARCODE, fromLast = TRUE), ]

  donor <- stats::setNames(info$BEST, info$BARCODE)[colnames(ifnb)]
  ifnb$donor_id <- unname(donor)
  ifnb[, !is.na(ifnb$donor_id)]
}

# -----------------------------------------------------------------------------
# Data Loading Functions
# -----------------------------------------------------------------------------

#' Load data based on input format
#'
#' @param input_path Path to input data
#' @param input_format Format of input data (10x_h5, 10x_mtx, rds, csv, seurat_data)
#' @param sample_name Sample identifier
#' @param min_cells Minimum cells per feature
#' @param min_features Minimum features per cell
#' @return Seurat object
load_seurat_data <- function(input_path,
                              input_format = "auto",
                              sample_name = "Sample",
                              min_cells = 3,
                              min_features = 200) {

  # Auto-detect format
  if (input_format == "auto") {
    input_format <- detect_input_format(input_path)
  }

  message(sprintf("Loading data from: %s (format: %s)", input_path, input_format))

  counts <- switch(input_format,
    "10x_h5" = {
      Read10X_h5(input_path)
    },
    "10x_mtx" = {
      Read10X(data.dir = input_path)
    },
    "rds" = {
      obj <- readRDS(input_path)
      if (inherits(obj, "Seurat")) {
        return(obj)  # Already a Seurat object
      }
      obj  # Assume it's a counts matrix
    },
    "csv" = {
      as.sparse(read.csv(input_path, row.names = 1))
    },
    "tsv" = {
      as.sparse(read.delim(input_path, row.names = 1))
    },
    "seurat_data" = {
      # Load from SeuratData package
      if (!requireNamespace("SeuratData", quietly = TRUE)) {
        stop("SeuratData package required for seurat_data format")
      }
      library(SeuratData)
      # input_path should be the dataset name
      InstallData(input_path)
      data(list = input_path, package = "SeuratData")
      return(get(input_path))
    },
    stop(sprintf("Unknown input format: %s", input_format))
  )

  # Handle multi-modal data (e.g., CITE-seq with Gene Expression and ADT)
  if (is.list(counts) && !inherits(counts, "dgCMatrix")) {
    message("Multi-modal data detected. Creating object with primary assay.")
    seurat_obj <- CreateSeuratObject(
      counts = counts[[1]],
      project = sample_name,
      min.cells = min_cells,
      min.features = min_features
    )
    # Add additional assays
    for (assay_name in names(counts)[-1]) {
      seurat_obj[[assay_name]] <- CreateAssayObject(counts = counts[[assay_name]])
    }
  } else {
    seurat_obj <- CreateSeuratObject(
      counts = counts,
      project = sample_name,
      min.cells = min_cells,
      min.features = min_features
    )
  }

  return(seurat_obj)
}

#' Detect input format from file path
#'
#' @param input_path Path to input data
#' @return String indicating format
detect_input_format <- function(input_path) {
  if (dir.exists(input_path)) {
    # Check if it's a 10x directory
    if (any(file.exists(file.path(input_path, c("matrix.mtx", "matrix.mtx.gz"))))) {
      return("10x_mtx")
    }
    stop("Directory provided but no recognized format found")
  }

  ext <- tolower(tools::file_ext(input_path))
  switch(ext,
    "h5" = "10x_h5",
    "rds" = "rds",
    "csv" = "csv",
    "tsv" = "tsv",
    "txt" = "tsv",
    stop(sprintf("Cannot auto-detect format for extension: %s", ext))
  )
}

# -----------------------------------------------------------------------------
# QC and Filtering Functions
# -----------------------------------------------------------------------------

#' Calculate common QC metrics
#'
#' @param seurat_obj Seurat object
#' @param species Species for mitochondrial gene pattern ("human" or "mouse")
#' @return Seurat object with QC metrics
calculate_qc_metrics <- function(seurat_obj, species = "human") {

  # Mitochondrial gene patterns
  mt_pattern <- switch(tolower(species),
    "human" = "^MT-",
    "mouse" = "^mt-",
    "^MT-|^mt-"  # Default: try both
  )

  # Calculate percent mitochondrial
  seurat_obj[["percent.mt"]] <- PercentageFeatureSet(seurat_obj, pattern = mt_pattern)

  # Calculate percent ribosomal (optional)
  rp_pattern <- switch(tolower(species),
    "human" = "^RP[SL]",
    "mouse" = "^Rp[sl]",
    "^RP[SL]|^Rp[sl]"
  )
  seurat_obj[["percent.ribo"]] <- PercentageFeatureSet(seurat_obj, pattern = rp_pattern)

  return(seurat_obj)
}

#' Filter cells based on QC metrics
#'
#' @param seurat_obj Seurat object
#' @param min_features Minimum features per cell
#' @param max_features Maximum features per cell
#' @param max_mt Maximum percent mitochondrial
#' @return Filtered Seurat object
filter_cells <- function(seurat_obj,
                          min_features = 200,
                          max_features = Inf,
                          max_mt = 5) {

  n_before <- ncol(seurat_obj)

  seurat_obj <- subset(
    seurat_obj,
    subset = nFeature_RNA > min_features &
             nFeature_RNA < max_features &
             percent.mt < max_mt
  )

  n_after <- ncol(seurat_obj)
  message(sprintf("Filtered cells: %d -> %d (removed %d)",
                  n_before, n_after, n_before - n_after))

  return(seurat_obj)
}

# -----------------------------------------------------------------------------
# Output Functions
# -----------------------------------------------------------------------------

#' Save Seurat object
#'
#' @param seurat_obj Seurat object
#' @param output_path Output file path
#' @param compress Logical, compress output
save_seurat <- function(seurat_obj, output_path, compress = TRUE) {
  message(sprintf("Saving Seurat object to: %s", output_path))
  saveRDS(seurat_obj, file = output_path, compress = compress)
}

#' Save plot to file
#'
#' @param plot ggplot object
#' @param filename Output filename
#' @param output_dir Output directory
#' @param width Plot width in inches
#' @param height Plot height in inches
#' @param dpi Resolution
save_plot <- function(plot,
                       filename,
                       output_dir = "output",
                       width = 10,
                       height = 8,
                       dpi = 300) {

  if (!dir.exists(output_dir)) {
    dir.create(output_dir, recursive = TRUE)
  }

  output_path <- file.path(output_dir, filename)
  message(sprintf("Saving plot to: %s", output_path))

  ggsave(
    filename = output_path,
    plot = plot,
    width = width,
    height = height,
    dpi = dpi
  )

  return(output_path)
}

#' Save a ComplexHeatmap object (e.g. from PlotAnnotatedHeatmap()) as a PNG
#'
#' ggsave() can't draw these, so this opens a png device and draws the
#' heatmap with its legends below it.
#'
#' @param ht ComplexHeatmap Heatmap object
#' @param filename Output file name
#' @param output_dir Output directory
#' @param width,height Size in inches
#' @param dpi Resolution
#' @return Path to the saved file
save_complex_heatmap <- function(ht, filename, output_dir = "output",
                                 width = 10, height = 8, dpi = 300) {
  if (!dir.exists(output_dir)) {
    dir.create(output_dir, recursive = TRUE)
  }
  output_path <- file.path(output_dir, filename)
  message(sprintf("Saving plot to: %s", output_path))

  grDevices::png(output_path, width = width, height = height, units = "in", res = dpi)
  on.exit(grDevices::dev.off())
  ComplexHeatmap::draw(ht, heatmap_legend_side = "bottom",
                       annotation_legend_side = "bottom")
  output_path
}

#' Add scaled values for features so DoHeatmap can plot them
#'
#' DoHeatmap only plots features in the scale.data layer, which usually holds
#' just the variable features, so top markers outside that set are dropped.
#' Returns a copy for plotting with the missing features added: SCT residuals
#' for SCT assays, otherwise ScaleData on only the missing features (existing
#' scaled values, including any regression, are kept).
#'
#' @param seurat_obj Seurat object
#' @param features Features to be plotted
#' @return Seurat object (use for plotting only)
scale_for_heatmap <- function(seurat_obj, features) {
  assay <- DefaultAssay(seurat_obj)
  features <- intersect(unique(features), rownames(seurat_obj[[assay]]))
  scaled <- tryCatch(
    rownames(LayerData(seurat_obj, assay = assay, layer = "scale.data")),
    error = function(e) character())
  missing <- setdiff(features, scaled)

  if (length(missing) == 0) {
    return(seurat_obj)
  }

  if (inherits(seurat_obj[[assay]], "SCTAssay")) {
    seurat_obj <- GetResidual(seurat_obj, features = missing, assay = assay, verbose = FALSE)
  } else if (length(scaled) == 0) {
    seurat_obj <- ScaleData(seurat_obj, features = missing, assay = assay, verbose = FALSE)
  } else {
    extra <- ScaleData(seurat_obj, features = missing, assay = assay, verbose = FALSE)
    combined <- rbind(LayerData(seurat_obj, assay = assay, layer = "scale.data"),
                      LayerData(extra, assay = assay, layer = "scale.data"))
    seurat_obj <- SetAssayData(seurat_obj, assay = assay, layer = "scale.data",
                               new.data = combined)
  }

  seurat_obj
}

#' Create timestamped output directory
#'
#' @param base_dir Base directory
#' @param prefix Prefix for directory name
#' @return Path to created directory
create_output_dir <- function(base_dir = "output", prefix = "run") {
  timestamp <- format(Sys.time(), "%Y%m%d_%H%M%S")
  output_dir <- file.path(base_dir, paste0(prefix, "_", timestamp))
  dir.create(output_dir, recursive = TRUE)
  message(sprintf("Created output directory: %s", output_dir))
  return(output_dir)
}

# -----------------------------------------------------------------------------
# Logging Functions
# -----------------------------------------------------------------------------

#' Initialize logging
#'
#' @param log_file Path to log file (NULL for no file logging)
#' @param verbose Logical, print messages to console
init_logging <- function(log_file = NULL, verbose = TRUE) {

  log_env <- new.env()
  log_env$log_file <- log_file
  log_env$verbose <- verbose
  log_env$start_time <- Sys.time()

  if (!is.null(log_file)) {
    log_dir <- dirname(log_file)
    if (!dir.exists(log_dir)) {
      dir.create(log_dir, recursive = TRUE)
    }
    writeLines(
      sprintf("=== Seurat Analysis Log ===\nStarted: %s\n", Sys.time()),
      log_file
    )
  }

  return(log_env)
}

#' Log a message
#'
#' @param msg Message to log
#' @param log_env Logging environment
log_message <- function(msg, log_env = NULL) {
  timestamp <- format(Sys.time(), "[%Y-%m-%d %H:%M:%S]")
  full_msg <- paste(timestamp, msg)

  if (is.null(log_env) || log_env$verbose) {
    message(full_msg)
  }

  if (!is.null(log_env) && !is.null(log_env$log_file)) {
    cat(full_msg, "\n", file = log_env$log_file, append = TRUE)
  }
}

# -----------------------------------------------------------------------------
# Argument Parsing Helpers
# -----------------------------------------------------------------------------

#' Parse common arguments from command line
#'
#' @return optparse OptionParser object with common options
get_common_parser <- function() {
  if (!requireNamespace("optparse", quietly = TRUE)) {
    install.packages("optparse", repos = "https://cloud.r-project.org", quiet = TRUE)
  }
  library(optparse)

  option_list <- list(
    make_option(c("-i", "--input"), type = "character", default = NULL,
                help = "Input file or directory path", metavar = "PATH"),
    make_option(c("-o", "--output"), type = "character", default = "output",
                help = "Output directory [default: %default]", metavar = "DIR"),
    make_option(c("-n", "--name"), type = "character", default = "Sample",
                help = "Sample name [default: %default]", metavar = "NAME"),
    make_option(c("-f", "--format"), type = "character", default = "auto",
                help = "Input format (auto, 10x_h5, 10x_mtx, rds, csv) [default: %default]",
                metavar = "FORMAT"),
    make_option(c("-s", "--species"), type = "character", default = "human",
                help = "Species (human, mouse) [default: %default]", metavar = "SPECIES"),
    make_option(c("-t", "--threads"), type = "integer", default = 1,
                help = "Number of threads [default: %default]", metavar = "INT"),
    make_option(c("-v", "--verbose"), action = "store_true", default = TRUE,
                help = "Verbose output"),
    make_option(c("-q", "--quiet"), action = "store_false", dest = "verbose",
                help = "Suppress messages")
  )

  return(OptionParser(option_list = option_list))
}

#' UMAP initialization methods accepted by --umap_init
#'
#' Passed to RunUMAP(uwot.init = ...), which hands them to uwot's `init`.
umap_init_methods <- c("spectral", "normlaplacian", "random", "lvrandom",
                       "laplacian", "pca", "spca", "agspectral")

#' Stop if --umap_init is not a known uwot initialization method
check_umap_init <- function(init) {
  if (!init %in% umap_init_methods) {
    stop(sprintf("Unknown --umap_init '%s'. Choose one of: %s",
                 init, paste(umap_init_methods, collapse = ", ")), call. = FALSE)
  }
  invisible(init)
}

#' Print analysis summary
#'
#' @param seurat_obj Seurat object
print_summary <- function(seurat_obj) {
  cat("\n=== Seurat Object Summary ===\n")
  cat(sprintf("Cells: %d\n", ncol(seurat_obj)))
  cat(sprintf("Features: %d\n", nrow(seurat_obj)))
  cat(sprintf("Assays: %s\n", paste(Assays(seurat_obj), collapse = ", ")))

  if (length(Reductions(seurat_obj)) > 0) {
    cat(sprintf("Reductions: %s\n", paste(Reductions(seurat_obj), collapse = ", ")))
  }

  if ("seurat_clusters" %in% colnames(seurat_obj@meta.data)) {
    cat(sprintf("Clusters: %d\n", length(unique(seurat_obj$seurat_clusters))))
  }
  cat("=============================\n\n")
}
