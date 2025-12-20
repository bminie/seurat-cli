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

#' Load core Seurat dependencies
#'
#' @param include_spatial Logical, include spatial analysis packages
#' @param include_integration Logical, include integration packages
#' @param quiet Logical, suppress messages
load_seurat_deps <- function(include_spatial = FALSE,
                              include_integration = FALSE,
                              quiet = FALSE) {

  core_packages <- c("Seurat", "dplyr", "ggplot2", "patchwork")

  if (include_spatial) {
    core_packages <- c(core_packages, "SeuratData")
  }

  if (include_integration) {
    core_packages <- c(core_packages, "SeuratData")
  }

  install_and_load(packages = core_packages, quiet = quiet)
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
