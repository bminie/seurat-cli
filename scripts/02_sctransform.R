#!/usr/bin/env Rscript
#' =============================================================================
#' SCTransform Normalization Vignette
#' =============================================================================
#'
#' This script implements the Seurat SCTransform vignette for improved
#' normalization of single-cell RNA-seq data using regularized negative
#' binomial regression.
#'
#' Key features:
#'   - Replaces NormalizeData, ScaleData, and FindVariableFeatures
#'   - Better handling of technical variation
#'   - Optional: regress out confounding variables (e.g., percent.mt)
#'
#' Based on: https://satijalab.org/seurat/articles/sctransform_vignette
#'
#' Usage:
#'   Rscript 02_sctransform.R --input /path/to/data --output ./results
#'
#' =============================================================================

# Source common utilities
get_script_dir <- function() {
  # Method 1: Rscript --file argument
  args <- commandArgs(trailingOnly = FALSE)
  file_arg <- grep("^--file=", args, value = TRUE)
  if (length(file_arg) > 0) {
    return(dirname(normalizePath(sub("^--file=", "", file_arg))))
  }
  # Method 2: source() - search through frames for ofile
  for (i in seq_len(sys.nframe())) {
    ofile <- sys.frame(i)$ofile
    if (!is.null(ofile)) {
      return(dirname(normalizePath(ofile)))
    }
  }
  # Fallback to current directory
  "."
}
script_dir <- get_script_dir()
source(file.path(script_dir, "..", "utils", "common.R"))

# -----------------------------------------------------------------------------
# Argument Parsing
# -----------------------------------------------------------------------------

library(optparse)

option_list <- list(
  # Input/Output
  make_option(c("-i", "--input"), type = "character", default = NULL,
              help = "Input file or directory path (required)", metavar = "PATH"),
  make_option(c("-o", "--output"), type = "character", default = "output/sctransform",
              help = "Output directory [default: %default]", metavar = "DIR"),
  make_option(c("-n", "--name"), type = "character", default = "Sample",
              help = "Project/sample name [default: %default]", metavar = "NAME"),
  make_option(c("-f", "--format"), type = "character", default = "auto",
              help = "Input format [default: %default]"),

  # Species
  make_option(c("-s", "--species"), type = "character", default = "human",
              help = "Species for QC (human, mouse) [default: %default]"),

  # QC Parameters
  make_option("--min_cells", type = "integer", default = 3,
              help = "Minimum cells per feature [default: %default]"),
  make_option("--min_features", type = "integer", default = 200,
              help = "Minimum features per cell [default: %default]"),

  # SCTransform Parameters
  make_option("--n_variable_features", type = "integer", default = 3000,
              help = "Number of variable features [default: %default]"),
  make_option("--vars_to_regress", type = "character", default = NULL,
              help = "Variables to regress out (comma-separated, e.g., 'percent.mt')"),
  make_option("--return_only_var_genes", action = "store_true", default = TRUE,
              help = "Return only variable genes [default: TRUE]"),
  make_option("--vst_flavor", type = "character", default = "v2",
              help = "SCTransform version: v1 or v2 [default: %default]"),

  # Dimensional Reduction
  make_option("--n_pcs", type = "integer", default = 50,
              help = "Number of PCs to compute [default: %default]"),
  make_option("--dims_use", type = "character", default = "1:30",
              help = "PCs to use for clustering/UMAP [default: %default]"),

  # Clustering
  make_option("--resolution", type = "double", default = 0.8,
              help = "Clustering resolution [default: %default]"),

  # Markers
  make_option("--find_markers", action = "store_true", default = TRUE,
              help = "Find cluster markers"),
  make_option("--no_markers", action = "store_false", dest = "find_markers",
              help = "Skip marker finding"),
  make_option("--min_pct", type = "double", default = 0.25,
              help = "Minimum fraction for markers [default: %default]"),
  make_option("--logfc_threshold", type = "double", default = 0.25,
              help = "Log fold-change threshold [default: %default]"),

  # General
  make_option(c("-t", "--threads"), type = "integer", default = 1,
              help = "Number of threads [default: %default]"),
  make_option(c("-v", "--verbose"), action = "store_true", default = TRUE,
              help = "Verbose output"),
  make_option(c("-q", "--quiet"), action = "store_false", dest = "verbose",
              help = "Suppress messages"),
  make_option("--seed", type = "integer", default = 42,
              help = "Random seed [default: %default]"),

  # Demo mode
  make_option("--demo", action = "store_true", default = FALSE,
              help = "Run with demo PBMC dataset from SeuratData")
)

parser <- OptionParser(
  usage = "%prog [options]",
  option_list = option_list,
  description = "SCTransform normalization workflow for single-cell RNA-seq data"
)

args <- parse_args(parser)

# -----------------------------------------------------------------------------
# Validate Arguments
# -----------------------------------------------------------------------------

if (!args$demo && is.null(args$input)) {
  print_help(parser)
  stop("Error: --input is required (or use --demo for demo dataset)")
}

# Parse dims_use
parse_dims <- function(dims_str) {
  if (grepl(":", dims_str)) {
    parts <- as.integer(strsplit(dims_str, ":")[[1]])
    seq(parts[1], parts[2])
  } else {
    as.integer(strsplit(dims_str, ",")[[1]])
  }
}
dims_use <- parse_dims(args$dims_use)

# Parse vars_to_regress
vars_to_regress <- NULL
if (!is.null(args$vars_to_regress)) {
  vars_to_regress <- trimws(strsplit(args$vars_to_regress, ",")[[1]])
}

set.seed(args$seed)

# -----------------------------------------------------------------------------
# Setup
# -----------------------------------------------------------------------------

if (!dir.exists(args$output)) {
  dir.create(args$output, recursive = TRUE)
}

log_file <- file.path(args$output, "analysis.log")
log_env <- init_logging(log_file = log_file, verbose = args$verbose)

log_message("Starting SCTransform analysis workflow", log_env)
log_message(sprintf("Output directory: %s", args$output), log_env)
log_message(sprintf("SCTransform version: %s", args$vst_flavor), log_env)

if (args$threads > 1) {
  library(future)
  plan(multisession, workers = args$threads)
  log_message(sprintf("Using %d threads", args$threads), log_env)
}

# -----------------------------------------------------------------------------
# Load Required Packages
# -----------------------------------------------------------------------------

log_message("Loading required packages...", log_env)
load_seurat_deps(quiet = !args$verbose)

# Ensure glmGamPoi is available for SCTransform v2
if (args$vst_flavor == "v2") {
  if (!requireNamespace("glmGamPoi", quietly = TRUE)) {
    log_message("Installing glmGamPoi for SCTransform v2...", log_env)
    install_and_load(bioc_packages = "glmGamPoi", quiet = !args$verbose)
  }
}

# -----------------------------------------------------------------------------
# Load Data
# -----------------------------------------------------------------------------

if (args$demo) {
  log_message("Loading demo PBMC dataset...", log_env)

  install_and_load(packages = "SeuratData", quiet = !args$verbose)

  if (!"pbmc3k" %in% rownames(installed.packages())) {
    InstallData("pbmc3k")
  }
  data("pbmc3k")
  seurat_obj <- pbmc3k
  rm(pbmc3k)

  # Update to current Seurat version format (v5 compatibility)
  seurat_obj <- UpdateSeuratObject(seurat_obj)

} else {
  log_message(sprintf("Loading data from: %s", args$input), log_env)

  seurat_obj <- load_seurat_data(
    input_path = args$input,
    input_format = args$format,
    sample_name = args$name,
    min_cells = args$min_cells,
    min_features = args$min_features
  )
}

log_message(sprintf("Loaded %d cells and %d features",
                    ncol(seurat_obj), nrow(seurat_obj)), log_env)

# -----------------------------------------------------------------------------
# Calculate QC Metrics
# -----------------------------------------------------------------------------

log_message("Calculating QC metrics...", log_env)
seurat_obj <- calculate_qc_metrics(seurat_obj, species = args$species)

# QC visualization
p_qc <- VlnPlot(seurat_obj,
                features = c("nFeature_RNA", "nCount_RNA", "percent.mt"),
                ncol = 3, pt.size = 0.1)
save_plot(p_qc, "01_qc_violin.png", output_dir = args$output, width = 12, height = 5)

# -----------------------------------------------------------------------------
# SCTransform Normalization
# -----------------------------------------------------------------------------

log_message("Running SCTransform...", log_env)
if (!is.null(vars_to_regress)) {
  log_message(sprintf("  Regressing out: %s", paste(vars_to_regress, collapse = ", ")), log_env)
}

seurat_obj <- SCTransform(
  seurat_obj,
  variable.features.n = args$n_variable_features,
  vars.to.regress = vars_to_regress,
  return.only.var.genes = args$return_only_var_genes,
  vst.flavor = args$vst_flavor,
  verbose = args$verbose
)

log_message(sprintf("SCTransform complete. %d variable features selected.",
                    length(VariableFeatures(seurat_obj))), log_env)

# Variable features plot
top10 <- head(VariableFeatures(seurat_obj), 10)
p_var <- VariableFeaturePlot(seurat_obj)
p_var <- LabelPoints(plot = p_var, points = top10, repel = TRUE)
save_plot(p_var, "02_variable_features.png", output_dir = args$output, width = 10, height = 6)

# -----------------------------------------------------------------------------
# Dimensional Reduction
# -----------------------------------------------------------------------------

log_message(sprintf("Running PCA (%d components)...", args$n_pcs), log_env)

seurat_obj <- RunPCA(seurat_obj, npcs = args$n_pcs, verbose = args$verbose)

# PCA visualization
p_pca <- DimPlot(seurat_obj, reduction = "pca")
save_plot(p_pca, "03_pca.png", output_dir = args$output, width = 8, height = 6)

# Elbow plot
p_elbow <- ElbowPlot(seurat_obj, ndims = min(args$n_pcs, 50))
save_plot(p_elbow, "04_elbow_plot.png", output_dir = args$output, width = 8, height = 5)

# -----------------------------------------------------------------------------
# Clustering and UMAP
# -----------------------------------------------------------------------------

log_message(sprintf("Finding neighbors using dims: %s", args$dims_use), log_env)

seurat_obj <- FindNeighbors(seurat_obj, dims = dims_use, verbose = args$verbose)

log_message(sprintf("Clustering (resolution: %.2f)...", args$resolution), log_env)

seurat_obj <- FindClusters(seurat_obj, resolution = args$resolution, verbose = args$verbose)

n_clusters <- length(unique(seurat_obj$seurat_clusters))
log_message(sprintf("Found %d clusters", n_clusters), log_env)

log_message("Running UMAP...", log_env)

seurat_obj <- RunUMAP(seurat_obj, dims = dims_use, verbose = args$verbose)

# UMAP visualization
p_umap <- DimPlot(seurat_obj, reduction = "umap", label = TRUE, label.size = 5) +
  ggtitle(sprintf("%s - SCTransform UMAP", args$name))
save_plot(p_umap, "05_umap_clusters.png", output_dir = args$output, width = 10, height = 8)

# Feature plots for QC metrics
p_qc_umap <- FeaturePlot(seurat_obj,
                          features = c("nFeature_RNA", "nCount_RNA", "percent.mt"),
                          ncol = 3)
save_plot(p_qc_umap, "06_umap_qc_features.png", output_dir = args$output,
          width = 15, height = 5)

# -----------------------------------------------------------------------------
# Cluster Markers (with PrepSCTFindMarkers)
# -----------------------------------------------------------------------------

if (args$find_markers) {
  log_message("Preparing for marker finding (PrepSCTFindMarkers)...", log_env)

  # This step is important for SCTransform-normalized data
  seurat_obj <- PrepSCTFindMarkers(seurat_obj, verbose = args$verbose)

  log_message("Finding cluster markers...", log_env)

  markers <- FindAllMarkers(
    seurat_obj,
    only.pos = TRUE,
    min.pct = args$min_pct,
    logfc.threshold = args$logfc_threshold,
    verbose = args$verbose
  )

  # Check if markers were found
  if (is.null(markers) || nrow(markers) == 0 || !"cluster" %in% colnames(markers)) {
    log_message("WARNING: No markers found. Skipping marker analysis.", log_env)
    # Create empty marker files
    write.csv(data.frame(), file.path(args$output, "cluster_markers_all.csv"), row.names = FALSE)
    write.csv(data.frame(), file.path(args$output, "cluster_markers_top10.csv"), row.names = FALSE)
  } else {
    # Save markers
    write.csv(markers, file.path(args$output, "cluster_markers_all.csv"), row.names = FALSE)

    # Top markers per cluster
    top_markers <- markers |>
      group_by(cluster) |>
      slice_max(n = 10, order_by = avg_log2FC)

    write.csv(top_markers, file.path(args$output, "cluster_markers_top10.csv"), row.names = FALSE)

    # Marker heatmap
    top5 <- markers |>
      group_by(cluster) |>
      slice_max(n = 5, order_by = avg_log2FC) |>
      pull(gene)

    if (length(top5) > 0) {
      p_heatmap <- DoHeatmap(seurat_obj, features = unique(top5)) +
        theme(text = element_text(size = 8))
      save_plot(p_heatmap, "07_marker_heatmap.png", output_dir = args$output,
                width = 12, height = max(8, n_clusters * 0.8))
    }

    log_message(sprintf("Found %d markers", nrow(markers)), log_env)
  }
}

# -----------------------------------------------------------------------------
# Save Results
# -----------------------------------------------------------------------------

log_message("Saving results...", log_env)

save_seurat(seurat_obj, file.path(args$output, "seurat_sctransform.rds"))

# Cluster assignments
cluster_df <- data.frame(
  cell = colnames(seurat_obj),
  cluster = seurat_obj$seurat_clusters
)
write.csv(cluster_df, file.path(args$output, "cluster_assignments.csv"), row.names = FALSE)

# UMAP coordinates
umap_coords <- as.data.frame(Embeddings(seurat_obj, "umap"))
umap_coords$cell <- rownames(umap_coords)
umap_coords$cluster <- seurat_obj$seurat_clusters
write.csv(umap_coords, file.path(args$output, "umap_coordinates.csv"), row.names = FALSE)

# Save variable features
write.csv(
  data.frame(feature = VariableFeatures(seurat_obj)),
  file.path(args$output, "variable_features.csv"),
  row.names = FALSE
)

# Summary and session info
print_summary(seurat_obj)
writeLines(capture.output(sessionInfo()), file.path(args$output, "session_info.txt"))

# Save parameters
params <- as.data.frame(t(as.data.frame(args)))
write.csv(params, file.path(args$output, "parameters.csv"))

log_message("SCTransform analysis complete!", log_env)
log_message(sprintf("Results saved to: %s", args$output), log_env)
