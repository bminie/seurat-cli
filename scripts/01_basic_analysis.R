#!/usr/bin/env Rscript
#' =============================================================================
#' Basic Analysis - Standard scRNA-seq Analysis Workflow
#' =============================================================================
#'
#' This script implements a standard single-cell RNA-seq analysis workflow
#' including:
#'   - QC and filtering
#'   - Normalization
#'   - Feature selection
#'   - Scaling
#'   - Dimensional reduction (PCA, UMAP)
#'   - Clustering
#'   - Marker identification
#'
#' Based on: https://satijalab.org/seurat/articles/pbmc3k_tutorial
#'
#' Usage:
#'   Rscript 01_basic_analysis.R --input /path/to/data --output ./results
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
  return(".")
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
  make_option(c("-o", "--output"), type = "character", default = "output/pbmc3k",
              help = "Output directory [default: %default]", metavar = "DIR"),
  make_option(c("-n", "--name"), type = "character", default = "PBMC",
              help = "Project/sample name [default: %default]", metavar = "NAME"),
  make_option(c("-f", "--format"), type = "character", default = "auto",
              help = "Input format (auto, 10x_h5, 10x_mtx, rds, csv) [default: %default]"),

  # Species
  make_option(c("-s", "--species"), type = "character", default = "human",
              help = "Species for QC (human, mouse) [default: %default]"),

  # QC Parameters
  make_option("--min_cells", type = "integer", default = 3,
              help = "Minimum cells per feature [default: %default]"),
  make_option("--min_features", type = "integer", default = 200,
              help = "Minimum features per cell [default: %default]"),
  make_option("--max_features", type = "integer", default = 2500,
              help = "Maximum features per cell [default: %default]"),
  make_option("--max_mt", type = "double", default = 5,
              help = "Maximum percent mitochondrial [default: %default]"),

  # Normalization
  make_option("--norm_method", type = "character", default = "LogNormalize",
              help = "Normalization method [default: %default]"),
  make_option("--scale_factor", type = "double", default = 10000,
              help = "Scale factor for normalization [default: %default]"),

  # Variable Features
  make_option("--n_variable_features", type = "integer", default = 2000,
              help = "Number of variable features [default: %default]"),
  make_option("--selection_method", type = "character", default = "vst",
              help = "Variable feature selection method [default: %default]"),

  # Dimensional Reduction
  make_option("--n_pcs", type = "integer", default = 50,
              help = "Number of PCs to compute [default: %default]"),
  make_option("--dims_use", type = "character", default = "1:10",
              help = "PCs to use for clustering/UMAP (e.g., '1:10' or '1,2,3,5') [default: %default]"),

  # Clustering
  make_option("--resolution", type = "double", default = 0.5,
              help = "Clustering resolution [default: %default]"),
  make_option("--algorithm", type = "integer", default = 1,
              help = "Clustering algorithm (1=Louvain, 2=Louvain/multilevel, 3=SLM, 4=Leiden) [default: %default]"),

  # Markers
  make_option("--find_markers", action = "store_true", default = TRUE,
              help = "Find cluster markers [default: TRUE]"),
  make_option("--no_markers", action = "store_false", dest = "find_markers",
              help = "Skip marker finding"),
  make_option("--min_pct", type = "double", default = 0.25,
              help = "Minimum fraction of cells expressing marker [default: %default]"),
  make_option("--logfc_threshold", type = "double", default = 0.25,
              help = "Log fold-change threshold for markers [default: %default]"),
  make_option("--test_use", type = "character", default = "wilcox",
              help = "Statistical test for markers [default: %default]"),

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
              help = "Run with demo PBMC3K dataset from SeuratData")
)

parser <- OptionParser(
  usage = "%prog [options]",
  option_list = option_list,
  description = "Standard scRNA-seq analysis workflow based on Seurat PBMC3K tutorial"
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
    return(seq(parts[1], parts[2]))
  } else {
    return(as.integer(strsplit(dims_str, ",")[[1]]))
  }
}
dims_use <- parse_dims(args$dims_use)

# Set random seed
set.seed(args$seed)

# -----------------------------------------------------------------------------
# Setup
# -----------------------------------------------------------------------------

# Create output directory
if (!dir.exists(args$output)) {
  dir.create(args$output, recursive = TRUE)
}

# Initialize logging
log_file <- file.path(args$output, "analysis.log")
log_env <- init_logging(log_file = log_file, verbose = args$verbose)

log_message("Starting PBMC3K-style analysis workflow", log_env)
log_message(sprintf("Output directory: %s", args$output), log_env)

# Set up parallelization
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

# -----------------------------------------------------------------------------
# Load Data
# -----------------------------------------------------------------------------

if (args$demo) {
  log_message("Loading demo PBMC3K dataset from SeuratData...", log_env)

  install_and_load(packages = "SeuratData", quiet = !args$verbose)

  # Install and load pbmc3k dataset
  if (!"pbmc3k" %in% rownames(installed.packages())) {
    InstallData("pbmc3k")
  }
  data("pbmc3k")
  seurat_obj <- pbmc3k
  rm(pbmc3k)

  # If using raw data
  if ("pbmc3k.final" %in% ls()) {
    seurat_obj <- pbmc3k.final
    rm(pbmc3k.final)
  }

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
# Quality Control
# -----------------------------------------------------------------------------

log_message("Calculating QC metrics...", log_env)

seurat_obj <- calculate_qc_metrics(seurat_obj, species = args$species)

# Generate QC violin plot
p_qc_violin <- VlnPlot(
  seurat_obj,
  features = c("nFeature_RNA", "nCount_RNA", "percent.mt"),
  ncol = 3,
  pt.size = 0.1
)

save_plot(p_qc_violin, "01_qc_violin.png", output_dir = args$output,
          width = 12, height = 5)

# Generate QC scatter plots
p_scatter1 <- FeatureScatter(seurat_obj, feature1 = "nCount_RNA", feature2 = "percent.mt")
p_scatter2 <- FeatureScatter(seurat_obj, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")

save_plot(p_scatter1 + p_scatter2, "02_qc_scatter.png", output_dir = args$output,
          width = 12, height = 5)

# Filter cells
log_message(sprintf("Filtering cells: features %d-%d, MT < %.1f%%",
                    args$min_features, args$max_features, args$max_mt), log_env)

seurat_obj <- filter_cells(
  seurat_obj,
  min_features = args$min_features,
  max_features = args$max_features,
  max_mt = args$max_mt
)

log_message(sprintf("After filtering: %d cells", ncol(seurat_obj)), log_env)

# -----------------------------------------------------------------------------
# Normalization
# -----------------------------------------------------------------------------

log_message(sprintf("Normalizing data (method: %s, scale_factor: %.0f)...",
                    args$norm_method, args$scale_factor), log_env)

seurat_obj <- NormalizeData(
  seurat_obj,
  normalization.method = args$norm_method,
  scale.factor = args$scale_factor,
  verbose = args$verbose
)

# -----------------------------------------------------------------------------
# Feature Selection
# -----------------------------------------------------------------------------

log_message(sprintf("Finding %d variable features (method: %s)...",
                    args$n_variable_features, args$selection_method), log_env)

seurat_obj <- FindVariableFeatures(
  seurat_obj,
  selection.method = args$selection_method,
  nfeatures = args$n_variable_features,
  verbose = args$verbose
)

# Plot variable features
top10 <- head(VariableFeatures(seurat_obj), 10)
p_var1 <- VariableFeaturePlot(seurat_obj)
p_var2 <- LabelPoints(plot = p_var1, points = top10, repel = TRUE, xnudge = 0, ynudge = 0)

save_plot(p_var2, "03_variable_features.png", output_dir = args$output,
          width = 10, height = 6)

# Save variable features list
write.csv(
  data.frame(feature = VariableFeatures(seurat_obj)),
  file.path(args$output, "variable_features.csv"),
  row.names = FALSE
)

# -----------------------------------------------------------------------------
# Scaling
# -----------------------------------------------------------------------------

log_message("Scaling data...", log_env)

seurat_obj <- ScaleData(
  seurat_obj,
  features = rownames(seurat_obj),
  verbose = args$verbose
)

# -----------------------------------------------------------------------------
# PCA
# -----------------------------------------------------------------------------

log_message(sprintf("Running PCA (%d components)...", args$n_pcs), log_env)

seurat_obj <- RunPCA(
  seurat_obj,
  features = VariableFeatures(seurat_obj),
  npcs = args$n_pcs,
  verbose = args$verbose
)

# PCA visualization
p_pca1 <- DimPlot(seurat_obj, reduction = "pca", dims = c(1, 2))
p_pca2 <- DimPlot(seurat_obj, reduction = "pca", dims = c(2, 3))
save_plot(p_pca1 + p_pca2, "04_pca.png", output_dir = args$output, width = 12, height = 5)

# Elbow plot
p_elbow <- ElbowPlot(seurat_obj, ndims = min(args$n_pcs, 30))
save_plot(p_elbow, "05_elbow_plot.png", output_dir = args$output, width = 8, height = 5)

# PCA loadings heatmap
p_heatmap <- DimHeatmap(seurat_obj, dims = 1:min(9, length(dims_use)),
                        cells = 500, balanced = TRUE)
save_plot(p_heatmap, "06_pca_heatmap.png", output_dir = args$output,
          width = 12, height = 12)

# -----------------------------------------------------------------------------
# Clustering
# -----------------------------------------------------------------------------

log_message(sprintf("Finding neighbors using dims: %s", args$dims_use), log_env)

seurat_obj <- FindNeighbors(
  seurat_obj,
  dims = dims_use,
  verbose = args$verbose
)

log_message(sprintf("Clustering (resolution: %.2f, algorithm: %d)...",
                    args$resolution, args$algorithm), log_env)

seurat_obj <- FindClusters(
  seurat_obj,
  resolution = args$resolution,
  algorithm = args$algorithm,
  verbose = args$verbose
)

n_clusters <- length(unique(seurat_obj$seurat_clusters))
log_message(sprintf("Found %d clusters", n_clusters), log_env)

# -----------------------------------------------------------------------------
# UMAP
# -----------------------------------------------------------------------------

log_message("Running UMAP...", log_env)

seurat_obj <- RunUMAP(
  seurat_obj,
  dims = dims_use,
  verbose = args$verbose
)

# UMAP visualization
p_umap <- DimPlot(seurat_obj, reduction = "umap", label = TRUE, label.size = 5) +
  ggtitle(sprintf("%s - UMAP (res=%.2f)", args$name, args$resolution))

save_plot(p_umap, "07_umap_clusters.png", output_dir = args$output,
          width = 10, height = 8)

# -----------------------------------------------------------------------------
# Cluster Markers
# -----------------------------------------------------------------------------

if (args$find_markers) {
  log_message("Finding cluster markers...", log_env)

  markers <- FindAllMarkers(
    seurat_obj,
    only.pos = TRUE,
    min.pct = args$min_pct,
    logfc.threshold = args$logfc_threshold,
    test.use = args$test_use,
    verbose = args$verbose
  )

  # Check if markers were found
  if (is.null(markers) || nrow(markers) == 0 || !"cluster" %in% colnames(markers)) {
    log_message("WARNING: No markers found. Skipping marker analysis.", log_env)
    # Create empty marker files
    write.csv(data.frame(), file.path(args$output, "cluster_markers_all.csv"), row.names = FALSE)
    write.csv(data.frame(), file.path(args$output, "cluster_markers_top10.csv"), row.names = FALSE)
  } else {
    # Save all markers
    write.csv(markers, file.path(args$output, "cluster_markers_all.csv"), row.names = FALSE)

    # Get top markers per cluster
    top_markers <- markers %>%
      group_by(cluster) %>%
      slice_max(n = 10, order_by = avg_log2FC)

    write.csv(top_markers, file.path(args$output, "cluster_markers_top10.csv"), row.names = FALSE)

    # Heatmap of top markers
    top5_per_cluster <- markers %>%
      group_by(cluster) %>%
      slice_max(n = 5, order_by = avg_log2FC) %>%
      pull(gene)

    if (length(top5_per_cluster) > 0) {
      p_marker_heatmap <- DoHeatmap(seurat_obj, features = unique(top5_per_cluster)) +
        theme(text = element_text(size = 8))
      save_plot(p_marker_heatmap, "08_marker_heatmap.png", output_dir = args$output,
                width = 12, height = max(8, n_clusters * 0.8))
    }

    log_message(sprintf("Found %d markers across %d clusters", nrow(markers), n_clusters), log_env)
  }
}

# -----------------------------------------------------------------------------
# Save Results
# -----------------------------------------------------------------------------

log_message("Saving results...", log_env)

# Save Seurat object
save_seurat(seurat_obj, file.path(args$output, "seurat_object.rds"))

# Save cluster assignments
cluster_df <- data.frame(
  cell = colnames(seurat_obj),
  cluster = seurat_obj$seurat_clusters,
  row.names = NULL
)
write.csv(cluster_df, file.path(args$output, "cluster_assignments.csv"), row.names = FALSE)

# Save UMAP coordinates
umap_coords <- as.data.frame(Embeddings(seurat_obj, "umap"))
umap_coords$cell <- rownames(umap_coords)
umap_coords$cluster <- seurat_obj$seurat_clusters
write.csv(umap_coords, file.path(args$output, "umap_coordinates.csv"), row.names = FALSE)

# Print summary
print_summary(seurat_obj)

# Save session info
writeLines(capture.output(sessionInfo()), file.path(args$output, "session_info.txt"))

# Save parameters
params <- as.data.frame(t(as.data.frame(args)))
write.csv(params, file.path(args$output, "parameters.csv"))

log_message("Analysis complete!", log_env)
log_message(sprintf("Results saved to: %s", args$output), log_env)
