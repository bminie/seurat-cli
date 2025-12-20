#!/usr/bin/env Rscript
#' =============================================================================
#' Seurat Visualization Vignette
#' =============================================================================
#'
#' This script provides comprehensive visualization capabilities for Seurat
#' objects, implementing the visualization vignette:
#'   - Dimensional reduction plots (PCA, UMAP, tSNE)
#'   - Feature plots
#'   - Violin/Ridge plots
#'   - Dot plots
#'   - Heatmaps
#'   - Feature scatter plots
#'
#' Based on: https://satijalab.org/seurat/articles/visualization_vignette
#'
#' Usage:
#'   Rscript 06_visualization.R --input seurat.rds --output ./figures
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
              help = "Input Seurat RDS file (required)", metavar = "FILE"),
  make_option(c("-o", "--output"), type = "character", default = "output/visualization",
              help = "Output directory [default: %default]", metavar = "DIR"),

  # Features to plot
  make_option("--features", type = "character", default = NULL,
              help = "Features to plot (comma-separated gene names)"),
  make_option("--marker_file", type = "character", default = NULL,
              help = "CSV file with marker genes (must have 'gene' column)"),
  make_option("--top_markers", type = "integer", default = 10,
              help = "Number of top markers per cluster to plot [default: %default]"),

  # Group/identity options
  make_option("--group_by", type = "character", default = "seurat_clusters",
              help = "Metadata column for grouping [default: %default]"),
  make_option("--split_by", type = "character", default = NULL,
              help = "Metadata column for splitting plots"),

  # Plot types to generate
  make_option("--dim_plots", action = "store_true", default = TRUE,
              help = "Generate dimensional reduction plots"),
  make_option("--feature_plots", action = "store_true", default = TRUE,
              help = "Generate feature plots"),
  make_option("--violin_plots", action = "store_true", default = TRUE,
              help = "Generate violin plots"),
  make_option("--dot_plot", action = "store_true", default = TRUE,
              help = "Generate dot plots"),
  make_option("--heatmap", action = "store_true", default = TRUE,
              help = "Generate heatmaps"),
  make_option("--ridge_plots", action = "store_true", default = FALSE,
              help = "Generate ridge plots"),
  make_option("--scatter_plots", action = "store_true", default = FALSE,
              help = "Generate feature scatter plots"),

  # Reduction to use
  make_option("--reduction", type = "character", default = "umap",
              help = "Reduction for dim/feature plots [default: %default]"),

  # Plot customization
  make_option("--pt_size", type = "double", default = 0.5,
              help = "Point size for scatter plots [default: %default]"),
  make_option("--label_clusters", action = "store_true", default = TRUE,
              help = "Label clusters on dim plots"),
  make_option("--ncol", type = "integer", default = 3,
              help = "Number of columns for multi-panel plots [default: %default]"),
  make_option("--width", type = "double", default = 10,
              help = "Base plot width in inches [default: %default]"),
  make_option("--height", type = "double", default = 8,
              help = "Base plot height in inches [default: %default]"),
  make_option("--dpi", type = "integer", default = 300,
              help = "Plot resolution [default: %default]"),
  make_option("--format", type = "character", default = "png",
              help = "Output format: png, pdf, svg [default: %default]"),

  # General
  make_option(c("-v", "--verbose"), action = "store_true", default = TRUE,
              help = "Verbose output"),
  make_option(c("-q", "--quiet"), action = "store_false", dest = "verbose",
              help = "Suppress messages"),
  make_option("--demo", action = "store_true", default = FALSE,
              help = "Run with demo PBMC3K dataset from SeuratData")
)

parser <- OptionParser(
  usage = "%prog [options]",
  option_list = option_list,
  description = "Generate comprehensive visualizations from a Seurat object"
)

args <- parse_args(parser)

# Validate
if (!args$demo && is.null(args$input)) {
  print_help(parser)
  stop("Error: --input is required (or use --demo for demo dataset)")
}

# -----------------------------------------------------------------------------
# Setup
# -----------------------------------------------------------------------------

if (!dir.exists(args$output)) {
  dir.create(args$output, recursive = TRUE)
}

log_file <- file.path(args$output, "visualization.log")
log_env <- init_logging(log_file = log_file, verbose = args$verbose)

log_message("Starting visualization workflow", log_env)

# -----------------------------------------------------------------------------
# Load Packages and Data
# -----------------------------------------------------------------------------

log_message("Loading packages...", log_env)
load_seurat_deps(quiet = !args$verbose)
library(viridis)

if (args$demo) {
  log_message("Loading demo PBMC3K dataset from SeuratData...", log_env)

  install_and_load(packages = "SeuratData", quiet = !args$verbose)

  if (!"pbmc3k" %in% rownames(installed.packages())) {
    InstallData("pbmc3k")
  }
  data("pbmc3k")
  seurat_obj <- pbmc3k
  rm(pbmc3k)

  # Update to current Seurat version format (v5 compatibility)
  seurat_obj <- UpdateSeuratObject(seurat_obj)

  # Run basic preprocessing for demo
  log_message("Running basic preprocessing for demo data...", log_env)
  seurat_obj <- NormalizeData(seurat_obj, verbose = FALSE)
  seurat_obj <- FindVariableFeatures(seurat_obj, verbose = FALSE)
  seurat_obj <- ScaleData(seurat_obj, verbose = FALSE)
  seurat_obj <- RunPCA(seurat_obj, verbose = FALSE)
  seurat_obj <- FindNeighbors(seurat_obj, dims = 1:10, verbose = FALSE)
  seurat_obj <- FindClusters(seurat_obj, resolution = 0.5, verbose = FALSE)
  seurat_obj <- RunUMAP(seurat_obj, dims = 1:10, verbose = FALSE)

} else {
  log_message(sprintf("Loading Seurat object: %s", args$input), log_env)
  seurat_obj <- readRDS(args$input)
}

log_message(sprintf("Loaded %d cells, %d features",
                    ncol(seurat_obj), nrow(seurat_obj)), log_env)

# Set default identity
if (args$group_by %in% colnames(seurat_obj@meta.data)) {
  Idents(seurat_obj) <- args$group_by
} else {
  log_message(sprintf("WARNING: group_by column '%s' not found in metadata. Using default identities.", args$group_by), log_env)
  log_message(sprintf("  Available columns: %s", paste(head(colnames(seurat_obj@meta.data), 10), collapse = ", ")), log_env)
}

# Validate split_by column if specified
if (!is.null(args$split_by) && !args$split_by %in% colnames(seurat_obj@meta.data)) {
  log_message(sprintf("WARNING: split_by column '%s' not found in metadata. Split plots will be skipped.", args$split_by), log_env)
  args$split_by <- NULL
}

# -----------------------------------------------------------------------------
# Prepare Features
# -----------------------------------------------------------------------------

features_to_plot <- c()

# From command line
if (!is.null(args$features)) {
  features_to_plot <- trimws(strsplit(args$features, ",")[[1]])
}

# From marker file
if (!is.null(args$marker_file) && file.exists(args$marker_file)) {
  markers_df <- read.csv(args$marker_file)
  if ("gene" %in% colnames(markers_df)) {
    marker_genes <- unique(markers_df$gene)
    features_to_plot <- c(features_to_plot, marker_genes)
  }
}

# Filter to genes in dataset
features_to_plot <- unique(features_to_plot)
features_to_plot <- features_to_plot[features_to_plot %in% rownames(seurat_obj)]

log_message(sprintf("Features to plot: %d", length(features_to_plot)), log_env)

# If no features specified, try to get top variable features
if (length(features_to_plot) == 0) {
  log_message("No features specified, using top variable features", log_env)
  features_to_plot <- head(VariableFeatures(seurat_obj), 12)

  # Warn if still no features available
  if (length(features_to_plot) == 0) {
    log_message("WARNING: No variable features found. Feature-based plots will be skipped.", log_env)
    log_message("  Consider running FindVariableFeatures() first or specifying --features", log_env)
  }
}

# Check for available reductions
available_reductions <- Reductions(seurat_obj)

if (length(available_reductions) == 0) {
  log_message("WARNING: No dimensional reductions found in object.", log_env)
  log_message("  Dimensional reduction plots will be skipped.", log_env)
  log_message("  Consider running RunPCA/RunUMAP first.", log_env)
  args$dim_plots <- FALSE
  args$feature_plots <- FALSE
} else {
  log_message(sprintf("Available reductions: %s",
                      paste(available_reductions, collapse = ", ")), log_env)

  if (!args$reduction %in% available_reductions) {
    log_message(sprintf("Warning: %s not found, using first available", args$reduction), log_env)
    args$reduction <- available_reductions[1]
  }
}

# Helper function to save plots
save_figure <- function(plot, name, w = args$width, h = args$height) {
  filename <- paste0(name, ".", args$format)
  filepath <- file.path(args$output, filename)
  ggsave(filepath, plot, width = w, height = h, dpi = args$dpi)
  log_message(sprintf("  Saved: %s", filename), log_env)
  return(filepath)
}

# -----------------------------------------------------------------------------
# Dimensional Reduction Plots
# -----------------------------------------------------------------------------

if (args$dim_plots) {
  log_message("Generating dimensional reduction plots...", log_env)

  # Basic DimPlot
  p_dim <- DimPlot(seurat_obj, reduction = args$reduction,
                   group.by = args$group_by,
                   label = args$label_clusters,
                   pt.size = args$pt_size) +
    ggtitle(sprintf("%s - %s", args$reduction, args$group_by))
  save_figure(p_dim, sprintf("01_dimplot_%s", args$reduction))

  # Split by variable if specified
  if (!is.null(args$split_by) && args$split_by %in% colnames(seurat_obj@meta.data)) {
    n_splits <- length(unique(seurat_obj@meta.data[[args$split_by]]))

    p_split <- DimPlot(seurat_obj, reduction = args$reduction,
                       group.by = args$group_by,
                       split.by = args$split_by,
                       label = args$label_clusters,
                       pt.size = args$pt_size,
                       ncol = min(n_splits, 4))
    save_figure(p_split, sprintf("02_dimplot_split_%s", args$split_by),
                w = min(n_splits, 4) * 5, h = 5 * ceiling(n_splits/4))
  }

  # Additional reductions if available
  for (red in c("pca", "tsne", "umap")) {
    if (red %in% available_reductions && red != args$reduction) {
      p <- DimPlot(seurat_obj, reduction = red,
                   group.by = args$group_by,
                   label = args$label_clusters,
                   pt.size = args$pt_size)
      save_figure(p, sprintf("01_dimplot_%s", red))
    }
  }
}

# -----------------------------------------------------------------------------
# Feature Plots
# -----------------------------------------------------------------------------

if (args$feature_plots && length(features_to_plot) > 0) {
  log_message("Generating feature plots...", log_env)

  # Split features into batches for readability
  batch_size <- 9
  n_batches <- ceiling(length(features_to_plot) / batch_size)

  for (i in seq_len(n_batches)) {
    start_idx <- (i - 1) * batch_size + 1
    end_idx <- min(i * batch_size, length(features_to_plot))
    batch_features <- features_to_plot[start_idx:end_idx]

    p_feat <- FeaturePlot(seurat_obj,
                          features = batch_features,
                          reduction = args$reduction,
                          pt.size = args$pt_size,
                          ncol = args$ncol)

    save_figure(p_feat, sprintf("03_featureplot_batch%d", i),
                w = args$ncol * 4, h = ceiling(length(batch_features)/args$ncol) * 4)
  }

  # Feature plot with split
  if (!is.null(args$split_by) && length(features_to_plot) > 0) {
    top_features <- head(features_to_plot, 4)

    p_feat_split <- FeaturePlot(seurat_obj,
                                 features = top_features,
                                 split.by = args$split_by,
                                 reduction = args$reduction,
                                 pt.size = args$pt_size * 0.5)
    n_splits <- length(unique(seurat_obj@meta.data[[args$split_by]]))
    save_figure(p_feat_split, "03_featureplot_split",
                w = n_splits * 4, h = length(top_features) * 3)
  }
}

# -----------------------------------------------------------------------------
# Violin Plots
# -----------------------------------------------------------------------------

if (args$violin_plots && length(features_to_plot) > 0) {
  log_message("Generating violin plots...", log_env)

  # Violin plots in batches
  batch_size <- 6
  n_batches <- ceiling(length(features_to_plot) / batch_size)

  for (i in seq_len(min(n_batches, 3))) {  # Max 3 batches
    start_idx <- (i - 1) * batch_size + 1
    end_idx <- min(i * batch_size, length(features_to_plot))
    batch_features <- features_to_plot[start_idx:end_idx]

    p_vln <- VlnPlot(seurat_obj,
                      features = batch_features,
                      group.by = args$group_by,
                      pt.size = 0,
                      ncol = min(length(batch_features), 3))

    save_figure(p_vln, sprintf("04_violinplot_batch%d", i),
                w = min(length(batch_features), 3) * 4,
                h = ceiling(length(batch_features)/3) * 4)
  }

  # Split violin
  if (!is.null(args$split_by) && length(features_to_plot) > 0) {
    top_features <- head(features_to_plot, 4)

    p_vln_split <- VlnPlot(seurat_obj,
                            features = top_features,
                            group.by = args$group_by,
                            split.by = args$split_by,
                            pt.size = 0,
                            ncol = 2)
    save_figure(p_vln_split, "04_violinplot_split",
                w = 10, h = ceiling(length(top_features)/2) * 4)
  }
}

# -----------------------------------------------------------------------------
# Dot Plots
# -----------------------------------------------------------------------------

if (args$dot_plot && length(features_to_plot) > 0) {
  log_message("Generating dot plots...", log_env)

  # Use top features for readability
  dot_features <- head(features_to_plot, 30)

  p_dot <- DotPlot(seurat_obj,
                   features = dot_features,
                   group.by = args$group_by) +
    RotatedAxis() +
    theme(axis.text.x = element_text(size = 8))

  save_figure(p_dot, "05_dotplot",
              w = max(8, length(dot_features) * 0.4),
              h = max(6, length(unique(Idents(seurat_obj))) * 0.4))

  # Clustered dot plot (if have cluster column)
  if (!is.null(args$split_by)) {
    p_dot_split <- DotPlot(seurat_obj,
                           features = head(dot_features, 15),
                           group.by = args$group_by,
                           split.by = args$split_by) +
      RotatedAxis()
    save_figure(p_dot_split, "05_dotplot_split",
                w = 12, h = 8)
  }
}

# -----------------------------------------------------------------------------
# Heatmaps
# -----------------------------------------------------------------------------

if (args$heatmap && length(features_to_plot) > 0) {
  log_message("Generating heatmaps...", log_env)

  # Heatmap features
  heatmap_features <- head(features_to_plot, 50)

  p_heat <- DoHeatmap(seurat_obj,
                       features = heatmap_features,
                       group.by = args$group_by,
                       size = 3) +
    theme(text = element_text(size = 8))

  save_figure(p_heat, "06_heatmap",
              w = 12, h = max(8, length(heatmap_features) * 0.2))

  # Scaled heatmap for top markers per cluster (if markers available)
  if (!is.null(args$marker_file) && file.exists(args$marker_file)) {
    markers_df <- read.csv(args$marker_file)
    if ("cluster" %in% colnames(markers_df) && "gene" %in% colnames(markers_df)) {
      top_by_cluster <- markers_df %>%
        group_by(cluster) %>%
        slice_head(n = args$top_markers) %>%
        pull(gene) %>%
        unique()

      top_by_cluster <- top_by_cluster[top_by_cluster %in% rownames(seurat_obj)]

      if (length(top_by_cluster) > 0) {
        p_heat_markers <- DoHeatmap(seurat_obj,
                                     features = top_by_cluster,
                                     group.by = args$group_by,
                                     size = 3) +
          theme(text = element_text(size = 7))
        save_figure(p_heat_markers, "06_heatmap_markers",
                    w = 14, h = max(10, length(top_by_cluster) * 0.15))
      }
    }
  }
}

# -----------------------------------------------------------------------------
# Ridge Plots
# -----------------------------------------------------------------------------

if (args$ridge_plots && length(features_to_plot) > 0) {
  log_message("Generating ridge plots...", log_env)

  ridge_features <- head(features_to_plot, 6)

  p_ridge <- RidgePlot(seurat_obj,
                        features = ridge_features,
                        group.by = args$group_by,
                        ncol = 2)

  save_figure(p_ridge, "07_ridgeplot",
              w = 10, h = ceiling(length(ridge_features)/2) * 4)
}

# -----------------------------------------------------------------------------
# Feature Scatter Plots
# -----------------------------------------------------------------------------

if (args$scatter_plots && length(features_to_plot) >= 2) {
  log_message("Generating feature scatter plots...", log_env)

  # Plot pairs of top features
  n_pairs <- min(4, length(features_to_plot) %/% 2)
  scatter_plots <- list()

  for (i in seq_len(n_pairs)) {
    f1 <- features_to_plot[2*i - 1]
    f2 <- features_to_plot[2*i]

    scatter_plots[[i]] <- FeatureScatter(seurat_obj,
                                          feature1 = f1,
                                          feature2 = f2,
                                          group.by = args$group_by)
  }

  if (length(scatter_plots) > 0) {
    p_scatter <- wrap_plots(scatter_plots, ncol = 2)
    save_figure(p_scatter, "08_feature_scatter",
                w = 10, h = ceiling(length(scatter_plots)/2) * 4)
  }
}

# -----------------------------------------------------------------------------
# QC Plots
# -----------------------------------------------------------------------------

log_message("Generating QC plots...", log_env)

# Available QC metrics
qc_features <- c("nCount_RNA", "nFeature_RNA", "percent.mt", "percent.ribo")
qc_features <- qc_features[qc_features %in% colnames(seurat_obj@meta.data)]

if (length(qc_features) > 0) {
  p_qc_vln <- VlnPlot(seurat_obj,
                       features = qc_features,
                       group.by = args$group_by,
                       pt.size = 0,
                       ncol = length(qc_features))
  save_figure(p_qc_vln, "09_qc_violin", w = length(qc_features) * 3, h = 5)

  p_qc_feat <- FeaturePlot(seurat_obj,
                            features = qc_features,
                            reduction = args$reduction,
                            pt.size = args$pt_size * 0.5,
                            ncol = 2)
  save_figure(p_qc_feat, "09_qc_featureplot",
              w = 8, h = ceiling(length(qc_features)/2) * 4)
}

# -----------------------------------------------------------------------------
# Summary
# -----------------------------------------------------------------------------

# Create summary document
summary_text <- c(
  "# Visualization Summary",
  "",
  sprintf("Input: %s", args$input),
  sprintf("Output: %s", args$output),
  sprintf("Cells: %d", ncol(seurat_obj)),
  sprintf("Features plotted: %d", length(features_to_plot)),
  sprintf("Reduction used: %s", args$reduction),
  sprintf("Group by: %s", args$group_by),
  "",
  "## Generated Figures:",
  list.files(args$output, pattern = paste0("\\.", args$format, "$"))
)

writeLines(summary_text, file.path(args$output, "visualization_summary.md"))

# Save parameters
params_df <- data.frame(
  parameter = names(args),
  value = sapply(args, function(x) paste(x, collapse = ","))
)
write.csv(params_df, file.path(args$output, "parameters.csv"), row.names = FALSE)

# Save session info
writeLines(capture.output(sessionInfo()), file.path(args$output, "session_info.txt"))

log_message("Visualization complete!", log_env)
log_message(sprintf("Figures saved to: %s", args$output), log_env)
