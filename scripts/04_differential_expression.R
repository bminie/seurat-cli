#!/usr/bin/env Rscript
#' =============================================================================
#' Differential Expression Testing Vignette
#' =============================================================================
#' 
#' This script implements the Seurat differential expression vignette with
#' flexible support for multiple DE testing frameworks:
#'   - Wilcoxon rank sum test (default, fast)
#'   - Student's t-test
#'   - Likelihood ratio test
#'   - MAST (model-based analysis)
#'   - DESeq2 (pseudobulk)
#'   - Negative binomial
#'
#' Supports:
#'   - Finding markers between clusters
#'   - Comparing conditions within clusters
#'   - Custom group comparisons
#'   - Pseudobulk DE analysis
#'
#' Based on: https://satijalab.org/seurat/articles/de_vignette
#'
#' Usage:
#'   Rscript 04_differential_expression.R --input seurat.rds --output ./de_results
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
  make_option(c("-o", "--output"), type = "character", default = "output/de_analysis",
              help = "Output directory [default: %default]", metavar = "DIR"),
  
  # Analysis mode
  make_option("--mode", type = "character", default = "all_markers",
              help = "Analysis mode: all_markers, between_clusters, between_conditions, custom [default: %default]"),
  
  # Cluster-based analysis
  make_option("--cluster_col", type = "character", default = "seurat_clusters",
              help = "Metadata column containing cluster identities [default: %default]"),
  make_option("--ident_1", type = "character", default = NULL,
              help = "Identity class 1 for comparison (for between_clusters mode)"),
  make_option("--ident_2", type = "character", default = NULL,
              help = "Identity class 2 for comparison (NULL = all other cells)"),
  
  # Condition-based analysis
  make_option("--condition_col", type = "character", default = NULL,
              help = "Metadata column for condition grouping"),
  make_option("--condition_1", type = "character", default = NULL,
              help = "Condition 1 (e.g., 'treatment')"),
  make_option("--condition_2", type = "character", default = NULL,
              help = "Condition 2 (e.g., 'control')"),
  make_option("--subset_cluster", type = "character", default = NULL,
              help = "Subset to specific cluster for condition comparison"),
  
  # DE testing parameters
  make_option("--test_use", type = "character", default = "wilcox",
              help = "Test to use: wilcox, bimod, t, poisson, negbinom, LR, MAST, DESeq2 [default: %default]"),
  make_option("--min_pct", type = "double", default = 0.1,
              help = "Minimum fraction of cells expressing gene [default: %default]"),
  make_option("--logfc_threshold", type = "double", default = 0.25,
              help = "Log fold-change threshold [default: %default]"),
  make_option("--only_pos", action = "store_true", default = FALSE,
              help = "Only return positive markers"),
  make_option("--max_cells_per_ident", type = "integer", default = NULL,
              help = "Maximum cells per identity (for downsampling)"),
  make_option("--min_cells_group", type = "integer", default = 3,
              help = "Minimum cells in either group [default: %default]"),
  make_option("--min_cells_feature", type = "integer", default = 3,
              help = "Minimum cells expressing feature [default: %default]"),
  
  # Pseudobulk options
  make_option("--pseudobulk", action = "store_true", default = FALSE,
              help = "Use pseudobulk DE analysis (recommended for conditions)"),
  make_option("--pseudobulk_group", type = "character", default = NULL,
              help = "Metadata column for pseudobulk grouping (e.g., 'sample')"),
  
  # Multiple testing
  make_option("--p_adjust_method", type = "character", default = "BH",
              help = "P-value adjustment method [default: %default]"),
  make_option("--p_threshold", type = "double", default = 0.05,
              help = "Adjusted p-value threshold for significance [default: %default]"),
  
  # Output options
  make_option("--top_n", type = "integer", default = 50,
              help = "Number of top markers to plot [default: %default]"),
  make_option("--features_plot", type = "character", default = NULL,
              help = "Specific features to plot (comma-separated)"),
  
  # General
  make_option(c("-t", "--threads"), type = "integer", default = 1,
              help = "Number of threads [default: %default]"),
  make_option(c("-v", "--verbose"), action = "store_true", default = TRUE,
              help = "Verbose output"),
  make_option(c("-q", "--quiet"), action = "store_false", dest = "verbose",
              help = "Suppress messages"),
  make_option("--seed", type = "integer", default = 42,
              help = "Random seed [default: %default]"),
  make_option("--demo", action = "store_true", default = FALSE,
              help = "Run with demo PBMC3K dataset from SeuratData")
)

parser <- OptionParser(
  usage = "%prog [options]",
  option_list = option_list,
  description = "Differential expression testing using Seurat"
)

args <- parse_args(parser)

# Validate arguments
if (!args$demo && is.null(args$input)) {
  print_help(parser)
  stop("Error: --input is required (or use --demo for demo dataset)")
}

if (!args$demo && !file.exists(args$input)) {
  stop(sprintf("Input file not found: %s", args$input))
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

log_message("Starting differential expression analysis", log_env)
log_message(sprintf("Mode: %s", args$mode), log_env)
log_message(sprintf("Test: %s", args$test_use), log_env)

if (args$threads > 1) {
  library(future)
  plan(multisession, workers = args$threads)
}

# -----------------------------------------------------------------------------
# Load Required Packages
# -----------------------------------------------------------------------------

log_message("Loading required packages...", log_env)
load_seurat_deps(quiet = !args$verbose)

# Install test-specific packages
if (args$test_use == "MAST") {
  install_and_load(bioc_packages = "MAST", quiet = !args$verbose)
}
if (args$test_use == "DESeq2" || args$pseudobulk) {
  install_and_load(bioc_packages = "DESeq2", quiet = !args$verbose)
}

# -----------------------------------------------------------------------------
# Load Data
# -----------------------------------------------------------------------------

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

} else {
  log_message(sprintf("Loading Seurat object: %s", args$input), log_env)
  seurat_obj <- readRDS(args$input)
}

log_message(sprintf("Loaded %d cells, %d features", 
                    ncol(seurat_obj), nrow(seurat_obj)), log_env)

# Set identities
if (args$cluster_col %in% colnames(seurat_obj@meta.data)) {
  Idents(seurat_obj) <- args$cluster_col
  log_message(sprintf("Using identity column: %s", args$cluster_col), log_env)
  log_message(sprintf("Identities: %s", 
                      paste(levels(Idents(seurat_obj)), collapse = ", ")), log_env)
}

# -----------------------------------------------------------------------------
# Differential Expression Analysis
# -----------------------------------------------------------------------------

de_results <- NULL

if (args$mode == "all_markers") {
  # Find markers for all clusters
  log_message("Finding markers for all clusters...", log_env)
  
  de_results <- FindAllMarkers(
    seurat_obj,
    only.pos = args$only_pos,
    min.pct = args$min_pct,
    logfc.threshold = args$logfc_threshold,
    test.use = args$test_use,
    max.cells.per.ident = args$max_cells_per_ident,
    min.cells.group = args$min_cells_group,
    verbose = args$verbose
  )
  
} else if (args$mode == "between_clusters") {
  # Compare two clusters
  if (is.null(args$ident_1)) {
    stop("--ident_1 is required for between_clusters mode")
  }
  
  log_message(sprintf("Finding markers: %s vs %s", 
                      args$ident_1, 
                      ifelse(is.null(args$ident_2), "all others", args$ident_2)), log_env)
  
  de_results <- FindMarkers(
    seurat_obj,
    ident.1 = args$ident_1,
    ident.2 = args$ident_2,
    only.pos = args$only_pos,
    min.pct = args$min_pct,
    logfc.threshold = args$logfc_threshold,
    test.use = args$test_use,
    max.cells.per.ident = args$max_cells_per_ident,
    verbose = args$verbose
  )
  
  de_results$gene <- rownames(de_results)
  de_results$cluster <- args$ident_1
  
} else if (args$mode == "between_conditions") {
  # Compare conditions (within a cluster or globally)
  if (is.null(args$condition_col) || is.null(args$condition_1)) {
    stop("--condition_col and --condition_1 required for between_conditions mode")
  }
  
  # Subset if needed
  obj_subset <- seurat_obj
  if (!is.null(args$subset_cluster)) {
    log_message(sprintf("Subsetting to cluster: %s", args$subset_cluster), log_env)
    obj_subset <- subset(seurat_obj, idents = args$subset_cluster)
  }
  
  # Set identities to condition
  Idents(obj_subset) <- args$condition_col
  
  log_message(sprintf("Comparing conditions: %s vs %s", 
                      args$condition_1,
                      ifelse(is.null(args$condition_2), "all others", args$condition_2)), log_env)
  
  if (args$pseudobulk && !is.null(args$pseudobulk_group)) {
    # Pseudobulk DE analysis
    log_message("Running pseudobulk analysis...", log_env)
    
    # Create pseudobulk counts
    bulk <- AggregateExpression(
      obj_subset,
      group.by = c(args$condition_col, args$pseudobulk_group),
      return.seurat = TRUE,
      verbose = args$verbose
    )
    
    Idents(bulk) <- args$condition_col
    
    de_results <- FindMarkers(
      bulk,
      ident.1 = args$condition_1,
      ident.2 = args$condition_2,
      test.use = "DESeq2",
      verbose = args$verbose
    )
    
    de_results$gene <- rownames(de_results)
    
  } else {
    # Standard DE
    de_results <- FindMarkers(
      obj_subset,
      ident.1 = args$condition_1,
      ident.2 = args$condition_2,
      only.pos = args$only_pos,
      min.pct = args$min_pct,
      logfc.threshold = args$logfc_threshold,
      test.use = args$test_use,
      verbose = args$verbose
    )
    
    de_results$gene <- rownames(de_results)
  }
  
  de_results$comparison <- sprintf("%s_vs_%s", 
                                    args$condition_1, 
                                    ifelse(is.null(args$condition_2), "other", args$condition_2))
  
} else if (args$mode == "custom") {
  # Custom comparison using cells.1 and cells.2 (advanced)
  log_message("Custom mode - please modify script for specific needs", log_env)
  stop("Custom mode requires script modification")
}

# -----------------------------------------------------------------------------
# Process Results
# -----------------------------------------------------------------------------

log_message("Processing DE results...", log_env)

# Check if any DE results were found
if (is.null(de_results) || nrow(de_results) == 0) {
  log_message("WARNING: No differential expression results found!", log_env)
  log_message("This may be due to:", log_env)
  log_message("  - Too stringent filtering parameters (min.pct, logfc.threshold)", log_env)
  log_message("  - Insufficient cells per cluster", log_env)
  log_message("  - Data preprocessing issues", log_env)

  # Create empty results files
  empty_df <- data.frame(
    gene = character(),
    p_val = numeric(),
    avg_log2FC = numeric(),
    pct.1 = numeric(),
    pct.2 = numeric(),
    p_val_adj = numeric()
  )
  write.csv(empty_df, file.path(args$output, "de_results_all.csv"), row.names = FALSE)
  write.csv(empty_df, file.path(args$output, "de_results_significant.csv"), row.names = FALSE)

  # Save parameters and session info
  params_df <- data.frame(
    parameter = names(args),
    value = sapply(args, function(x) paste(x, collapse = ","))
  )
  write.csv(params_df, file.path(args$output, "parameters.csv"), row.names = FALSE)
  writeLines(capture.output(sessionInfo()), file.path(args$output, "session_info.txt"))

  log_message("Empty results files saved", log_env)
  log_message("Analysis complete (no DE genes found)", log_env)
  quit(status = 0)
}

# Add adjusted p-values if not present
if (!"p_val_adj" %in% colnames(de_results) && "p_val" %in% colnames(de_results)) {
  de_results$p_val_adj <- p.adjust(de_results$p_val, method = args$p_adjust_method)
}

# Filter significant results
sig_results <- de_results %>%
  filter(p_val_adj < args$p_threshold)

log_message(sprintf("Total DE genes: %d", nrow(de_results)), log_env)
log_message(sprintf("Significant (adj.p < %.2f): %d",
                    args$p_threshold, nrow(sig_results)), log_env)

# -----------------------------------------------------------------------------
# Save Results
# -----------------------------------------------------------------------------

log_message("Saving results...", log_env)

# All results
write.csv(de_results, file.path(args$output, "de_results_all.csv"), row.names = FALSE)

# Significant results
write.csv(sig_results, file.path(args$output, "de_results_significant.csv"), row.names = FALSE)

# Top markers per cluster (if applicable)
if ("cluster" %in% colnames(de_results)) {
  top_markers <- de_results %>%
    filter(p_val_adj < args$p_threshold) %>%
    group_by(cluster) %>%
    slice_max(n = 10, order_by = avg_log2FC)
  
  write.csv(top_markers, file.path(args$output, "top_markers_per_cluster.csv"), row.names = FALSE)
}

# -----------------------------------------------------------------------------
# Visualization
# -----------------------------------------------------------------------------

log_message("Generating visualizations...", log_env)

# Volcano plot function
make_volcano <- function(df, title = "Volcano Plot") {
  df$significant <- df$p_val_adj < args$p_threshold
  df$label <- ifelse(df$significant & abs(df$avg_log2FC) > 1, df$gene, "")
  
  ggplot(df, aes(x = avg_log2FC, y = -log10(p_val_adj))) +
    geom_point(aes(color = significant), alpha = 0.6, size = 1) +
    geom_text(aes(label = label), size = 2.5, hjust = 0, vjust = 0, 
              check_overlap = TRUE, nudge_x = 0.1) +
    scale_color_manual(values = c("gray60", "red3")) +
    geom_hline(yintercept = -log10(args$p_threshold), linetype = "dashed", color = "gray40") +
    geom_vline(xintercept = c(-args$logfc_threshold, args$logfc_threshold), 
               linetype = "dashed", color = "gray40") +
    theme_minimal() +
    labs(title = title, x = "Log2 Fold Change", y = "-Log10 Adjusted P-value") +
    theme(legend.position = "none")
}

# Generate volcano plots
if (args$mode == "all_markers") {
  for (clust in unique(de_results$cluster)) {
    clust_results <- de_results %>% filter(cluster == clust)
    p_volcano <- make_volcano(clust_results, title = sprintf("Cluster %s Markers", clust))
    save_plot(p_volcano, sprintf("volcano_cluster_%s.png", clust), 
              output_dir = args$output, width = 8, height = 6)
  }
} else {
  p_volcano <- make_volcano(de_results, title = "Differential Expression")
  save_plot(p_volcano, "volcano_plot.png", output_dir = args$output, width = 8, height = 6)
}

# Heatmap of top markers
if (args$mode == "all_markers" && nrow(sig_results) > 0) {
  top_genes <- sig_results %>%
    group_by(cluster) %>%
    slice_max(n = 5, order_by = avg_log2FC) %>%
    pull(gene) %>%
    unique()
  
  if (length(top_genes) > 0) {
    p_heatmap <- DoHeatmap(seurat_obj, features = top_genes) +
      theme(text = element_text(size = 8))
    save_plot(p_heatmap, "heatmap_top_markers.png", output_dir = args$output,
              width = 12, height = max(8, length(top_genes) * 0.25))
  }
}

# Feature plots for top markers
if (nrow(sig_results) > 0) {
  top_features <- head(sig_results$gene[order(sig_results$p_val_adj)], 
                       min(args$top_n, 9))
  
  # Check for UMAP
  if ("umap" %in% Reductions(seurat_obj)) {
    p_features <- FeaturePlot(seurat_obj, features = top_features, 
                               ncol = 3, reduction = "umap")
    save_plot(p_features, "feature_plots_top_markers.png", output_dir = args$output,
              width = 12, height = ceiling(length(top_features)/3) * 4)
  }
  
  # Violin plots
  p_violin <- VlnPlot(seurat_obj, features = head(top_features, 6), 
                       pt.size = 0, ncol = 3)
  save_plot(p_violin, "violin_plots_top_markers.png", output_dir = args$output,
            width = 12, height = 8)
}

# Specific features if provided
if (!is.null(args$features_plot)) {
  features <- trimws(strsplit(args$features_plot, ",")[[1]])
  features <- features[features %in% rownames(seurat_obj)]
  
  if (length(features) > 0) {
    p_custom <- FeaturePlot(seurat_obj, features = features, ncol = 3)
    save_plot(p_custom, "feature_plots_custom.png", output_dir = args$output,
              width = 4 * min(length(features), 3), 
              height = 4 * ceiling(length(features)/3))
  }
}

# -----------------------------------------------------------------------------
# Summary Statistics
# -----------------------------------------------------------------------------

# Summary by cluster (for all_markers mode)
if (args$mode == "all_markers") {
  summary_df <- de_results %>%
    group_by(cluster) %>%
    summarise(
      total_de_genes = n(),
      significant = sum(p_val_adj < args$p_threshold),
      upregulated = sum(avg_log2FC > 0 & p_val_adj < args$p_threshold),
      downregulated = sum(avg_log2FC < 0 & p_val_adj < args$p_threshold),
      mean_logfc = mean(avg_log2FC),
      .groups = "drop"
    )
  
  write.csv(summary_df, file.path(args$output, "de_summary_by_cluster.csv"), row.names = FALSE)
  
  log_message("\nDE Summary by Cluster:", log_env)
  print(summary_df)
}

# Overall summary
overall_summary <- list(
  mode = args$mode,
  test_used = args$test_use,
  total_genes_tested = nrow(de_results),
  significant_genes = nrow(sig_results),
  upregulated = sum(sig_results$avg_log2FC > 0, na.rm = TRUE),
  downregulated = sum(sig_results$avg_log2FC < 0, na.rm = TRUE),
  logfc_threshold = args$logfc_threshold,
  pvalue_threshold = args$p_threshold,
  min_pct = args$min_pct
)

writeLines(
  capture.output(print(as.data.frame(overall_summary))),
  file.path(args$output, "analysis_summary.txt")
)

# Session info
writeLines(capture.output(sessionInfo()), file.path(args$output, "session_info.txt"))

# Parameters
params <- as.data.frame(t(as.data.frame(args)))
write.csv(params, file.path(args$output, "parameters.csv"))

log_message("Differential expression analysis complete!", log_env)
log_message(sprintf("Results saved to: %s", args$output), log_env)
