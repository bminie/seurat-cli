#!/usr/bin/env Rscript
#' =============================================================================
#' Cell Cycle Scoring and Regression Vignette
#' =============================================================================
#'
#' This script implements the Seurat cell cycle vignette for:
#'   - Scoring cells for cell cycle phase (G1, S, G2M)
#'   - Optionally regressing out cell cycle effects
#'   - Comparing cell cycle effects before/after regression
#'
#' Based on: https://satijalab.org/seurat/articles/cell_cycle_vignette
#'
#' Usage:
#'   Rscript 05_cell_cycle.R --input seurat.rds --output ./cell_cycle
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
              help = "Input: Seurat RDS file, 10x directory, or raw data", metavar = "PATH"),
  make_option(c("-o", "--output"), type = "character", default = "output/cell_cycle",
              help = "Output directory [default: %default]", metavar = "DIR"),
  make_option(c("-n", "--name"), type = "character", default = "Sample",
              help = "Project/sample name [default: %default]"),
  make_option(c("-f", "--format"), type = "character", default = "auto",
              help = "Input format [default: %default]"),

  # Species
  make_option(c("-s", "--species"), type = "character", default = "human",
              help = "Species (human, mouse) for gene lists [default: %default]"),

  # Cell cycle options
  make_option("--s_genes", type = "character", default = NULL,
              help = "Custom S phase genes (comma-separated) or file path"),
  make_option("--g2m_genes", type = "character", default = NULL,
              help = "Custom G2M phase genes (comma-separated) or file path"),
  make_option("--gene_set_version", type = "character", default = "2019",
              help = "Gene set version: 2019 (updated) or original [default: %default]"),

  # Regression options
  make_option("--regress_cc", action = "store_true", default = FALSE,
              help = "Regress out cell cycle effects"),
  make_option("--regress_cc_difference", action = "store_true", default = FALSE,
              help = "Regress out difference (S - G2M) instead of both scores"),
  make_option("--vars_to_regress", type = "character", default = NULL,
              help = "Additional variables to regress (comma-separated)"),

  # Analysis parameters
  make_option("--n_variable_features", type = "integer", default = 2000,
              help = "Number of variable features [default: %default]"),
  make_option("--n_pcs", type = "integer", default = 50,
              help = "Number of PCs [default: %default]"),
  make_option("--dims_use", type = "character", default = "1:10",
              help = "PCs for UMAP [default: %default]"),
  make_option("--resolution", type = "double", default = 0.5,
              help = "Clustering resolution [default: %default]"),

  # Pre-processing (if starting from raw data)
  make_option("--min_cells", type = "integer", default = 3,
              help = "Minimum cells per feature [default: %default]"),
  make_option("--min_features", type = "integer", default = 200,
              help = "Minimum features per cell [default: %default]"),
  make_option("--max_mt", type = "double", default = 5,
              help = "Maximum percent mitochondrial [default: %default]"),

  # General
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
  description = "Cell cycle scoring and regression for scRNA-seq data"
)

args <- parse_args(parser)

# Validate
if (!args$demo && is.null(args$input)) {
  print_help(parser)
  stop("Error: --input is required (or use --demo for demo dataset)")
}

# Parse dims
parse_dims <- function(dims_str) {
  if (grepl(":", dims_str)) {
    parts <- as.integer(strsplit(dims_str, ":")[[1]])
    seq(parts[1], parts[2])
  } else {
    as.integer(strsplit(dims_str, ",")[[1]])
  }
}
dims_use <- parse_dims(args$dims_use)

set.seed(args$seed)

# -----------------------------------------------------------------------------
# Setup
# -----------------------------------------------------------------------------

if (!dir.exists(args$output)) {
  dir.create(args$output, recursive = TRUE)
}

log_file <- file.path(args$output, "analysis.log")
log_env <- init_logging(log_file = log_file, verbose = args$verbose)

log_message("Starting cell cycle analysis", log_env)

# -----------------------------------------------------------------------------
# Load Required Packages
# -----------------------------------------------------------------------------

log_message("Loading required packages...", log_env)
load_seurat_deps(quiet = !args$verbose)

# -----------------------------------------------------------------------------
# Load or Prepare Cell Cycle Genes
# -----------------------------------------------------------------------------

log_message("Preparing cell cycle gene sets...", log_env)

# Load gene lists
load_gene_list <- function(input, default_genes) {
  if (is.null(input)) {
    return(default_genes)
  }

  if (file.exists(input)) {
    genes <- readLines(input)
  } else {
    genes <- trimws(strsplit(input, ",")[[1]])
  }

  return(genes)
}

# Get default gene sets
if (args$gene_set_version == "2019") {
  s_genes_default <- cc.genes.updated.2019$s.genes
  g2m_genes_default <- cc.genes.updated.2019$g2m.genes
} else {
  s_genes_default <- cc.genes$s.genes
  g2m_genes_default <- cc.genes$g2m.genes
}

# Convert to mouse if needed
if (tolower(args$species) == "mouse") {
  # Convert to mouse gene names (capitalize first letter only)
  convert_to_mouse <- function(genes) {
    sapply(genes, function(g) {
      paste0(toupper(substr(g, 1, 1)), tolower(substr(g, 2, nchar(g))))
    })
  }
  s_genes_default <- convert_to_mouse(s_genes_default)
  g2m_genes_default <- convert_to_mouse(g2m_genes_default)
}

s_genes <- load_gene_list(args$s_genes, s_genes_default)
g2m_genes <- load_gene_list(args$g2m_genes, g2m_genes_default)

log_message(sprintf("S phase genes: %d", length(s_genes)), log_env)
log_message(sprintf("G2M phase genes: %d", length(g2m_genes)), log_env)

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
  needs_preprocessing <- TRUE

} else {
  log_message(sprintf("Loading data from: %s", args$input), log_env)

  # Check if input is already a Seurat object
  if (grepl("\\.rds$", tolower(args$input))) {
    seurat_obj <- readRDS(args$input)
    needs_preprocessing <- !("RNA" %in% Assays(seurat_obj) &&
                             length(GetAssayData(seurat_obj, slot = "data")) > 0)
  } else {
    seurat_obj <- load_seurat_data(
      input_path = args$input,
      input_format = args$format,
      sample_name = args$name,
      min_cells = args$min_cells,
      min_features = args$min_features
    )
    needs_preprocessing <- TRUE
  }
}

log_message(sprintf("Loaded %d cells, %d features",
                    ncol(seurat_obj), nrow(seurat_obj)), log_env)

# -----------------------------------------------------------------------------
# Pre-processing (if needed)
# -----------------------------------------------------------------------------

if (needs_preprocessing) {
  log_message("Running preprocessing...", log_env)

  # QC
  seurat_obj <- calculate_qc_metrics(seurat_obj, species = args$species)

  # Filter
  seurat_obj <- filter_cells(seurat_obj,
                              min_features = args$min_features,
                              max_mt = args$max_mt)

  # Normalize
  seurat_obj <- NormalizeData(seurat_obj, verbose = args$verbose)

  # Variable features
  seurat_obj <- FindVariableFeatures(seurat_obj,
                                      nfeatures = args$n_variable_features,
                                      verbose = args$verbose)
}

# Check which genes are in the dataset
s_genes_present <- s_genes[s_genes %in% rownames(seurat_obj)]
g2m_genes_present <- g2m_genes[g2m_genes %in% rownames(seurat_obj)]

log_message(sprintf("S genes found in data: %d/%d",
                    length(s_genes_present), length(s_genes)), log_env)
log_message(sprintf("G2M genes found in data: %d/%d",
                    length(g2m_genes_present), length(g2m_genes)), log_env)

# Warn if insufficient genes for reliable scoring
if (length(s_genes_present) < 5) {
  log_message("WARNING: Very few S phase genes found in dataset. Cell cycle scoring may be unreliable.", log_env)
}
if (length(g2m_genes_present) < 5) {
  log_message("WARNING: Very few G2M phase genes found in dataset. Cell cycle scoring may be unreliable.", log_env)
}
if (length(s_genes_present) == 0 || length(g2m_genes_present) == 0) {
  log_message("ERROR: No cell cycle genes found for one or both phases. Check species setting or provide custom gene lists.", log_env)
  stop("Insufficient cell cycle genes for scoring")
}

# Save gene lists
writeLines(s_genes_present, file.path(args$output, "s_genes_used.txt"))
writeLines(g2m_genes_present, file.path(args$output, "g2m_genes_used.txt"))

# -----------------------------------------------------------------------------
# Cell Cycle Scoring
# -----------------------------------------------------------------------------

log_message("Scoring cell cycle phases...", log_env)

seurat_obj <- CellCycleScoring(
  seurat_obj,
  s.features = s_genes,
  g2m.features = g2m_genes,
  set.ident = FALSE
)

# Add cell cycle difference score
seurat_obj$CC.Difference <- seurat_obj$S.Score - seurat_obj$G2M.Score

# Summarize - handle cases where a phase might have zero cells
cc_summary <- table(seurat_obj$Phase)
g1_count <- ifelse("G1" %in% names(cc_summary), cc_summary["G1"], 0)
s_count <- ifelse("S" %in% names(cc_summary), cc_summary["S"], 0)
g2m_count <- ifelse("G2M" %in% names(cc_summary), cc_summary["G2M"], 0)
log_message(sprintf("Cell cycle phases: G1=%d, S=%d, G2M=%d", g1_count, s_count, g2m_count), log_env)

# Save cell cycle scores
cc_scores <- data.frame(
  cell = colnames(seurat_obj),
  S.Score = seurat_obj$S.Score,
  G2M.Score = seurat_obj$G2M.Score,
  Phase = seurat_obj$Phase,
  CC.Difference = seurat_obj$CC.Difference
)
write.csv(cc_scores, file.path(args$output, "cell_cycle_scores.csv"), row.names = FALSE)

# -----------------------------------------------------------------------------
# Visualize Before Regression
# -----------------------------------------------------------------------------

log_message("Generating pre-regression visualizations...", log_env)

# Scale and run PCA
seurat_obj <- ScaleData(seurat_obj, verbose = args$verbose)
seurat_obj <- RunPCA(seurat_obj, features = VariableFeatures(seurat_obj),
                      npcs = args$n_pcs, verbose = args$verbose)

# PCA colored by phase
p_pca_phase <- DimPlot(seurat_obj, reduction = "pca", group.by = "Phase") +
  ggtitle("PCA - Before CC Regression")
save_plot(p_pca_phase, "01_pca_by_phase_before.png", output_dir = args$output)

# PCA colored by scores
p_pca_s <- FeaturePlot(seurat_obj, features = "S.Score", reduction = "pca") +
  ggtitle("S Score - Before Regression")
p_pca_g2m <- FeaturePlot(seurat_obj, features = "G2M.Score", reduction = "pca") +
  ggtitle("G2M Score - Before Regression")
save_plot(p_pca_s + p_pca_g2m, "02_pca_scores_before.png",
          output_dir = args$output, width = 12, height = 5)

# UMAP
seurat_obj <- RunUMAP(seurat_obj, dims = dims_use, reduction.name = "umap_before",
                       verbose = args$verbose)

p_umap_phase <- DimPlot(seurat_obj, reduction = "umap_before", group.by = "Phase") +
  ggtitle("UMAP - Before CC Regression")
save_plot(p_umap_phase, "03_umap_by_phase_before.png", output_dir = args$output)

# Score distribution
p_scores <- ggplot(cc_scores, aes(x = S.Score, y = G2M.Score, color = Phase)) +
  geom_point(alpha = 0.5, size = 0.5) +
  theme_minimal() +
  labs(title = "Cell Cycle Scores", x = "S Score", y = "G2M Score")
save_plot(p_scores, "04_cc_score_distribution.png", output_dir = args$output)

# -----------------------------------------------------------------------------
# Cell Cycle Regression (Optional)
# -----------------------------------------------------------------------------

if (args$regress_cc || args$regress_cc_difference) {
  log_message("Regressing out cell cycle effects...", log_env)

  # Determine what to regress
  vars_to_regress <- c()

  if (args$regress_cc_difference) {
    vars_to_regress <- c("CC.Difference")
    log_message("  Regressing CC.Difference (S.Score - G2M.Score)", log_env)
  } else if (args$regress_cc) {
    vars_to_regress <- c("S.Score", "G2M.Score")
    log_message("  Regressing S.Score and G2M.Score", log_env)
  }

  # Add additional variables if specified
  if (!is.null(args$vars_to_regress)) {
    extra_vars <- trimws(strsplit(args$vars_to_regress, ",")[[1]])
    vars_to_regress <- c(vars_to_regress, extra_vars)
    log_message(sprintf("  Also regressing: %s", paste(extra_vars, collapse = ", ")), log_env)
  }

  # Re-scale with regression
  seurat_obj <- ScaleData(
    seurat_obj,
    vars.to.regress = vars_to_regress,
    features = rownames(seurat_obj),
    verbose = args$verbose
  )

  # Re-run PCA
  seurat_obj <- RunPCA(seurat_obj, features = VariableFeatures(seurat_obj),
                        npcs = args$n_pcs, reduction.name = "pca_regressed",
                        verbose = args$verbose)

  # Visualize after regression
  log_message("Generating post-regression visualizations...", log_env)

  p_pca_phase_after <- DimPlot(seurat_obj, reduction = "pca_regressed", group.by = "Phase") +
    ggtitle("PCA - After CC Regression")
  save_plot(p_pca_phase_after, "05_pca_by_phase_after.png", output_dir = args$output)

  # UMAP after regression
  seurat_obj <- RunUMAP(seurat_obj, dims = dims_use,
                         reduction = "pca_regressed",
                         reduction.name = "umap_after",
                         verbose = args$verbose)

  p_umap_phase_after <- DimPlot(seurat_obj, reduction = "umap_after", group.by = "Phase") +
    ggtitle("UMAP - After CC Regression")
  save_plot(p_umap_phase_after, "06_umap_by_phase_after.png", output_dir = args$output)

  # Comparison plot
  p_comparison <- (p_pca_phase + p_pca_phase_after) / (p_umap_phase + p_umap_phase_after)
  save_plot(p_comparison, "07_before_after_comparison.png",
            output_dir = args$output, width = 12, height = 10)

  # Use regressed reduction for downstream
  default_reduction <- "pca_regressed"
  default_umap <- "umap_after"
} else {
  default_reduction <- "pca"
  default_umap <- "umap_before"
}

# -----------------------------------------------------------------------------
# Clustering
# -----------------------------------------------------------------------------

log_message("Clustering...", log_env)

seurat_obj <- FindNeighbors(seurat_obj, reduction = default_reduction,
                             dims = dims_use, verbose = args$verbose)
seurat_obj <- FindClusters(seurat_obj, resolution = args$resolution, verbose = args$verbose)

n_clusters <- length(unique(seurat_obj$seurat_clusters))
log_message(sprintf("Found %d clusters", n_clusters), log_env)

# Cluster visualization
p_clusters <- DimPlot(seurat_obj, reduction = default_umap, label = TRUE) +
  ggtitle(sprintf("Clusters (res=%.2f)", args$resolution))
save_plot(p_clusters, "08_umap_clusters.png", output_dir = args$output)

# Cell cycle by cluster
cc_by_cluster <- as.data.frame(table(seurat_obj$seurat_clusters, seurat_obj$Phase))
colnames(cc_by_cluster) <- c("Cluster", "Phase", "Count")

p_cc_cluster <- ggplot(cc_by_cluster, aes(x = Cluster, y = Count, fill = Phase)) +
  geom_bar(stat = "identity", position = "fill") +
  theme_minimal() +
  labs(title = "Cell Cycle Phase by Cluster", y = "Proportion") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
save_plot(p_cc_cluster, "09_cc_phase_by_cluster.png", output_dir = args$output)

# Save cluster composition
cluster_cc_comp <- as.data.frame.matrix(table(seurat_obj$seurat_clusters, seurat_obj$Phase))
cluster_cc_comp$Cluster <- rownames(cluster_cc_comp)
write.csv(cluster_cc_comp, file.path(args$output, "cluster_cc_composition.csv"), row.names = FALSE)

# -----------------------------------------------------------------------------
# Save Results
# -----------------------------------------------------------------------------

log_message("Saving results...", log_env)

# Save cluster assignments
cluster_assignments <- data.frame(
  cell = colnames(seurat_obj),
  cluster = as.integer(as.character(seurat_obj$seurat_clusters))
)
write.csv(cluster_assignments, file.path(args$output, "cluster_assignments.csv"), row.names = FALSE)

save_seurat(seurat_obj, file.path(args$output, "seurat_cell_cycle.rds"))

# Summary statistics
summary_stats <- list(
  total_cells = ncol(seurat_obj),
  g1_cells = sum(seurat_obj$Phase == "G1"),
  s_cells = sum(seurat_obj$Phase == "S"),
  g2m_cells = sum(seurat_obj$Phase == "G2M"),
  g1_pct = round(100 * mean(seurat_obj$Phase == "G1"), 1),
  s_pct = round(100 * mean(seurat_obj$Phase == "S"), 1),
  g2m_pct = round(100 * mean(seurat_obj$Phase == "G2M"), 1),
  mean_s_score = round(mean(seurat_obj$S.Score), 4),
  mean_g2m_score = round(mean(seurat_obj$G2M.Score), 4),
  s_genes_used = length(s_genes_present),
  g2m_genes_used = length(g2m_genes_present),
  cc_regressed = args$regress_cc || args$regress_cc_difference
)

write.csv(as.data.frame(summary_stats),
          file.path(args$output, "analysis_summary.csv"), row.names = FALSE)

# Session info
writeLines(capture.output(sessionInfo()), file.path(args$output, "session_info.txt"))

# Parameters
params <- as.data.frame(t(as.data.frame(args)))
write.csv(params, file.path(args$output, "parameters.csv"))

print_summary(seurat_obj)

log_message("Cell cycle analysis complete!", log_env)
log_message(sprintf("Results saved to: %s", args$output), log_env)
