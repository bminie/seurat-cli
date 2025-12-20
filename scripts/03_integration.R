#!/usr/bin/env Rscript
#' =============================================================================
#' scRNA-seq Data Integration Vignette
#' =============================================================================
#'
#' This script implements the Seurat integration vignettes for combining
#' multiple single-cell datasets:
#'   - Standard integration (CCA-based anchors)
#'   - RPCA-based integration (faster, more conservative)
#'   - Harmony integration (optional)
#'
#' Supports:
#'   - Integration across conditions/batches
#'   - Integration across technologies
#'   - Integration across species (with appropriate gene mapping)
#'
#' Based on: https://satijalab.org/seurat/articles/integration_introduction
#'           https://satijalab.org/seurat/articles/seurat5_integration
#'
#' Usage:
#'   Rscript 03_integration.R --input_list samples.txt --output ./integrated
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
  make_option(c("-i", "--input_list"), type = "character", default = NULL,
              help = "File with list of input paths (one per line) or comma-separated paths",
              metavar = "PATH"),
  make_option(c("-o", "--output"), type = "character", default = "output/integration",
              help = "Output directory [default: %default]", metavar = "DIR"),
  make_option(c("-n", "--name"), type = "character", default = "Integrated",
              help = "Project name [default: %default]", metavar = "NAME"),
  make_option(c("-f", "--format"), type = "character", default = "auto",
              help = "Input format [default: %default]"),

  # Sample info
  make_option("--sample_names", type = "character", default = NULL,
              help = "Comma-separated sample names (in same order as inputs)"),
  make_option("--conditions", type = "character", default = NULL,
              help = "Comma-separated condition labels for samples"),

  # Species
  make_option(c("-s", "--species"), type = "character", default = "human",
              help = "Species (human, mouse) [default: %default]"),

  # Pre-processing
  make_option("--min_cells", type = "integer", default = 3,
              help = "Minimum cells per feature [default: %default]"),
  make_option("--min_features", type = "integer", default = 200,
              help = "Minimum features per cell [default: %default]"),
  make_option("--max_features", type = "integer", default = 5000,
              help = "Maximum features per cell [default: %default]"),
  make_option("--max_mt", type = "double", default = 10,
              help = "Maximum percent mitochondrial [default: %default]"),

  # Normalization
  make_option("--normalization", type = "character", default = "LogNormalize",
              help = "Normalization method (LogNormalize, SCT) [default: %default]"),
  make_option("--n_variable_features", type = "integer", default = 2000,
              help = "Number of variable features per sample [default: %default]"),

  # Integration method
  make_option("--method", type = "character", default = "CCAIntegration",
              help = "Integration method: CCAIntegration, RPCAIntegration, HarmonyIntegration, FastMNNIntegration [default: %default]"),
  make_option("--reference", type = "integer", default = NULL,
              help = "Index of reference sample (1-based) for RPCA. NULL = no reference"),
  make_option("--k_anchor", type = "integer", default = 5,
              help = "Number of anchors for integration [default: %default]"),

  # Dimensional Reduction
  make_option("--n_pcs", type = "integer", default = 50,
              help = "Number of PCs to compute [default: %default]"),
  make_option("--dims_use", type = "character", default = "1:30",
              help = "PCs to use for integration/UMAP [default: %default]"),

  # Clustering
  make_option("--resolution", type = "double", default = 0.5,
              help = "Clustering resolution [default: %default]"),
  make_option("--resolutions", type = "character", default = NULL,
              help = "Multiple resolutions (comma-separated) for testing"),

  # Markers
  make_option("--find_markers", action = "store_true", default = TRUE,
              help = "Find cluster markers"),
  make_option("--no_markers", action = "store_false", dest = "find_markers",
              help = "Skip marker finding"),
  make_option("--find_conserved", action = "store_true", default = FALSE,
              help = "Find conserved markers across conditions"),

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
              help = "Run with demo ifnb dataset from SeuratData")
)

parser <- OptionParser(
  usage = "%prog [options]",
  option_list = option_list,
  description = "Integrate multiple scRNA-seq datasets using Seurat v5"
)

args <- parse_args(parser)

# -----------------------------------------------------------------------------
# Validate Arguments
# -----------------------------------------------------------------------------

if (!args$demo && is.null(args$input_list)) {
  print_help(parser)
  stop("Error: --input_list is required (or use --demo)")
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

# Parse resolutions
resolutions <- args$resolution
if (!is.null(args$resolutions)) {
  resolutions <- as.numeric(strsplit(args$resolutions, ",")[[1]])
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

log_message("Starting integration workflow", log_env)
log_message(sprintf("Integration method: %s", args$method), log_env)
log_message(sprintf("Normalization: %s", args$normalization), log_env)

if (args$threads > 1) {
  library(future)
  plan(multisession, workers = args$threads)
  log_message(sprintf("Using %d threads", args$threads), log_env)
}

# -----------------------------------------------------------------------------
# Load Required Packages
# -----------------------------------------------------------------------------

log_message("Loading required packages...", log_env)
load_seurat_deps(include_integration = TRUE, quiet = !args$verbose)

# Load additional packages for specific methods
if (args$method == "HarmonyIntegration") {
  install_and_load(packages = "harmony", quiet = !args$verbose)
}

# -----------------------------------------------------------------------------
# Load Data
# -----------------------------------------------------------------------------

if (args$demo) {
  log_message("Loading demo ifnb dataset (stimulated vs control PBMCs)...", log_env)

  install_and_load(packages = "SeuratData", quiet = !args$verbose)

  if (!"ifnb" %in% rownames(installed.packages())) {
    InstallData("ifnb")
  }
  data("ifnb")

  # The ifnb dataset is already a merged object with 'stim' metadata
  seurat_obj <- ifnb
  rm(ifnb)

  # Update to current Seurat version format (v5 compatibility)
  seurat_obj <- UpdateSeuratObject(seurat_obj)

  # For demo, split by condition
  seurat_list <- SplitObject(seurat_obj, split.by = "stim")
  sample_names <- names(seurat_list)

} else {
  # Parse input list
  if (file.exists(args$input_list)) {
    input_paths <- readLines(args$input_list)
    input_paths <- input_paths[nchar(trimws(input_paths)) > 0]
  } else {
    input_paths <- trimws(strsplit(args$input_list, ",")[[1]])
  }

  # Parse sample names
  if (!is.null(args$sample_names)) {
    sample_names <- trimws(strsplit(args$sample_names, ",")[[1]])
    if (length(sample_names) != length(input_paths)) {
      stop("Number of sample names must match number of inputs")
    }
  } else {
    sample_names <- paste0("Sample", seq_along(input_paths))
  }

  log_message(sprintf("Loading %d samples...", length(input_paths)), log_env)

  seurat_list <- lapply(seq_along(input_paths), function(i) {
    log_message(sprintf("  Loading %s: %s", sample_names[i], input_paths[i]), log_env)

    obj <- load_seurat_data(
      input_path = input_paths[i],
      input_format = args$format,
      sample_name = sample_names[i],
      min_cells = args$min_cells,
      min_features = args$min_features
    )

    # Add sample identifier to metadata
    obj$sample <- sample_names[i]
    obj$orig.ident <- sample_names[i]

    return(obj)
  })

  names(seurat_list) <- sample_names
}

# Add conditions if provided
if (!is.null(args$conditions)) {
  conditions <- trimws(strsplit(args$conditions, ",")[[1]])
  if (length(conditions) == length(seurat_list)) {
    for (i in seq_along(seurat_list)) {
      seurat_list[[i]]$condition <- conditions[i]
    }
  }
}

log_message(sprintf("Loaded %d samples", length(seurat_list)), log_env)
for (nm in names(seurat_list)) {
  log_message(sprintf("  %s: %d cells", nm, ncol(seurat_list[[nm]])), log_env)
}

# -----------------------------------------------------------------------------
# QC and Filtering
# -----------------------------------------------------------------------------

log_message("Performing QC and filtering...", log_env)

seurat_list <- lapply(seurat_list, function(obj) {
  obj <- calculate_qc_metrics(obj, species = args$species)

  obj <- subset(
    obj,
    subset = nFeature_RNA > args$min_features &
             nFeature_RNA < args$max_features &
             percent.mt < args$max_mt
  )

  return(obj)
})

for (nm in names(seurat_list)) {
  log_message(sprintf("  %s after QC: %d cells", nm, ncol(seurat_list[[nm]])), log_env)
}

# -----------------------------------------------------------------------------
# Merge and Normalize
# -----------------------------------------------------------------------------

log_message("Merging datasets...", log_env)

seurat_obj <- merge(
  seurat_list[[1]],
  y = seurat_list[-1],
  add.cell.ids = names(seurat_list),
  project = args$name
)

log_message(sprintf("Merged object: %d cells, %d features",
                    ncol(seurat_obj), nrow(seurat_obj)), log_env)

# Split layers by sample for Seurat v5 integration
seurat_obj[["RNA"]] <- split(seurat_obj[["RNA"]], f = seurat_obj$orig.ident)

# QC visualization before integration
p_qc_pre <- VlnPlot(seurat_obj, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"),
                    group.by = "orig.ident", pt.size = 0, ncol = 3)
save_plot(p_qc_pre, "01_qc_by_sample.png", output_dir = args$output, width = 15, height = 5)

# Normalize
log_message(sprintf("Normalizing (%s)...", args$normalization), log_env)

if (args$normalization == "SCT") {
  seurat_obj <- SCTransform(
    seurat_obj,
    variable.features.n = args$n_variable_features,
    verbose = args$verbose
  )
} else {
  seurat_obj <- NormalizeData(seurat_obj, verbose = args$verbose)
  seurat_obj <- FindVariableFeatures(seurat_obj,
                                      nfeatures = args$n_variable_features,
                                      verbose = args$verbose)
  seurat_obj <- ScaleData(seurat_obj, verbose = args$verbose)
}

# -----------------------------------------------------------------------------
# Run PCA (before integration)
# -----------------------------------------------------------------------------

log_message("Running PCA...", log_env)

seurat_obj <- RunPCA(seurat_obj, npcs = args$n_pcs, verbose = args$verbose)

# UMAP before integration
seurat_obj <- RunUMAP(seurat_obj, dims = dims_use, reduction.name = "umap.unintegrated",
                       verbose = args$verbose)

p_pre_integration <- DimPlot(seurat_obj, reduction = "umap.unintegrated",
                              group.by = "orig.ident") +
  ggtitle("Before Integration")

save_plot(p_pre_integration, "02_umap_before_integration.png",
          output_dir = args$output, width = 10, height = 8)

# -----------------------------------------------------------------------------
# Integration
# -----------------------------------------------------------------------------

log_message(sprintf("Running integration (%s)...", args$method), log_env)

# Seurat v5 integration
seurat_obj <- IntegrateLayers(
  object = seurat_obj,
  method = get(args$method),
  orig.reduction = "pca",
  new.reduction = "integrated",
  verbose = args$verbose
)

# Re-join layers after integration
seurat_obj <- JoinLayers(seurat_obj)

log_message("Integration complete", log_env)

# -----------------------------------------------------------------------------
# Post-Integration Analysis
# -----------------------------------------------------------------------------

log_message("Running UMAP on integrated data...", log_env)

seurat_obj <- RunUMAP(
  seurat_obj,
  reduction = "integrated",
  dims = dims_use,
  reduction.name = "umap.integrated",
  verbose = args$verbose
)

# Comparison plot
p_post_integration <- DimPlot(seurat_obj, reduction = "umap.integrated",
                               group.by = "orig.ident") +
  ggtitle("After Integration")

p_comparison <- p_pre_integration + p_post_integration
save_plot(p_comparison, "03_integration_comparison.png",
          output_dir = args$output, width = 16, height = 6)

# Split by sample
p_split <- DimPlot(seurat_obj, reduction = "umap.integrated",
                    split.by = "orig.ident", ncol = min(length(sample_names), 4))
save_plot(p_split, "04_umap_split_by_sample.png",
          output_dir = args$output,
          width = min(length(sample_names), 4) * 5, height = 5)

# -----------------------------------------------------------------------------
# Clustering
# -----------------------------------------------------------------------------

log_message("Clustering integrated data...", log_env)

# Use explicit graph name for consistency
graph_name <- "integrated_snn"
seurat_obj <- FindNeighbors(seurat_obj, reduction = "integrated",
                             dims = dims_use,
                             graph.name = c(paste0(graph_name, "_nn"), graph_name),
                             verbose = args$verbose)

# Run clustering at specified resolution(s)
for (res in resolutions) {
  log_message(sprintf("  Clustering at resolution %.2f", res), log_env)

  seurat_obj <- FindClusters(seurat_obj, resolution = res,
                              graph.name = graph_name,
                              verbose = args$verbose)
}

# Use primary resolution for default clustering
cluster_col <- paste0(graph_name, "_res.", args$resolution)
Idents(seurat_obj) <- seurat_obj[[cluster_col]][[1]]
seurat_obj$seurat_clusters <- Idents(seurat_obj)

n_clusters <- length(unique(seurat_obj$seurat_clusters))
log_message(sprintf("Found %d clusters at resolution %.2f",
                    n_clusters, args$resolution), log_env)

# Cluster visualization
p_clusters <- DimPlot(seurat_obj, reduction = "umap.integrated",
                       label = TRUE, label.size = 4) +
  ggtitle(sprintf("Integrated Clusters (res=%.2f)", args$resolution))

save_plot(p_clusters, "05_umap_clusters.png",
          output_dir = args$output, width = 10, height = 8)

# Clusters by sample
p_clusters_split <- DimPlot(seurat_obj, reduction = "umap.integrated",
                             split.by = "orig.ident", label = TRUE, label.size = 3,
                             ncol = min(length(sample_names), 4))
save_plot(p_clusters_split, "06_clusters_by_sample.png",
          output_dir = args$output,
          width = min(length(sample_names), 4) * 5, height = 5)

# Cluster composition by sample
cluster_comp <- table(seurat_obj$seurat_clusters, seurat_obj$orig.ident)
write.csv(as.data.frame.matrix(cluster_comp),
          file.path(args$output, "cluster_composition.csv"))

# -----------------------------------------------------------------------------
# Marker Finding
# -----------------------------------------------------------------------------

if (args$find_markers) {
  log_message("Finding cluster markers...", log_env)

  if (args$normalization == "SCT") {
    seurat_obj <- PrepSCTFindMarkers(seurat_obj, verbose = args$verbose)
  }

  markers <- FindAllMarkers(
    seurat_obj,
    only.pos = TRUE,
    min.pct = 0.25,
    logfc.threshold = 0.25,
    verbose = args$verbose
  )

  # Check if markers were found
  if (is.null(markers) || nrow(markers) == 0 || !"cluster" %in% colnames(markers)) {
    log_message("WARNING: No markers found. Skipping marker analysis.", log_env)
    # Create empty marker files
    write.csv(data.frame(), file.path(args$output, "cluster_markers_all.csv"), row.names = FALSE)
    write.csv(data.frame(), file.path(args$output, "cluster_markers_top10.csv"), row.names = FALSE)
  } else {
    write.csv(markers, file.path(args$output, "cluster_markers_all.csv"), row.names = FALSE)

    top_markers <- markers |>
      group_by(cluster) |>
      slice_max(n = 10, order_by = avg_log2FC)

    write.csv(top_markers, file.path(args$output, "cluster_markers_top10.csv"), row.names = FALSE)

    log_message(sprintf("Found %d markers", nrow(markers)), log_env)

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
  }
}

# Conserved markers (across conditions)
if (args$find_conserved && "condition" %in% colnames(seurat_obj@meta.data)) {
  log_message("Finding conserved markers across conditions...", log_env)

  conserved_markers_list <- list()

  for (cluster_id in levels(seurat_obj$seurat_clusters)) {
    tryCatch({
      conserved <- FindConservedMarkers(
        seurat_obj,
        ident.1 = cluster_id,
        grouping.var = "condition",
        verbose = args$verbose
      )
      conserved$cluster <- cluster_id
      conserved_markers_list[[cluster_id]] <- conserved
    }, error = function(e) {
      log_message(sprintf("  Warning: Could not find conserved markers for cluster %s",
                          cluster_id), log_env)
    })
  }

  if (length(conserved_markers_list) > 0) {
    conserved_markers <- do.call(rbind, conserved_markers_list)
    write.csv(conserved_markers,
              file.path(args$output, "conserved_markers.csv"))
    log_message("Conserved markers saved", log_env)
  }
}

# -----------------------------------------------------------------------------
# Save Results
# -----------------------------------------------------------------------------

log_message("Saving results...", log_env)

save_seurat(seurat_obj, file.path(args$output, "seurat_integrated.rds"))

# Cluster assignments
cluster_df <- data.frame(
  cell = colnames(seurat_obj),
  cluster = seurat_obj$seurat_clusters,
  sample = seurat_obj$orig.ident
)
if ("condition" %in% colnames(seurat_obj@meta.data)) {
  cluster_df$condition <- seurat_obj$condition
}
write.csv(cluster_df, file.path(args$output, "cluster_assignments.csv"), row.names = FALSE)

# UMAP coordinates
umap_coords <- as.data.frame(Embeddings(seurat_obj, "umap.integrated"))
umap_coords$cell <- rownames(umap_coords)
umap_coords$cluster <- seurat_obj$seurat_clusters
umap_coords$sample <- seurat_obj$orig.ident
write.csv(umap_coords, file.path(args$output, "umap_coordinates.csv"), row.names = FALSE)

# Summary
print_summary(seurat_obj)
writeLines(capture.output(sessionInfo()), file.path(args$output, "session_info.txt"))

# Parameters
params <- as.data.frame(t(as.data.frame(args)))
write.csv(params, file.path(args$output, "parameters.csv"))

log_message("Integration analysis complete!", log_env)
log_message(sprintf("Results saved to: %s", args$output), log_env)
