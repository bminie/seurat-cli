#!/usr/bin/env Rscript
#' =============================================================================
#' Imaging Spatial Transcriptomics (10x Xenium and Atera)
#' =============================================================================
#'
#' This script implements the Xenium section of Seurat's spatial analysis
#' vignette and the Atera tutorial (Seurat >= 6.0) for output folders from
#' 10x Genomics Xenium and Atera instruments.
#'
#' Workflow:
#'   1. Load counts, cell centroids, optional cell/nucleus outlines, and
#'      transcript positions for selected genes (LoadXenium / LoadAtera)
#'   2. QC: remove cells with too few transcripts; plot counts in space
#'   3. Normalize (SCTransform by default), PCA, UMAP, clustering
#'   4. Cluster markers, clusters and genes on the tissue image
#'   5. Optional zoomed-in view with cell outlines and transcripts
#'   6. Atera: plot the Azimuth cell type annotations shipped with the data
#'
#' Based on:
#'   https://satijalab.org/seurat/articles/spatial_vignette_2 (Xenium)
#'   https://satijalab.org/seurat/articles/atera_tutorial
#'
#' Usage:
#'   Rscript 08_spatial.R --input /path/to/xenium_outs --molecules Gad1,Sst \
#'     --output ./results
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
              help = "Xenium or Atera output folder (required)", metavar = "DIR"),
  make_option(c("-o", "--output"), type = "character", default = "output/spatial",
              help = "Output directory [default: %default]", metavar = "DIR"),
  make_option("--platform", type = "character", default = "auto",
              help = "auto, xenium or atera [default: %default]"),

  # Loading
  make_option("--molecules", type = "character", default = NULL,
              help = paste("Comma-separated genes whose transcript positions to load and plot.",
                           "Only these are read, which keeps memory down")),
  make_option("--segmentations", type = "character", default = "auto",
              help = paste("Outlines to load: cell, nucleus, none, or auto",
                           "(cell for Xenium, none for Atera) [default: %default]")),
  make_option("--qv_threshold", type = "double", default = 20,
              help = "Minimum transcript quality score for --molecules [default: %default]"),
  make_option("--bpcells_dir", type = "character", default = NULL,
              help = paste("Atera only: keep the counts on disk in this folder (BPCells).",
                           "Default: load into memory")),

  # QC and analysis
  make_option("--min_counts", type = "integer", default = 1,
              help = "Remove cells with fewer transcripts than this [default: %default]"),
  make_option("--normalization", type = "character", default = "SCT",
              help = "SCT or LogNormalize [default: %default]"),
  make_option("--npcs", type = "integer", default = 30,
              help = "Number of PCs to compute [default: %default]"),
  make_option("--dims_use", type = "character", default = "1:30",
              help = "PCs to use for neighbors and UMAP [default: %default]"),
  make_option("--resolution", type = "double", default = 0.3,
              help = "Clustering resolution [default: %default]"),
  make_option("--features", type = "character", default = NULL,
              help = "Comma-separated genes to plot on the tissue. Default: top cluster markers"),
  make_option("--crop", type = "character", default = NULL,
              help = paste("Zoomed view: xmin,xmax,ymin,ymax in the data's micron coordinates",
                           "(see cell_coordinates.csv)")),
  make_option("--no_markers", action = "store_true", default = FALSE,
              help = "Skip finding cluster markers"),

  # General
  make_option(c("-t", "--threads"), type = "integer", default = 1,
              help = "Number of threads [default: %default]"),
  make_option(c("-v", "--verbose"), action = "store_true", default = TRUE,
              help = "Verbose output"),
  make_option(c("-q", "--quiet"), action = "store_false", dest = "verbose",
              help = "Suppress messages"),
  make_option("--seed", type = "integer", default = 42,
              help = "Random seed [default: %default]")
)

parser <- OptionParser(
  option_list = option_list,
  description = "\nAnalyze 10x Xenium or Atera imaging spatial data: QC, clustering, and spatial plots.",
  epilogue = paste(
    "Examples:",
    "  Rscript 08_spatial.R --input xenium_outs --molecules Gad1,Sst,Pvalb,Gfap \\",
    "    --crop 1200,2900,3750,4550 --output ./results",
    "  Rscript 08_spatial.R --input atera_outs --platform atera --threads 8 --output ./results",
    sep = "\n"
  )
)
args <- parse_args(parser)

# -----------------------------------------------------------------------------
# Validate Arguments
# -----------------------------------------------------------------------------

if (is.null(args$input)) {
  print_help(parser)
  stop("Error: --input is required")
}
if (!dir.exists(args$input)) {
  stop(sprintf("Input folder not found: %s", args$input))
}
if (!args$platform %in% c("auto", "xenium", "atera")) {
  stop(sprintf("Unknown --platform '%s'. Choose auto, xenium or atera", args$platform))
}
if (!args$segmentations %in% c("auto", "cell", "nucleus", "none")) {
  stop(sprintf("Unknown --segmentations '%s'. Choose auto, cell, nucleus or none",
               args$segmentations))
}
if (!args$normalization %in% c("SCT", "LogNormalize")) {
  stop(sprintf("Unknown --normalization '%s'. Choose SCT or LogNormalize", args$normalization))
}

parse_dims <- function(dims_str) {
  if (grepl(":", dims_str)) {
    parts <- as.integer(strsplit(dims_str, ":")[[1]])
    seq(parts[1], parts[2])
  } else {
    as.integer(strsplit(dims_str, ",")[[1]])
  }
}
split_list <- function(x) if (is.null(x)) character() else trimws(strsplit(x, ",")[[1]])

dims_use <- parse_dims(args$dims_use)
if (max(dims_use) > args$npcs) {
  stop(sprintf("--dims_use goes up to %d but --npcs is %d", max(dims_use), args$npcs))
}
molecule_genes <- split_list(args$molecules)

crop <- NULL
if (!is.null(args$crop)) {
  crop <- as.numeric(split_list(args$crop))
  if (length(crop) != 4 || anyNA(crop) || crop[1] >= crop[2] || crop[3] >= crop[4]) {
    stop("--crop must be xmin,xmax,ymin,ymax with xmin < xmax and ymin < ymax")
  }
}

# Xenium folders have experiment.xenium and a count matrix in .h5 or MEX
# format; Atera folders have experiment.spatial and zarr archives only
detect_platform <- function(dir) {
  has <- function(f) file.exists(file.path(dir, f))
  if (has("experiment.xenium") || has("cell_feature_matrix.h5") ||
      dir.exists(file.path(dir, "cell_feature_matrix"))) {
    return("xenium")
  }
  if (has("experiment.spatial") || has("csc_cell_feature_matrix.zarr.zip") ||
      has("cell_feature_matrix.zarr.zip")) {
    return("atera")
  }
  stop(sprintf(paste("Could not recognize %s as Xenium or Atera output",
                     "(no experiment.xenium, experiment.spatial or cell feature matrix).",
                     "Set --platform if the files were renamed."), dir))
}
platform <- if (args$platform == "auto") detect_platform(args$input) else args$platform
assay_name <- if (platform == "xenium") "Xenium" else "Atera"
segmentations <- switch(args$segmentations,
                        auto = if (platform == "xenium") "cell" else NULL,
                        none = NULL,
                        args$segmentations)

set.seed(args$seed)

# -----------------------------------------------------------------------------
# Setup
# -----------------------------------------------------------------------------

if (!dir.exists(args$output)) {
  dir.create(args$output, recursive = TRUE)
}

log_file <- file.path(args$output, "analysis.log")
log_env <- init_logging(log_file = log_file, verbose = args$verbose)

log_message("Starting imaging spatial analysis", log_env)
log_message(sprintf("Input: %s (%s)", args$input, platform), log_env)
log_message(sprintf("Output directory: %s", args$output), log_env)

log_message("Loading required packages...", log_env)
load_seurat_deps(quiet = !args$verbose)
setup_threads(args$threads, log_env)

if (platform == "atera") {
  missing_pkgs <- c("blosc", "jsonlite", "data.table")[
    !vapply(c("blosc", "jsonlite", "data.table"), requireNamespace, logical(1), quietly = TRUE)]
  if (length(missing_pkgs) > 0) {
    stop(sprintf("Reading Atera data needs: %s. Install with: Rscript setup.R --all",
                 paste(missing_pkgs, collapse = ", ")), call. = FALSE)
  }
}
if (!is.null(crop) && !requireNamespace("sf", quietly = TRUE)) {
  stop("--crop needs the sf package: install.packages('sf') or Rscript setup.R --all",
       call. = FALSE)
}
if (platform == "xenium" && file.exists(file.path(args$input, "transcripts.parquet")) &&
    !file.exists(file.path(args$input, "transcripts.csv.gz")) &&
    length(molecule_genes) > 0 && !requireNamespace("arrow", quietly = TRUE)) {
  stop("This Xenium folder stores transcripts as Parquet only; install the arrow package")
}

# -----------------------------------------------------------------------------
# Load Data
# -----------------------------------------------------------------------------

log_message("Loading spatial data...", log_env)
load_start <- Sys.time()

if (platform == "xenium") {
  seurat_obj <- LoadXenium(
    data.dir = args$input,
    fov = "fov",
    assay = assay_name,
    mols.qv.threshold = args$qv_threshold,
    cell.centroids = TRUE,
    molecule.coordinates = length(molecule_genes) > 0,
    segmentations = segmentations,
    genes = if (length(molecule_genes) > 0) molecule_genes else NULL
  )
} else {
  # The Atera reader parallelizes over future workers
  if (args$threads > 1) {
    future::plan(future::multisession, workers = args$threads)
  }
  seurat_obj <- LoadAtera(
    data.dir = args$input,
    fov = "fov",
    assay = assay_name,
    mols.qv.threshold = args$qv_threshold,
    cell.centroids = TRUE,
    molecule.coordinates = length(molecule_genes) > 0,
    segmentations = segmentations,
    cell.annotations = TRUE,
    genes = if (length(molecule_genes) > 0) molecule_genes else NULL,
    bpcells.dir = if (is.null(args$bpcells_dir)) FALSE else args$bpcells_dir
  )
  if (args$threads > 1) future::plan(future::sequential)
}

log_message(sprintf("Loaded %d cells x %d genes in %.0f s", ncol(seurat_obj),
                    nrow(seurat_obj[[assay_name]]),
                    as.numeric(difftime(Sys.time(), load_start, units = "secs"))), log_env)
DefaultAssay(seurat_obj) <- assay_name

# Transcript positions use assay feature names ("_" becomes "-")
molecule_genes <- gsub("_", "-", molecule_genes)
unknown_mols <- setdiff(molecule_genes, rownames(seurat_obj[[assay_name]]))
if (length(unknown_mols) > 0) {
  log_message(sprintf("WARNING: --molecules genes not in the panel: %s",
                      paste(unknown_mols, collapse = ", ")), log_env)
  molecule_genes <- setdiff(molecule_genes, unknown_mols)
}

# -----------------------------------------------------------------------------
# Quality Control
# -----------------------------------------------------------------------------

count_col <- paste0("nCount_", assay_name)
feature_col <- paste0("nFeature_", assay_name)

p_qc <- VlnPlot(seurat_obj, features = c(feature_col, count_col), ncol = 2, pt.size = 0)
save_plot(p_qc, "01_qc_violin.png", output_dir = args$output, width = 10, height = 5)

n_before <- ncol(seurat_obj)
seurat_obj <- seurat_obj[, seurat_obj[[count_col]][, 1] >= args$min_counts]
log_message(sprintf("Cells with >= %d transcripts: %d of %d", args$min_counts,
                    ncol(seurat_obj), n_before), log_env)
if (ncol(seurat_obj) < 50) {
  stop(sprintf("Only %d cells pass --min_counts %d; too few to analyze",
               ncol(seurat_obj), args$min_counts))
}

p_spatial_qc <- ImageFeaturePlot(seurat_obj, fov = "fov", features = count_col, size = 0.5) +
  ggtitle("Transcripts per cell") |
  ImageFeaturePlot(seurat_obj, fov = "fov", features = feature_col, size = 0.5) +
  ggtitle("Genes per cell")
save_plot(p_spatial_qc, "02_spatial_qc.png", output_dir = args$output, width = 14, height = 7)

if (length(molecule_genes) > 0) {
  p_mols <- ImageDimPlot(seurat_obj, fov = "fov", boundaries = "centroids", cols = "grey80",
                         size = 0.1, molecules = molecule_genes, nmols = 20000,
                         mols.size = 0.1, axes = TRUE) +
    ggtitle("Transcript positions")
  save_plot(p_mols, "03_molecules.png", output_dir = args$output, width = 10, height = 8)
}

# -----------------------------------------------------------------------------
# Normalization, Dimensional Reduction, Clustering
# -----------------------------------------------------------------------------

n_genes <- nrow(seurat_obj[[assay_name]])
npcs <- min(args$npcs, n_genes - 1, ncol(seurat_obj) - 1)
dims_use <- dims_use[dims_use <= npcs]

if (args$normalization == "SCT") {
  log_message("Running SCTransform...", log_env)
  seurat_obj <- SCTransform(seurat_obj, assay = assay_name, verbose = args$verbose)
} else {
  log_message("Running NormalizeData and ScaleData...", log_env)
  seurat_obj <- NormalizeData(seurat_obj, verbose = args$verbose)
  seurat_obj <- ScaleData(seurat_obj, features = rownames(seurat_obj), verbose = args$verbose)
}

# Imaging panels are small and targeted, so all genes are used for PCA
log_message(sprintf("Running PCA (%d PCs), UMAP and clustering...", npcs), log_env)
seurat_obj <- RunPCA(seurat_obj, npcs = npcs, features = rownames(seurat_obj),
                     verbose = args$verbose)
seurat_obj <- RunUMAP(seurat_obj, dims = dims_use, verbose = args$verbose)
seurat_obj <- FindNeighbors(seurat_obj, reduction = "pca", dims = dims_use,
                            verbose = args$verbose)
seurat_obj <- FindClusters(seurat_obj, resolution = args$resolution, verbose = args$verbose)

n_clusters <- length(unique(seurat_obj$seurat_clusters))
log_message(sprintf("Found %d clusters", n_clusters), log_env)

p_umap <- DimPlot(seurat_obj, reduction = "umap", label = TRUE, repel = TRUE) +
  ggtitle("Clusters")
save_plot(p_umap, "04_umap_clusters.png", output_dir = args$output, width = 9, height = 7)

p_space <- ImageDimPlot(seurat_obj, fov = "fov", cols = "polychrome", size = 0.75) +
  ggtitle("Clusters in space")
save_plot(p_space, "05_spatial_clusters.png", output_dir = args$output, width = 10, height = 8)

# -----------------------------------------------------------------------------
# Markers and Genes in Space
# -----------------------------------------------------------------------------

markers <- NULL
if (!args$no_markers) {
  log_message("Finding cluster markers...", log_env)
  if (args$normalization == "SCT") {
    seurat_obj <- PrepSCTFindMarkers(seurat_obj, verbose = FALSE)
  }
  markers <- FindAllMarkers(seurat_obj, only.pos = TRUE, min.pct = 0.1,
                            logfc.threshold = 0.25, verbose = FALSE)
  if (nrow(markers) > 0) {
    write.csv(markers, file.path(args$output, "cluster_markers_all.csv"), row.names = FALSE)
  }
}

plot_features <- split_list(args$features)
if (length(plot_features) == 0) {
  # Default: the strongest marker of each of the first clusters
  plot_features <- if (!is.null(markers) && nrow(markers) > 0) {
    top <- markers[order(markers$cluster, -markers$avg_log2FC), ]
    head(unique(top$gene[!duplicated(top$cluster)]), 6)
  } else {
    head(rownames(seurat_obj[[assay_name]]), 4)
  }
}
plot_features <- intersect(plot_features, rownames(seurat_obj))

if (length(plot_features) > 0) {
  p_feat_space <- ImageFeaturePlot(seurat_obj, fov = "fov", features = plot_features,
                                   max.cutoff = "q95", size = 0.5, cols = c("white", "red"))
  save_plot(p_feat_space, "06_spatial_features.png", output_dir = args$output,
            width = 12, height = 4 * ceiling(length(plot_features) / 3))
  p_feat_umap <- FeaturePlot(seurat_obj, features = plot_features, ncol = 3)
  save_plot(p_feat_umap, "07_umap_features.png", output_dir = args$output,
            width = 12, height = 4 * ceiling(length(plot_features) / 3))
}

# -----------------------------------------------------------------------------
# Zoomed View
# -----------------------------------------------------------------------------

if (!is.null(crop)) {
  log_message(sprintf("Zoomed view: x %g-%g, y %g-%g", crop[1], crop[2], crop[3], crop[4]),
              log_env)
  zoom <- Crop(seurat_obj[["fov"]], x = crop[1:2], y = crop[3:4], coords = "tissue")
  seurat_obj[["zoom"]] <- zoom
  if (!is.null(segmentations)) {
    DefaultBoundary(seurat_obj[["zoom"]]) <- "segmentation"
  }
  n_zoom <- length(Cells(seurat_obj[["zoom"]]))
  if (n_zoom == 0) {
    log_message("WARNING: no cells inside --crop; skipping the zoomed view", log_env)
  } else {
    p_zoom <- ImageDimPlot(seurat_obj, fov = "zoom", axes = TRUE, border.color = "white",
                           border.size = 0.1, cols = "polychrome", coord.fixed = FALSE,
                           molecules = if (length(molecule_genes) > 0) molecule_genes else NULL,
                           nmols = 10000) +
      ggtitle(sprintf("Zoomed view (%d cells)", n_zoom))
    save_plot(p_zoom, "08_zoom.png", output_dir = args$output, width = 10, height = 8)
  }
}

# -----------------------------------------------------------------------------
# Atera Cell Type Annotations
# -----------------------------------------------------------------------------

azimuth_cols <- grep("^azimuth_", colnames(seurat_obj[[]]), value = TRUE)
if (length(azimuth_cols) > 0) {
  label_col <- if ("azimuth_broad" %in% azimuth_cols) "azimuth_broad" else azimuth_cols[1]
  log_message(sprintf("Plotting cell type annotations (%s)", label_col), log_env)
  p_anno <- ImageDimPlot(seurat_obj, fov = "fov", group.by = label_col, size = 0.3) +
    ggtitle(label_col)
  save_plot(p_anno, "09_spatial_annotations.png", output_dir = args$output,
            width = 10, height = 8)
}

# -----------------------------------------------------------------------------
# Save Results
# -----------------------------------------------------------------------------

log_message("Saving results...", log_env)

save_seurat(seurat_obj, file.path(args$output, "seurat_spatial.rds"))

centroids <- GetTissueCoordinates(seurat_obj[["fov"]], which = "centroids")
cell_df <- data.frame(cell = colnames(seurat_obj), cluster = seurat_obj$seurat_clusters)
cell_df <- merge(cell_df, centroids[, c("cell", "x", "y")], by = "cell", all.x = TRUE,
                 sort = FALSE)
write.csv(cell_df, file.path(args$output, "cell_coordinates.csv"), row.names = FALSE)
write.csv(cell_df[, c("cell", "cluster")], file.path(args$output, "cluster_assignments.csv"),
          row.names = FALSE)

print_summary(seurat_obj)
writeLines(capture.output(sessionInfo()), file.path(args$output, "session_info.txt"))
params <- as.data.frame(t(as.data.frame(args)))
write.csv(params, file.path(args$output, "parameters.csv"))

log_message("Spatial analysis complete!", log_env)
log_message(sprintf("Results saved to: %s", args$output), log_env)
