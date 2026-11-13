#!/usr/bin/env Rscript
#' =============================================================================
#' Sample-Level Analysis (scSLIDE)
#' =============================================================================
#'
#' This script implements Seurat's scSLIDE vignette (Single-cell Sample-Level
#' Inference using Density Estimation, Seurat >= 6.0): it summarizes each
#' sample (e.g. donor) by how its cells are spread over a set of landmark
#' cells, finds the main axes of variation between samples, and tests which
#' genes change along them.
#'
#' Workflow:
#'   1. PrepareSampleObject: sketch a training set balanced across samples,
#'      pick landmark cells, learn a response-oriented PLS embedding, and
#'      connect every cell to its nearest landmarks (WNN of a cell-type
#'      embedding and the PLS embedding)
#'   2. GenerateSampleObject: count each sample's cells near each landmark
#'      (landmark x sample matrix, chi-squared normalized)
#'   3. RunDiffusionMap: sample-level trajectory components (DC1, DC2, ...)
#'   4. TrajDETest: pseudobulk differential expression along each component,
#'      per cell type
#'   5. PlotAnnotatedHeatmap: top trajectory genes across samples
#'
#' Differences from the vignette:
#'   - The cell-type embedding is an existing reduction (--reduction, e.g.
#'     the integrated reduction from 03_integration.R) instead of the cloud
#'     Azimuth annotation, which uploads the data to a web service.
#'   - The demo uses the ifnb dataset (8 donors x control/stimulated = 16
#'     samples) instead of the COVID-19 dataset. Donor IDs are downloaded once
#'     from the demuxlet calls published with the data (~5 MB).
#'
#' Based on: https://satijalab.org/seurat/articles/scSLIDE_vignette
#'
#' Usage:
#'   Rscript 07_sample_level.R --input seurat.rds --sample_col donor \
#'     --response disease --celltype_col cell_type --reduction integrated \
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
              help = "Seurat object (.rds) with normalized data from several samples (required)",
              metavar = "PATH"),
  make_option(c("-o", "--output"), type = "character", default = "output/sample_level",
              help = "Output directory [default: %default]", metavar = "DIR"),
  make_option("--assay", type = "character", default = "RNA",
              help = "Assay with the expression data [default: %default]"),

  # Sample design
  make_option("--sample_col", type = "character", default = NULL,
              help = "Metadata column with sample IDs, e.g. donor (required)"),
  make_option("--response", type = "character", default = NULL,
              help = "Sample-level metadata column the PLS embedding should capture, e.g. disease (required)"),
  make_option("--celltype_col", type = "character", default = NULL,
              help = paste("Metadata column with cell types. Used to find response-related genes",
                           "within cell types, to name landmarks, and to run trajectory DE per",
                           "cell type. Without it, DE uses all cells of each sample")),

  # Cell-level embedding and neighbor graph
  make_option("--reduction", type = "character", default = "pca",
              help = paste("Cell-type embedding combined with the PLS embedding.",
                           "pca is computed if missing [default: %default]")),
  make_option("--dims", type = "character", default = "1:30",
              help = "Dimensions of --reduction to use [default: %default]"),
  make_option("--modality", type = "character", default = "multi",
              help = paste("multi: neighbors from --reduction and PLS combined (WNN);",
                           "single: PLS embedding only [default: %default]")),
  make_option("--pls_method", type = "character", default = "cppls",
              help = "PLS algorithm: cppls, plsr or spls [default: %default]"),
  make_option("--ncomp", type = "integer", default = 10,
              help = "Number of PLS components [default: %default]"),
  make_option("--ncells_per_sample", type = "integer", default = 1000,
              help = "Cells per sample in the PLS training set [default: %default]"),
  make_option("--ncells_landmark", type = "integer", default = 2000,
              help = "Number of landmark cells [default: %default]"),
  make_option("--k_nn", type = "integer", default = 5,
              help = "Nearest landmarks per cell [default: %default]"),
  make_option("--no_cor_genes", action = "store_true", default = FALSE,
              help = "Don't add response-correlated genes to the variable features"),

  # Sample-level trajectory and DE
  make_option("--n_dcs", type = "integer", default = 10,
              help = "Diffusion map components to compute [default: %default]"),
  make_option("--traj_dcs", type = "character", default = "1",
              help = "Comma-separated diffusion components to test for DE [default: %default]"),
  make_option("--de_celltypes", type = "character", default = "all",
              help = "Comma-separated cell types for trajectory DE, or 'all' [default: %default]"),
  make_option("--latent_vars", type = "character", default = NULL,
              help = paste("Comma-separated sample-level covariates for trajectory DE",
                           "(log10 library size is always included)")),
  make_option("--min_samples", type = "integer", default = 6,
              help = "Minimum pseudobulk samples for a cell type to be tested [default: %default]"),
  make_option("--fdr", type = "double", default = 0.05,
              help = "Adjusted p-value cutoff for trajectory genes [default: %default]"),
  make_option("--top_genes", type = "integer", default = 20,
              help = "Top up- and down-regulated genes per heatmap [default: %default]"),

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
              help = paste("Run on the ifnb demo dataset with donor IDs",
                           "(16 donor x condition samples; downloads ~5 MB once)"))
)

parser <- OptionParser(
  option_list = option_list,
  description = "\nSample-level analysis with scSLIDE (Seurat >= 6): landmark density per sample, diffusion map trajectories, and trajectory DE.",
  epilogue = paste(
    "Examples:",
    "  Rscript 07_sample_level.R --demo --output ./demo_results",
    "  Rscript 07_sample_level.R --input integrated.rds --sample_col donor \\",
    "    --response disease --celltype_col cell_type --reduction integrated \\",
    "    --latent_vars age,sex --output ./results",
    sep = "\n"
  )
)
args <- parse_args(parser)

# -----------------------------------------------------------------------------
# Validate Arguments
# -----------------------------------------------------------------------------

if (args$demo) {
  # ifnb with donor IDs: one sample per donor and condition
  args$sample_col <- "sample"
  args$response <- "stim"
  args$celltype_col <- "seurat_annotations"
  args$reduction <- "integrated"
}

if (!args$demo && is.null(args$input)) {
  print_help(parser)
  stop("Error: --input is required (or use --demo for demo dataset)")
}
if (is.null(args$sample_col) || is.null(args$response)) {
  stop("Error: --sample_col and --response are required")
}
if (!args$modality %in% c("multi", "single")) {
  stop(sprintf("Unknown --modality '%s'. Choose multi or single", args$modality))
}
if (!args$pls_method %in% c("cppls", "plsr", "spls")) {
  stop(sprintf("Unknown --pls_method '%s'. Choose cppls, plsr or spls", args$pls_method))
}
if (args$pls_method == "spls" && !requireNamespace("spls", quietly = TRUE)) {
  stop("--pls_method spls needs the spls package: install.packages('spls')")
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

dims_use <- parse_dims(args$dims)
traj_dcs <- as.integer(split_list(args$traj_dcs))
if (length(traj_dcs) == 0 || anyNA(traj_dcs) || any(traj_dcs < 1)) {
  stop("--traj_dcs must be a comma-separated list of component numbers, e.g. 1,2")
}
latent_vars <- split_list(args$latent_vars)

set.seed(args$seed)

# -----------------------------------------------------------------------------
# Setup
# -----------------------------------------------------------------------------

if (!dir.exists(args$output)) {
  dir.create(args$output, recursive = TRUE)
}

log_file <- file.path(args$output, "analysis.log")
log_env <- init_logging(log_file = log_file, verbose = args$verbose)

log_message("Starting sample-level (scSLIDE) analysis", log_env)
log_message(sprintf("Output directory: %s", args$output), log_env)

log_message("Loading required packages...", log_env)
load_seurat_deps(quiet = !args$verbose)
setup_threads(args$threads, log_env)

if (!requireNamespace("destiny", quietly = TRUE)) {
  stop(paste("The diffusion map step needs the destiny package (Bioconductor).",
             "Install it with: Rscript setup.R --all"), call. = FALSE)
}
has_complexheatmap <- requireNamespace("ComplexHeatmap", quietly = TRUE)
if (!has_complexheatmap) {
  log_message("ComplexHeatmap not installed: trajectory heatmaps will be skipped", log_env)
}

# -----------------------------------------------------------------------------
# Load Data
# -----------------------------------------------------------------------------

if (args$demo) {
  log_message("Loading ifnb demo dataset with donor IDs...", log_env)
  seurat_obj <- add_ifnb_donor_ids(load_demo_dataset("ifnb"))
  seurat_obj$sample <- paste0(seurat_obj$stim, "-", seurat_obj$donor_id)

  # Cell-type embedding: CCA integration across conditions, so the embedding
  # reflects cell type rather than stimulation (as in 03_integration.R)
  log_message("Preparing demo: normalization, PCA and CCA integration...", log_env)
  seurat_obj[["RNA"]] <- as(seurat_obj[["RNA"]], "Assay5")
  seurat_obj[["RNA"]] <- split(seurat_obj[["RNA"]], f = seurat_obj$stim)
  seurat_obj <- NormalizeData(seurat_obj, verbose = FALSE)
  seurat_obj <- FindVariableFeatures(seurat_obj, verbose = FALSE)
  seurat_obj <- ScaleData(seurat_obj, verbose = FALSE)
  seurat_obj <- RunPCA(seurat_obj, verbose = FALSE)
  seurat_obj <- IntegrateLayers(seurat_obj, method = CCAIntegration, orig.reduction = "pca",
                                new.reduction = "integrated", verbose = FALSE)
  seurat_obj[["RNA"]] <- JoinLayers(seurat_obj[["RNA"]])
} else {
  log_message(sprintf("Loading Seurat object: %s", args$input), log_env)
  seurat_obj <- readRDS(args$input)
}

# -----------------------------------------------------------------------------
# Check Inputs
# -----------------------------------------------------------------------------

meta <- seurat_obj[[]]
for (col in c(args$sample_col, args$response, args$celltype_col, latent_vars)) {
  if (!col %in% colnames(meta)) {
    stop(sprintf("Column '%s' not found in metadata. Available: %s",
                 col, paste(colnames(meta), collapse = ", ")))
  }
}
if (!args$assay %in% Assays(seurat_obj)) {
  stop(sprintf("Assay '%s' not found. Available: %s",
               args$assay, paste(Assays(seurat_obj), collapse = ", ")))
}
DefaultAssay(seurat_obj) <- args$assay

if (anyNA(meta[[args$response]])) {
  stop(sprintf("--response column '%s' has missing values; remove those cells first", args$response))
}
# The response and covariates describe samples, so each must have one value
# per sample (GenerateSampleObject only carries such columns over)
for (col in c(args$response, latent_vars)) {
  values_per_sample <- tapply(as.character(meta[[col]]), meta[[args$sample_col]],
                              function(x) length(unique(x)))
  if (any(values_per_sample > 1)) {
    option <- if (col == args$response) "--response" else "--latent_vars"
    stop(sprintf(paste("%s '%s' must be a sample-level variable (one value per sample),",
                       "but it varies within: %s"),
                 option, col, paste(head(names(values_per_sample)[values_per_sample > 1], 5),
                                    collapse = ", ")))
  }
}

n_samples <- length(unique(meta[[args$sample_col]]))
log_message(sprintf("Cells: %d | Samples: %d", ncol(seurat_obj), n_samples), log_env)
if (n_samples < 4) {
  stop(sprintf("Sample-level analysis needs at least 4 samples in '%s' (found %d)",
               args$sample_col, n_samples))
}
if (n_samples < 10) {
  log_message(sprintf(paste("NOTE: only %d samples. Diffusion components beyond the main",
                            "difference between samples may be uninformative."), n_samples), log_env)
}

# Sketching needs a v5 assay with one counts layer and normalized data
if (inherits(seurat_obj[[args$assay]], "Assay")) {
  log_message("Converting assay to the Seurat v5 format (needed for sketching)", log_env)
  seurat_obj[[args$assay]] <- as(seurat_obj[[args$assay]], "Assay5")
}
if (length(Layers(seurat_obj[[args$assay]], search = "counts")) > 1) {
  seurat_obj[[args$assay]] <- JoinLayers(seurat_obj[[args$assay]])
}
if (length(Layers(seurat_obj[[args$assay]], search = "data")) == 0) {
  log_message("No normalized data found: running NormalizeData", log_env)
  seurat_obj <- NormalizeData(seurat_obj, verbose = FALSE)
}
if (length(VariableFeatures(seurat_obj)) == 0) {
  seurat_obj <- FindVariableFeatures(seurat_obj, nfeatures = 2000, verbose = FALSE)
}

if (args$modality == "multi") {
  if (!args$reduction %in% Reductions(seurat_obj)) {
    if (args$reduction != "pca") {
      stop(sprintf("Reduction '%s' not found. Available: %s",
                   args$reduction, paste(Reductions(seurat_obj), collapse = ", ")))
    }
    log_message("No PCA found: running ScaleData and RunPCA", log_env)
    seurat_obj <- ScaleData(seurat_obj, verbose = FALSE)
    seurat_obj <- RunPCA(seurat_obj, verbose = FALSE)
  }
  n_red <- ncol(Embeddings(seurat_obj, args$reduction))
  if (max(dims_use) > n_red) {
    stop(sprintf("--dims goes up to %d but '%s' has %d dimensions",
                 max(dims_use), args$reduction, n_red))
  }
}

# Landmarks are sketched from the training set, so there can't be more
cells_per_sample <- table(meta[[args$sample_col]])
n_training <- sum(pmin(cells_per_sample, args$ncells_per_sample))
ncells_landmark <- min(args$ncells_landmark, n_training)
if (ncells_landmark < args$ncells_landmark) {
  log_message(sprintf("Using %d landmark cells (the training set has %d cells)",
                      ncells_landmark, n_training), log_env)
}

# -----------------------------------------------------------------------------
# Landmarks and Cell-to-Landmark Neighbors
# -----------------------------------------------------------------------------

log_message("Preparing sample object: sketching, PLS, and landmark neighbors...", log_env)

nn_name <- if (args$modality == "multi") "weighted.nn" else "single.nn"
seurat_obj <- PrepareSampleObject(
  object = seurat_obj,
  assay = args$assay,
  add.hvg = !args$no_cor_genes,
  group.by.CorTest = args$celltype_col,
  Y = args$response,
  sketch.training = TRUE,
  group.by.Sketch = args$sample_col,
  ncells.per.group = args$ncells_per_sample,
  ncells.landmark = ncells_landmark,
  ncomp = args$ncomp,
  pls.function = args$pls_method,
  modality = args$modality,
  k.nn = args$k_nn,
  name.reduction.1 = args$reduction,
  dims.reduction.1 = dims_use,
  nn.name = nn_name,
  fix.wnn.weights = c(0.5, 0.5),
  verbose = args$verbose
)

# -----------------------------------------------------------------------------
# Sample-Level Object
# -----------------------------------------------------------------------------

log_message("Building the landmark x sample object...", log_env)

sample_obj <- GenerateSampleObject(
  object = seurat_obj,
  sketch.assay = "LANDMARK",
  nn.name = nn_name,
  group.by = args$sample_col,
  rename.group.by = args$celltype_col,
  k.nn = args$k_nn,
  normalization.method = "ChiSquared",
  add.meta.data = TRUE,
  return.seurat = TRUE,
  verbose = args$verbose
)
sample_obj <- ScaleData(sample_obj, features = rownames(sample_obj),
                        do.scale = FALSE, do.center = TRUE, verbose = FALSE)
log_message(sprintf("Sample object: %d landmarks x %d samples",
                    nrow(sample_obj), ncol(sample_obj)), log_env)

# -----------------------------------------------------------------------------
# Diffusion Map
# -----------------------------------------------------------------------------

# destiny needs fewer components than samples
n_dcs <- min(args$n_dcs, ncol(sample_obj) - 2)
log_message(sprintf("Running diffusion map (%d components)...", n_dcs), log_env)
sample_obj <- RunDiffusionMap(sample_obj, assay = "LMC", layer = "scale.data",
                              features = rownames(sample_obj), ncomp = n_dcs,
                              seed.use = args$seed, verbose = args$verbose)

# Components are weighted by their eigenvalues; ones with non-positive
# eigenvalues come out as all zeros and carry no information
dc_emb <- Embeddings(sample_obj, "DiffMap")
informative <- which(apply(dc_emb, 2, stats::sd) > 1e-8)
informative <- informative[informative == seq_along(informative)]
if (length(informative) == 0) {
  stop("The diffusion map found no varying components; check the sample design")
}
log_message(sprintf("Informative diffusion components: %d", length(informative)), log_env)

for (k in informative) {
  sample_obj[[paste0("traj_DC", k)]] <- dc_emb[, k]
}

skipped_dcs <- setdiff(traj_dcs, informative)
if (length(skipped_dcs) > 0) {
  log_message(sprintf("Skipping uninformative components: DC%s",
                      paste(skipped_dcs, collapse = ", DC")), log_env)
}
traj_dcs <- intersect(traj_dcs, informative)

# Per-sample table: metadata plus the informative components (without the
# leverage scores left over from sketching)
sample_meta <- sample_obj[[]]
sample_table <- data.frame(sample_id = colnames(sample_obj),
                           sample_meta[, !grepl("^leverage\\.score", colnames(sample_meta)),
                                       drop = FALSE],
                           check.names = FALSE)
write.csv(sample_table, file.path(args$output, "sample_embedding.csv"), row.names = FALSE)

# -----------------------------------------------------------------------------
# Sample-Level Plots
# -----------------------------------------------------------------------------

log_message("Generating sample-level plots...", log_env)

response_is_numeric <- is.numeric(sample_meta[[args$response]])
dim_pairs <- list(c(1, 2), c(1, 3))
dim_pairs <- dim_pairs[vapply(dim_pairs, function(d) all(d %in% informative), logical(1))]

if (length(dim_pairs) > 0) {
  p_dm <- wrap_plots(lapply(dim_pairs, function(d) {
    DimPlot(sample_obj, reduction = "DiffMap", group.by = args$response, dims = d,
            pt.size = 3)
  }), nrow = 1)
  save_plot(p_dm, "01_diffusion_map.png", output_dir = args$output,
            width = 6 * length(dim_pairs), height = 5)
}

# Each informative component (up to 3) against the response
dc_long <- do.call(rbind, lapply(head(informative, 3), function(k) {
  data.frame(component = paste0("DC", k), value = dc_emb[, k],
             response = sample_meta[[args$response]])
}))
p_resp <- ggplot(dc_long, aes(x = response, y = value)) +
  facet_wrap(~component, scales = "free_y") +
  theme_minimal() +
  labs(x = args$response, y = "Diffusion component value",
       title = sprintf("Sample positions by %s", args$response))
p_resp <- if (response_is_numeric) {
  p_resp + geom_point(size = 2) + geom_smooth(method = "lm", formula = y ~ x, se = FALSE)
} else {
  p_resp + geom_boxplot(outlier.shape = NA) + geom_jitter(width = 0.15, height = 0, size = 2)
}
save_plot(p_resp, "02_components_by_response.png", output_dir = args$output,
          width = 4 * min(3, length(informative)) + 1, height = 5)

# -----------------------------------------------------------------------------
# Pseudobulk Trajectory DE
# -----------------------------------------------------------------------------

log_message("Pseudobulking by sample and cell type...", log_env)

# Group by copies with "_" replaced: AggregateExpression joins group values
# with "_" and replaces "_" inside them with "-"
seurat_obj$pb_sample <- gsub("_", "-", as.character(seurat_obj[[args$sample_col]][, 1]))
seurat_obj$pb_celltype <- if (is.null(args$celltype_col)) {
  rep("all cells", ncol(seurat_obj))
} else {
  gsub("_", "-", as.character(seurat_obj[[args$celltype_col]][, 1]))
}
bulk <- AggregateExpression(seurat_obj, assays = args$assay,
                            group.by = c("pb_sample", "pb_celltype"),
                            return.seurat = TRUE, verbose = FALSE)
# Pseudobulk names are "<sample>_<cell type>" (neither contains "_"), or just
# "<sample>" when there is a single cell type. Read the groups back from the
# names, since a constant group isn't kept as metadata.
name_parts <- strsplit(colnames(bulk), "_", fixed = TRUE)
only_celltype <- unique(seurat_obj$pb_celltype)
bulk$pb_sample <- vapply(name_parts, `[`, character(1), 1)
bulk$pb_celltype <- vapply(name_parts, function(p) {
  if (length(p) > 1) p[2] else only_celltype[1]
}, character(1))

# Sample-level values (components and covariates) for each pseudobulk sample
sample_key <- gsub("_", "-", as.character(sample_meta[[args$sample_col]]))
row_for <- match(bulk$pb_sample, sample_key)
for (col in c(args$response, latent_vars, paste0("traj_DC", traj_dcs))) {
  bulk[[col]] <- sample_meta[[col]][row_for]
}
bulk$log_count <- log10(Matrix::colSums(LayerData(bulk, assay = args$assay, layer = "counts")))

celltypes <- sort(unique(bulk$pb_celltype))
if (args$de_celltypes != "all") {
  requested <- gsub("_", "-", split_list(args$de_celltypes))
  unknown <- setdiff(requested, celltypes)
  if (length(unknown) > 0) {
    stop(sprintf("--de_celltypes not found: %s. Available: %s",
                 paste(unknown, collapse = ", "), paste(celltypes, collapse = ", ")))
  }
  celltypes <- requested
}

safe_name <- function(x) gsub("[^A-Za-z0-9.-]+", "_", x)
de_summary <- list()
heatmap_index <- 0

for (ct in celltypes) {
  ct_bulk <- bulk[, bulk$pb_celltype == ct]
  if (ncol(ct_bulk) < args$min_samples) {
    log_message(sprintf("  %s: %d samples (< --min_samples %d), skipped",
                        ct, ncol(ct_bulk), args$min_samples), log_env)
    de_summary[[length(de_summary) + 1]] <- data.frame(
      celltype = ct, component = NA, n_samples = ncol(ct_bulk), n_tested = 0,
      n_significant = 0, status = "too few samples")
    next
  }

  for (k in traj_dcs) {
    traj_var <- paste0("traj_DC", k)
    log_message(sprintf("  %s: trajectory DE along DC%d (%d samples)", ct, k, ncol(ct_bulk)),
                log_env)
    de <- tryCatch(
      TrajDETest(ct_bulk, assay = args$assay, layer = "counts", traj.var = traj_var,
                 latent.vars = c(latent_vars, "log_count"), verbose = FALSE),
      error = function(e) {
        log_message(sprintf("    failed: %s", conditionMessage(e)), log_env)
        NULL
      })
    if (is.null(de)) {
      de_summary[[length(de_summary) + 1]] <- data.frame(
        celltype = ct, component = paste0("DC", k), n_samples = ncol(ct_bulk), n_tested = 0,
        n_significant = 0, status = "failed")
      next
    }

    de <- de[order(de$p_Traj), ]
    write.csv(de, file.path(args$output, sprintf("traj_de_%s_DC%d.csv", safe_name(ct), k)),
              row.names = FALSE)
    sig <- de[!is.na(de$adj_pval) & de$adj_pval < args$fdr, ]
    de_summary[[length(de_summary) + 1]] <- data.frame(
      celltype = ct, component = paste0("DC", k), n_samples = ncol(ct_bulk),
      n_tested = nrow(de), n_significant = nrow(sig), status = "tested")
    log_message(sprintf("    %d genes tested, %d with adjusted p < %.2f",
                        nrow(de), nrow(sig), args$fdr), log_env)

    # Heatmap of the top genes going up and down along the component
    if (has_complexheatmap && nrow(sig) > 0) {
      up <- head(sig$gene_ID[sig$beta_Traj > 0], args$top_genes)
      down <- head(sig$gene_ID[sig$beta_Traj < 0], args$top_genes)
      hm_genes <- c(up, rev(down))
      hm_bulk <- NormalizeData(ct_bulk, verbose = FALSE)
      hm_bulk <- ScaleData(hm_bulk, features = hm_genes, verbose = FALSE)
      ht <- PlotAnnotatedHeatmap(hm_bulk, features = hm_genes, sort.by = traj_var,
                                 anno.var = c(args$response, traj_var),
                                 row.font.size = 8)
      heatmap_index <- heatmap_index + 1
      save_complex_heatmap(
        ht, sprintf("%02d_traj_heatmap_%s_DC%d.png", heatmap_index + 2, safe_name(ct), k),
        output_dir = args$output, width = max(7, ncol(hm_bulk) * 0.3 + 4),
        height = max(5, length(hm_genes) * 0.18 + 3))
    }
  }
}

de_summary <- do.call(rbind, de_summary)
write.csv(de_summary, file.path(args$output, "traj_de_summary.csv"), row.names = FALSE)

# -----------------------------------------------------------------------------
# Save Results
# -----------------------------------------------------------------------------

log_message("Saving results...", log_env)

save_seurat(sample_obj, file.path(args$output, "sample_object.rds"))

writeLines(capture.output(sessionInfo()), file.path(args$output, "session_info.txt"))
params <- as.data.frame(t(as.data.frame(args)))
write.csv(params, file.path(args$output, "parameters.csv"))

log_message("Sample-level analysis complete!", log_env)
log_message(sprintf("Results saved to: %s", args$output), log_env)
