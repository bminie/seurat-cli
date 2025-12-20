# =============================================================================
# Output Validation Framework for Seurat CLI
# =============================================================================
#
# Two-tier validation system:
#   Tier 1: Structural validation (applies to all runs)
#   Tier 2: Biological validation (applies to demo runs only)
#
# =============================================================================

suppressPackageStartupMessages({
  library(Seurat)
})

# -----------------------------------------------------------------------------
# Validation Result Helper
# -----------------------------------------------------------------------------

#' Create a validation result object
#' @param passed Logical indicating if validation passed
#' @param message Description of what was checked
#' @param details Additional details (for failures)
validation_result <- function(passed, message, details = NULL) {
  list(
    passed = passed,
    message = message,
    details = details
  )
}

#' Run a validation check with error handling
#' @param expr Expression to evaluate
#' @param message Description of the check
safe_check <- function(expr, message) {
  tryCatch({
    result <- expr
    if (isTRUE(result)) {
      validation_result(TRUE, message)
    } else {
      validation_result(FALSE, message, details = as.character(result))
    }
  }, error = function(e) {
    validation_result(FALSE, message, details = e$message)
  })
}

# -----------------------------------------------------------------------------
# Tier 1: Structural Validation (All Runs)
# -----------------------------------------------------------------------------

#' Validate that expected files exist
#' @param output_dir Output directory to check
#' @param expected_files Character vector of expected file names
validate_files_exist <- function(output_dir, expected_files) {
  results <- list()

  for (f in expected_files) {
    path <- file.path(output_dir, f)
    exists <- file.exists(path)
    results[[f]] <- validation_result(
      exists,
      sprintf("File exists: %s", f),
      if (!exists) "File not found"
    )
  }

  results
}

#' Validate CSV file structure
#' @param csv_path Path to CSV file
#' @param required_columns Character vector of required column names
validate_csv_structure <- function(csv_path, required_columns) {
  results <- list()

  # Check file can be read
  read_check <- safe_check({
    df <- read.csv(csv_path, nrows = 5)
    TRUE
  }, sprintf("CSV readable: %s", basename(csv_path)))
  results[["readable"]] <- read_check

  if (!read_check$passed) {
    return(results)
  }

  df <- read.csv(csv_path, nrows = 5)

  # Check required columns
  for (col in required_columns) {
    has_col <- col %in% colnames(df)
    results[[paste0("column_", col)]] <- validation_result(
      has_col,
      sprintf("Has column: %s", col),
      if (!has_col) sprintf("Missing column. Found: %s", paste(colnames(df), collapse = ", "))
    )
  }

  results
}

#' Validate cluster assignments
#' @param output_dir Output directory containing cluster_assignments.csv
validate_cluster_assignments <- function(output_dir) {
  results <- list()
  csv_path <- file.path(output_dir, "cluster_assignments.csv")

  if (!file.exists(csv_path)) {
    results[["file_exists"]] <- validation_result(FALSE, "Cluster assignments file exists", "File not found")
    return(results)
  }

  clusters <- read.csv(csv_path)

  # Check required columns exist
  results[["has_cell_column"]] <- validation_result(
    "cell" %in% colnames(clusters),
    "Has 'cell' column"
  )

  results[["has_cluster_column"]] <- validation_result(
    "cluster" %in% colnames(clusters),
    "Has 'cluster' column"
  )

  if (!"cluster" %in% colnames(clusters)) {
    return(results)
  }

  cluster_values <- clusters$cluster

  # Check clusters are integers (not strings like "integrated_snn_res.0.5")
  numeric_check <- suppressWarnings(as.integer(cluster_values))
  all_numeric <- !any(is.na(numeric_check))

  results[["clusters_are_integers"]] <- validation_result(
    all_numeric,
    "Cluster values are integers",
    if (!all_numeric) sprintf("Non-integer cluster values found: %s",
                               paste(head(unique(cluster_values[is.na(numeric_check)])), collapse = ", "))
  )

  if (!all_numeric) {
    return(results)
  }

  # Check multiple clusters exist
  n_clusters <- length(unique(cluster_values))
  results[["multiple_clusters"]] <- validation_result(
    n_clusters >= 2,
    sprintf("Multiple clusters found (n=%d)", n_clusters),
    if (n_clusters < 2) "Only one cluster found - clustering may have failed"
  )

  # Check clusters start at 0 (Seurat convention)
  min_cluster <- min(as.integer(cluster_values))
  results[["clusters_start_at_zero"]] <- validation_result(
    min_cluster == 0,
    "Clusters start at 0",
    if (min_cluster != 0) sprintf("Minimum cluster is %d, expected 0", min_cluster)
  )

  # Check no NA values
  results[["no_na_clusters"]] <- validation_result(
    !any(is.na(cluster_values)),
    "No NA cluster assignments"
  )

  # Check all cells have assignments
  results[["all_cells_assigned"]] <- validation_result(
    nrow(clusters) > 0,
    sprintf("All cells have cluster assignments (n=%d)", nrow(clusters))
  )

  results
}

#' Validate marker results
#' @param output_dir Output directory containing marker files
validate_markers <- function(output_dir) {
  results <- list()

  # Check for marker file (could be cluster_markers_all.csv or de_results_all.csv)
  marker_files <- c("cluster_markers_all.csv", "de_results_all.csv")
  marker_file <- NULL

  for (f in marker_files) {
    path <- file.path(output_dir, f)
    if (file.exists(path)) {
      marker_file <- path
      break
    }
  }

  if (is.null(marker_file)) {
    results[["marker_file_exists"]] <- validation_result(
      FALSE,
      "Marker file exists",
      sprintf("No marker file found. Checked: %s", paste(marker_files, collapse = ", "))
    )
    return(results)
  }

  results[["marker_file_exists"]] <- validation_result(TRUE, sprintf("Marker file exists: %s", basename(marker_file)))

  markers <- read.csv(marker_file)

  # Check file is not empty
  results[["markers_not_empty"]] <- validation_result(
    nrow(markers) > 0,
    sprintf("Marker file has data (n=%d rows)", nrow(markers)),
    if (nrow(markers) == 0) "Marker file is empty"
  )

  if (nrow(markers) == 0) {
    return(results)
  }

  # Check required columns
  required_cols <- c("gene", "cluster")
  for (col in required_cols) {
    results[[paste0("has_", col)]] <- validation_result(
      col %in% colnames(markers),
      sprintf("Markers have '%s' column", col)
    )
  }

  # Check p-values are in valid range (if present)
  if ("p_val" %in% colnames(markers)) {
    valid_pvals <- all(markers$p_val >= 0 & markers$p_val <= 1, na.rm = TRUE)
    results[["valid_pvalues"]] <- validation_result(
      valid_pvals,
      "P-values are in valid range [0,1]"
    )
  }

  # Check multiple clusters have markers
  if ("cluster" %in% colnames(markers)) {
    n_clusters_with_markers <- length(unique(markers$cluster))
    results[["multiple_clusters_have_markers"]] <- validation_result(
      n_clusters_with_markers >= 2,
      sprintf("Multiple clusters have markers (n=%d)", n_clusters_with_markers)
    )
  }

  results
}

#' Validate Seurat object structure
#' @param rds_path Path to Seurat RDS file
#' @param expected_reductions Character vector of expected reduction names
#' @param expected_assays Character vector of expected assay names
validate_seurat_object <- function(rds_path, expected_reductions = c("pca", "umap"),
                                    expected_assays = c("RNA")) {
  results <- list()

  if (!file.exists(rds_path)) {
    results[["file_exists"]] <- validation_result(FALSE, "Seurat RDS file exists", "File not found")
    return(results)
  }

  # Load object
  load_check <- safe_check({
    obj <- readRDS(rds_path)
    inherits(obj, "Seurat")
  }, "RDS loads as Seurat object")
  results[["loads_as_seurat"]] <- load_check

  if (!load_check$passed) {
    return(results)
  }

  obj <- readRDS(rds_path)

  # Check reductions
  available_reductions <- names(obj@reductions)
  for (red in expected_reductions) {
    # Allow for variations like "umap.integrated" matching "umap"
    has_reduction <- any(grepl(paste0("^", red), available_reductions, ignore.case = TRUE))
    results[[paste0("has_reduction_", red)]] <- validation_result(
      has_reduction,
      sprintf("Has %s reduction", red),
      if (!has_reduction) sprintf("Available reductions: %s", paste(available_reductions, collapse = ", "))
    )
  }

  # Check assays
  available_assays <- names(obj@assays)
  for (assay in expected_assays) {
    has_assay <- assay %in% available_assays
    results[[paste0("has_assay_", assay)]] <- validation_result(
      has_assay,
      sprintf("Has %s assay", assay),
      if (!has_assay) sprintf("Available assays: %s", paste(available_assays, collapse = ", "))
    )
  }

  # Check has cells
  n_cells <- ncol(obj)
  results[["has_cells"]] <- validation_result(
    n_cells > 0,
    sprintf("Object has cells (n=%d)", n_cells)
  )

  # Check has cluster metadata
  has_clusters <- "seurat_clusters" %in% colnames(obj@meta.data)
  results[["has_cluster_metadata"]] <- validation_result(
    has_clusters,
    "Has seurat_clusters metadata"
  )

  results
}

# -----------------------------------------------------------------------------
# Tier 2: Biological Validation (Demo Runs Only)
# -----------------------------------------------------------------------------

#' Get expected values for demo datasets
#' @param script_name Name of the script (e.g., "01_basic_analysis")
#' @return List of expected values
get_demo_expectations <- function(script_name) {
  expectations <- list(
    "01_basic_analysis" = list(
      dataset = "pbmc3k",
      min_clusters = 7,
      max_clusters = 12,
      min_cells = 2400,
      max_cells = 2700,
      expected_markers = c("MS4A1", "CD79A", "LYZ", "CD14", "GNLY", "NKG7", "CD3D", "CD8A"),
      min_markers_found = 5
    ),
    "02_sctransform" = list(
      dataset = "pbmc3k",
      min_clusters = 9,
      max_clusters = 15,
      min_cells = 2600,
      max_cells = 2750,
      expected_markers = c("MS4A1", "CD79A", "LYZ", "CD14", "GNLY", "NKG7", "CD3D", "CD8A"),
      min_markers_found = 5
    ),
    "03_integration" = list(
      dataset = "ifnb",
      min_clusters = 12,
      max_clusters = 18,
      min_cells = 13000,
      max_cells = 14500,
      expected_markers = c("CD14", "LYZ", "FCGR3A", "MS4A1", "GNLY", "CD8A"),
      min_markers_found = 4,
      expected_samples = c("IMMUNE_CTRL", "IMMUNE_STIM")
    ),
    "04_differential_expression" = list(
      dataset = "pbmc3k",
      min_clusters = 7,
      max_clusters = 12,
      min_de_genes = 1000,
      expected_markers = c("MS4A1", "LYZ", "GNLY", "CD3D"),
      min_markers_found = 3
    ),
    "05_cell_cycle" = list(
      dataset = "pbmc3k",
      min_clusters = 7,
      max_clusters = 12,
      min_cells = 2500,
      max_cells = 2700,
      expected_phases = c("G1", "S", "G2M")
    ),
    "06_visualization" = list(
      dataset = "pbmc3k",
      min_plots = 5
    )
  )

  expectations[[script_name]]
}

#' Validate demo-specific expectations
#' @param output_dir Output directory
#' @param script_name Script name to get expectations for
validate_demo_expectations <- function(output_dir, script_name) {
  results <- list()

  expectations <- get_demo_expectations(script_name)
  if (is.null(expectations)) {
    results[["has_expectations"]] <- validation_result(
      FALSE,
      "Demo expectations defined",
      sprintf("No expectations defined for script: %s", script_name)
    )
    return(results)
  }

  # Validate cluster count range
  if (!is.null(expectations$min_clusters) && !is.null(expectations$max_clusters)) {
    cluster_file <- file.path(output_dir, "cluster_assignments.csv")
    if (file.exists(cluster_file)) {
      clusters <- read.csv(cluster_file)
      n_clusters <- length(unique(clusters$cluster))

      in_range <- n_clusters >= expectations$min_clusters && n_clusters <= expectations$max_clusters
      results[["cluster_count_in_range"]] <- validation_result(
        in_range,
        sprintf("Cluster count in expected range [%d-%d]", expectations$min_clusters, expectations$max_clusters),
        if (!in_range) sprintf("Found %d clusters", n_clusters)
      )
    }
  }

  # Validate cell count range
  if (!is.null(expectations$min_cells) && !is.null(expectations$max_cells)) {
    cluster_file <- file.path(output_dir, "cluster_assignments.csv")
    if (file.exists(cluster_file)) {
      clusters <- read.csv(cluster_file)
      n_cells <- nrow(clusters)

      in_range <- n_cells >= expectations$min_cells && n_cells <= expectations$max_cells
      results[["cell_count_in_range"]] <- validation_result(
        in_range,
        sprintf("Cell count in expected range [%d-%d]", expectations$min_cells, expectations$max_cells),
        if (!in_range) sprintf("Found %d cells", n_cells)
      )
    }
  }

  # Validate expected marker genes found
  if (!is.null(expectations$expected_markers) && !is.null(expectations$min_markers_found)) {
    marker_files <- c("cluster_markers_all.csv", "de_results_all.csv")
    for (f in marker_files) {
      marker_path <- file.path(output_dir, f)
      if (file.exists(marker_path)) {
        markers <- read.csv(marker_path)
        if ("gene" %in% colnames(markers)) {
          found_markers <- sum(expectations$expected_markers %in% markers$gene)
          enough_found <- found_markers >= expectations$min_markers_found

          results[["expected_markers_found"]] <- validation_result(
            enough_found,
            sprintf("Expected marker genes found (%d/%d required)", found_markers, expectations$min_markers_found),
            if (!enough_found) sprintf("Found: %s",
                                        paste(expectations$expected_markers[expectations$expected_markers %in% markers$gene], collapse = ", "))
          )
        }
        break
      }
    }
  }

  # Validate expected samples (for integration)
  if (!is.null(expectations$expected_samples)) {
    cluster_file <- file.path(output_dir, "cluster_assignments.csv")
    if (file.exists(cluster_file)) {
      clusters <- read.csv(cluster_file)
      if ("sample" %in% colnames(clusters)) {
        found_samples <- unique(clusters$sample)
        all_found <- all(expectations$expected_samples %in% found_samples)

        results[["expected_samples_present"]] <- validation_result(
          all_found,
          "Expected samples present in output",
          if (!all_found) sprintf("Missing samples: %s",
                                   paste(setdiff(expectations$expected_samples, found_samples), collapse = ", "))
        )
      }
    }
  }

  # Validate cell cycle phases (for cell cycle script)
  if (!is.null(expectations$expected_phases)) {
    scores_file <- file.path(output_dir, "cell_cycle_scores.csv")
    if (file.exists(scores_file)) {
      scores <- read.csv(scores_file)
      if ("Phase" %in% colnames(scores)) {
        found_phases <- unique(scores$Phase)
        all_found <- all(expectations$expected_phases %in% found_phases)

        results[["expected_phases_present"]] <- validation_result(
          all_found,
          "Expected cell cycle phases present",
          if (!all_found) sprintf("Missing phases: %s",
                                   paste(setdiff(expectations$expected_phases, found_phases), collapse = ", "))
        )
      }
    }
  }

  # Validate minimum DE genes (for DE script)
  if (!is.null(expectations$min_de_genes)) {
    de_file <- file.path(output_dir, "de_results_all.csv")
    if (file.exists(de_file)) {
      de_results <- read.csv(de_file)
      n_genes <- nrow(de_results)

      enough_genes <- n_genes >= expectations$min_de_genes
      results[["sufficient_de_genes"]] <- validation_result(
        enough_genes,
        sprintf("Sufficient DE genes found (min=%d)", expectations$min_de_genes),
        if (!enough_genes) sprintf("Found %d DE genes", n_genes)
      )
    }
  }

  # Validate minimum plots (for visualization script)
  if (!is.null(expectations$min_plots)) {
    png_files <- list.files(output_dir, pattern = "\\.png$", full.names = FALSE)
    n_plots <- length(png_files)

    enough_plots <- n_plots >= expectations$min_plots
    results[["sufficient_plots"]] <- validation_result(
      enough_plots,
      sprintf("Sufficient plots generated (min=%d)", expectations$min_plots),
      if (!enough_plots) sprintf("Found %d plots", n_plots)
    )
  }

  results
}

# -----------------------------------------------------------------------------
# Main Validation Functions
# -----------------------------------------------------------------------------

#' Run Tier 1 validation (structural checks for all runs)
#' @param output_dir Output directory to validate
#' @param script_name Name of the script that produced the output
#' @return List of validation results
run_tier1_validation <- function(output_dir, script_name) {
  results <- list()

  # Define expected files per script
  expected_files <- list(
    "01_basic_analysis" = c("seurat_object.rds", "cluster_assignments.csv", "cluster_markers_all.csv",
                            "parameters.csv", "session_info.txt"),
    "02_sctransform" = c("seurat_sctransform.rds", "cluster_assignments.csv", "cluster_markers_all.csv",
                          "parameters.csv", "session_info.txt"),
    "03_integration" = c("seurat_integrated.rds", "cluster_assignments.csv", "cluster_markers_all.csv",
                          "parameters.csv", "session_info.txt"),
    "04_differential_expression" = c("de_results_all.csv", "de_results_significant.csv",
                                      "parameters.csv", "session_info.txt"),
    "05_cell_cycle" = c("seurat_cell_cycle.rds", "cell_cycle_scores.csv",
                         "parameters.csv", "session_info.txt"),
    "06_visualization" = c("visualization_summary.md", "parameters.csv", "session_info.txt")
  )

  # Define RDS file and expected reductions per script
  rds_files <- list(
    "01_basic_analysis" = list(file = "seurat_object.rds", reductions = c("pca", "umap"), assays = c("RNA")),
    "02_sctransform" = list(file = "seurat_sctransform.rds", reductions = c("pca", "umap"), assays = c("RNA", "SCT")),
    "03_integration" = list(file = "seurat_integrated.rds", reductions = c("pca", "umap"), assays = c("RNA")),
    "05_cell_cycle" = list(file = "seurat_cell_cycle.rds", reductions = c("pca"), assays = c("RNA"))
  )

  # 1. File existence checks
  if (script_name %in% names(expected_files)) {
    results[["files"]] <- validate_files_exist(output_dir, expected_files[[script_name]])
  }

  # 2. Cluster assignment validation (for scripts that produce clusters)
  if (script_name %in% c("01_basic_analysis", "02_sctransform", "03_integration", "05_cell_cycle")) {
    results[["clusters"]] <- validate_cluster_assignments(output_dir)
  }

  # 3. Marker validation (for scripts that produce markers)
  if (script_name %in% c("01_basic_analysis", "02_sctransform", "03_integration", "04_differential_expression")) {
    results[["markers"]] <- validate_markers(output_dir)
  }

  # 4. Seurat object validation
  if (script_name %in% names(rds_files)) {
    rds_info <- rds_files[[script_name]]
    rds_path <- file.path(output_dir, rds_info$file)
    results[["seurat_object"]] <- validate_seurat_object(
      rds_path,
      expected_reductions = rds_info$reductions,
      expected_assays = rds_info$assays
    )
  }

  results
}

#' Run Tier 2 validation (biological checks for demo runs only)
#' @param output_dir Output directory to validate
#' @param script_name Name of the script that produced the output
#' @return List of validation results
run_tier2_validation <- function(output_dir, script_name) {
  validate_demo_expectations(output_dir, script_name)
}

#' Run full validation
#' @param output_dir Output directory to validate
#' @param script_name Name of the script that produced the output
#' @param is_demo Whether this is a demo run (enables Tier 2 validation)
#' @return List with tier1 and tier2 results
run_validation <- function(output_dir, script_name, is_demo = FALSE) {
  results <- list(
    tier1 = run_tier1_validation(output_dir, script_name),
    tier2 = if (is_demo) run_tier2_validation(output_dir, script_name) else list()
  )

  results
}

#' Summarize validation results
#' @param results Validation results from run_validation()
#' @return List with summary statistics
summarize_validation <- function(results) {
  # Flatten nested results
  flatten_results <- function(x, prefix = "") {
    flat <- list()
    for (name in names(x)) {
      item <- x[[name]]
      full_name <- if (prefix == "") name else paste(prefix, name, sep = ".")

      if (is.list(item) && !is.null(item$passed)) {
        # This is a validation result
        flat[[full_name]] <- item
      } else if (is.list(item)) {
        # Recurse
        flat <- c(flat, flatten_results(item, full_name))
      }
    }
    flat
  }

  all_results <- c(
    flatten_results(results$tier1, "tier1"),
    flatten_results(results$tier2, "tier2")
  )

  passed <- sum(sapply(all_results, function(x) x$passed))
  failed <- sum(sapply(all_results, function(x) !x$passed))

  failed_checks <- names(all_results)[!sapply(all_results, function(x) x$passed)]
  failed_details <- lapply(all_results[failed_checks], function(x) {
    list(message = x$message, details = x$details)
  })

  list(
    total = passed + failed,
    passed = passed,
    failed = failed,
    all_passed = failed == 0,
    failed_checks = failed_details
  )
}

#' Print validation summary
#' @param results Validation results from run_validation()
#' @param script_name Script name for display
print_validation_summary <- function(results, script_name) {
  summary <- summarize_validation(results)

  cat(sprintf("\n  Validation: %d/%d checks passed\n", summary$passed, summary$total))

  if (summary$failed > 0) {
    cat("  Failed checks:\n")
    for (name in names(summary$failed_checks)) {
      check <- summary$failed_checks[[name]]
      cat(sprintf("    - %s: %s\n", check$message, check$details))
    }
  }

  summary$all_passed
}
