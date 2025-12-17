#!/usr/bin/env Rscript

# =============================================================================
# Setup Script: Install all required packages for Seurat CLI
# Usage: Rscript setup.R [--all] [--bioc] [--minimal]
# =============================================================================

cat("\n")
cat("=", rep("=", 60), "\n", sep = "")
cat("  SEURAT CLI - DEPENDENCY INSTALLER\n")
cat("=", rep("=", 60), "\n\n", sep = "")

# -----------------------------------------------------------------------------
# Parse Arguments
# -----------------------------------------------------------------------------

args <- commandArgs(trailingOnly = TRUE)
install_all <- "--all" %in% args
install_bioc <- "--bioc" %in% args || install_all
minimal <- "--minimal" %in% args

if (length(args) == 0) {
  cat("Usage: Rscript setup.R [options]\n\n")
  cat("Options:\n")
  cat("  --minimal  Install only core packages (Seurat + basics)\n")
  cat("  --bioc     Also install Bioconductor packages (MAST, DESeq2)\n")
  cat("  --all      Install everything including optional packages\n")
  cat("\nNo options specified - installing recommended packages...\n\n")
}

# Increase timeout for large downloads (datasets can be 400MB+)
options(timeout = 600)

# -----------------------------------------------------------------------------
# Helper Functions
# -----------------------------------------------------------------------------

install_if_missing <- function(pkg, source = "CRAN", repo = NULL) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    cat(sprintf("  Installing %s from %s...\n", pkg, source))

    tryCatch({
      if (source == "CRAN") {
        install.packages(pkg, quiet = TRUE)
      } else if (source == "Bioconductor") {
        if (!requireNamespace("BiocManager", quietly = TRUE)) {
          install.packages("BiocManager", quiet = TRUE)
        }
        BiocManager::install(pkg, quiet = TRUE, update = FALSE)
      } else if (source == "GitHub") {
        if (!requireNamespace("remotes", quietly = TRUE)) {
          install.packages("remotes", quiet = TRUE)
        }
        remotes::install_github(repo, quiet = TRUE, upgrade = "never")
      }
    }, error = function(e) {
      cat(sprintf("    Error: %s\n", e$message))
    })

    if (requireNamespace(pkg, quietly = TRUE)) {
      cat(sprintf("    ✓ %s installed successfully\n", pkg))
      return(TRUE)
    } else {
      cat(sprintf("    ✗ %s installation failed\n", pkg))
      return(FALSE)
    }
  } else {
    cat(sprintf("  ✓ %s already installed\n", pkg))
    return(TRUE)
  }
}

install_seurat_dataset <- function(dataset_name) {
  # Check if dataset package is installed
  pkg_name <- paste0(dataset_name, ".SeuratData")
  if (!requireNamespace(pkg_name, quietly = TRUE)) {
    cat(sprintf("  Installing %s dataset...\n", dataset_name))
    tryCatch({
      SeuratData::InstallData(dataset_name)
      if (requireNamespace(pkg_name, quietly = TRUE)) {
        cat(sprintf("    ✓ %s installed\n", dataset_name))
        return(TRUE)
      } else {
        cat(sprintf("    ✗ %s installation failed\n", dataset_name))
        return(FALSE)
      }
    }, error = function(e) {
      cat(sprintf("    ✗ %s installation failed: %s\n", dataset_name, e$message))
      return(FALSE)
    })
  } else {
    cat(sprintf("  ✓ %s already installed\n", dataset_name))
    return(TRUE)
  }
}

# -----------------------------------------------------------------------------
# Results Tracking
# -----------------------------------------------------------------------------

results <- list()

# =============================================================================
# STEP 1: Core Packages (always installed)
# =============================================================================

cat("Step 1: Installing core packages...\n")
core_packages <- c(
  "Seurat",
  "dplyr",
  "ggplot2",
  "patchwork",
  "optparse",
  "Matrix",
  "future",
  "viridis"
)

for (pkg in core_packages) {
  results[[pkg]] <- install_if_missing(pkg, "CRAN")
}

# =============================================================================
# STEP 2: Recommended Packages (unless --minimal)
# =============================================================================

if (!minimal) {
  cat("\nStep 2: Installing recommended packages...\n")

  # GitHub packages
  cat("\n  From GitHub:\n")
  results[["SeuratData"]] <- install_if_missing(
    "SeuratData", "GitHub", "satijalab/seurat-data"
  )
  results[["presto"]] <- install_if_missing(
    "presto", "GitHub", "immunogenomics/presto"
  )

  # CRAN packages
  cat("\n  From CRAN:\n")
  cran_recommended <- c("hdf5r", "sctransform")
  for (pkg in cran_recommended) {
    results[[pkg]] <- install_if_missing(pkg, "CRAN")
  }

  # Bioconductor packages (for SCTransform v2)
  cat("\n  From Bioconductor:\n")
  results[["glmGamPoi"]] <- install_if_missing("glmGamPoi", "Bioconductor")
}

# =============================================================================
# STEP 3: Bioconductor DE Packages (if --bioc or --all)
# =============================================================================

if (install_bioc) {
  cat("\nStep 3: Installing Bioconductor DE packages...\n")
  bioc_packages <- c(
    "MAST",
    "DESeq2",
    "SingleCellExperiment"
  )

  for (pkg in bioc_packages) {
    results[[pkg]] <- install_if_missing(pkg, "Bioconductor")
  }
}

# =============================================================================
# STEP 4: Optional Integration Packages (if --all)
# =============================================================================

if (install_all) {
  cat("\nStep 4: Installing optional integration packages...\n")

  # CRAN
  cat("\n  From CRAN:\n")
  results[["harmony"]] <- install_if_missing("harmony", "CRAN")

  # Bioconductor
  cat("\n  From Bioconductor:\n")
  bioc_optional <- c("batchelor", "scran", "scater")
  for (pkg in bioc_optional) {
    results[[pkg]] <- install_if_missing(pkg, "Bioconductor")
  }
}

# =============================================================================
# STEP 5: Demo Datasets (if SeuratData is available)
# =============================================================================

if (!minimal) {
  cat("\nStep 5: Installing demo datasets...\n")
  cat("  (This may take a few minutes for large datasets)\n\n")

  if (requireNamespace("SeuratData", quietly = TRUE)) {
    # pbmc3k - used by scripts 01, 02, 04, 05, 06
    install_seurat_dataset("pbmc3k")

    # ifnb - used by script 03 (integration)
    install_seurat_dataset("ifnb")
  } else {
    cat("  ⚠ Skipping datasets (SeuratData not installed)\n")
    cat("    Run setup.R again after installing SeuratData\n")
  }
}

# =============================================================================
# Summary
# =============================================================================

cat("\n")
cat("=", rep("=", 60), "\n", sep = "")
cat("  INSTALLATION SUMMARY\n")
cat("=", rep("=", 60), "\n\n", sep = "")

installed <- sum(unlist(results))
total <- length(results)

cat(sprintf("Packages: %d/%d installed successfully\n", installed, total))

if (installed < total) {
  failed <- names(results)[!unlist(results)]
  cat("\nFailed packages:\n")
  for (pkg in failed) {
    cat(sprintf("  - %s\n", pkg))
  }
  cat("\nYou may need to install these manually.\n")
}

# Version info
cat("\nKey package versions:\n")
if (requireNamespace("Seurat", quietly = TRUE)) {
  cat(sprintf("  Seurat: %s\n", packageVersion("Seurat")))
}
if (requireNamespace("SeuratData", quietly = TRUE)) {
  cat(sprintf("  SeuratData: %s\n", packageVersion("SeuratData")))
}
cat(sprintf("  R: %s\n", R.version.string))

# Dataset status
if (!minimal && requireNamespace("SeuratData", quietly = TRUE)) {
  cat("\nDemo datasets:\n")
  if (requireNamespace("pbmc3k.SeuratData", quietly = TRUE)) {
    cat("  ✓ pbmc3k\n")
  } else {
    cat("  ✗ pbmc3k (not installed)\n")
  }
  if (requireNamespace("ifnb.SeuratData", quietly = TRUE)) {
    cat("  ✓ ifnb\n")
  } else {
    cat("  ✗ ifnb (not installed)\n")
  }
}

cat("\nSetup complete!\n")
cat("Run 'Rscript tests/run_demo_tests.R --quick' to verify installation.\n\n")
