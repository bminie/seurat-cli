#!/usr/bin/env Rscript
# =============================================================================
# Install the Seurat release seurat-cli needs (>= seurat_min_version in
# utils/common.R) from CRAN, when the installed copy is older.
#
# conda-forge's r-seurat can lag behind CRAN. In a conda environment this
# script builds Seurat from source against the environment's R, on top of the
# r-seurat that conda installed (which supplies nearly all dependencies).
# Once conda-forge has the new release, this script has nothing to do.
#
# Usage:
#   Rscript install_seurat.R           # install or upgrade if needed
#   Rscript install_seurat.R --check   # only report; exit 1 if too old
# =============================================================================

args <- commandArgs(trailingOnly = TRUE)
check_only <- "--check" %in% args

get_script_dir <- function() {
  file_arg <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
  if (length(file_arg) > 0) dirname(normalizePath(sub("^--file=", "", file_arg))) else "."
}
source(file.path(get_script_dir(), "utils", "common.R"))
min_version <- seurat_min_version

installed_version <- function(pkg) {
  # Read DESCRIPTION directly so an already-loaded old copy doesn't mislead us
  path <- system.file("DESCRIPTION", package = pkg)
  if (!nzchar(path)) return(NULL)
  package_version(read.dcf(path, fields = "Version")[1, 1])
}

current <- installed_version("Seurat")
cat(sprintf("Installed Seurat: %s (need >= %s)\n",
            if (is.null(current)) "none" else format(current), min_version))

if (!is.null(current) && current >= min_version) {
  cat("Seurat is up to date.\n")
  quit(status = 0)
}
if (check_only) quit(status = 1)

# -----------------------------------------------------------------------------
# Where to install from, and whether we must compile
# -----------------------------------------------------------------------------

repos <- getOption("repos")
if (is.null(repos) || !nzchar(repos["CRAN"]) || repos["CRAN"] == "@CRAN@") {
  repos <- c(CRAN = "https://cloud.r-project.org")
}

conda_prefix <- Sys.getenv("CONDA_PREFIX")
in_conda <- nzchar(conda_prefix) &&
  startsWith(normalizePath(R.home()), normalizePath(conda_prefix))

# CRAN's binary packages are built for CRAN's R, not conda's, so conda
# environments build from source. Elsewhere use the platform default
# (binaries on macOS and Windows, source on Linux).
pkg_type <- if (in_conda) "source" else getOption("pkgType")
needs_compiler <- pkg_type == "source" || .Platform$pkgType == "source"

if (in_conda) {
  cat(sprintf("conda environment: %s\n", conda_prefix))
}

if (needs_compiler) {
  r_bin <- file.path(R.home("bin"), "R")
  cxx <- suppressWarnings(system2(r_bin, c("CMD", "config", "CXX"), stdout = TRUE))
  cxx_cmd <- strsplit(trimws(cxx[1]), "\\s+")[[1]][1]
  if (!nzchar(Sys.which(cxx_cmd))) {
    hint <- if (in_conda) {
      sprintf("mamba install -p %s -c conda-forge compilers", conda_prefix)
    } else if (Sys.info()[["sysname"]] == "Darwin") {
      "xcode-select --install"
    } else {
      "install your distribution's build tools (e.g. apt-get install build-essential)"
    }
    stop(sprintf("Seurat must be compiled, but the C++ compiler '%s' was not found.\nInstall one with: %s",
                 cxx_cmd, hint), call. = FALSE)
  }
}

# -----------------------------------------------------------------------------
# Dependencies that are missing or older than the new Seurat requires
# -----------------------------------------------------------------------------

db <- utils::available.packages(repos = repos, type = "source")
if (!"Seurat" %in% rownames(db)) stop("Seurat not found in ", repos["CRAN"], call. = FALSE)
cat(sprintf("Installing Seurat %s from %s\n", db["Seurat", "Version"], repos["CRAN"]))

parse_deps <- function(field) {
  if (is.na(field)) return(list())
  lapply(trimws(strsplit(gsub("\\s+", " ", field), ",")[[1]]), function(x) {
    m <- regmatches(x, regexec("^([A-Za-z0-9.]+)\\s*(\\(\\s*>=\\s*([^)]+)\\))?", x))[[1]]
    list(pkg = m[2], min = if (nzchar(m[4])) trimws(m[4]) else NA_character_)
  })
}
base_pkgs <- rownames(utils::installed.packages(priority = "base"))
deps <- unlist(lapply(c("Depends", "Imports", "LinkingTo"),
                      function(f) parse_deps(db["Seurat", f])), recursive = FALSE)
to_install <- character()
for (d in deps) {
  if (d$pkg %in% c("R", base_pkgs)) next
  have <- installed_version(d$pkg)
  if (is.null(have) || (!is.na(d$min) && have < d$min)) {
    to_install <- c(to_install, d$pkg)
    cat(sprintf("  dependency %s: %s (need %s)\n", d$pkg,
                if (is.null(have)) "missing" else format(have),
                if (is.na(d$min)) "any" else paste(">=", d$min)))
  }
}
to_install <- unique(to_install)
if (length(to_install) > 0 && in_conda) {
  cat("  These come from CRAN now. Where conda-forge has builds for your platform,\n",
      "  you can use those instead:\n", sep = "")
  cat(sprintf("    mamba install -p %s -c conda-forge %s\n", conda_prefix,
              paste0("r-", tolower(to_install), collapse = " ")))
}

# -----------------------------------------------------------------------------
# Install
# -----------------------------------------------------------------------------

options(Ncpus = max(1L, parallel::detectCores() - 1L))
utils::install.packages(c(to_install, "Seurat"), repos = repos, type = pkg_type,
                        dependencies = NA)

# Check in a fresh R process, since this one may still have the old Seurat loaded
rscript <- file.path(R.home("bin"), "Rscript")
new_version <- suppressWarnings(system2(
  rscript, c("-e", shQuote("cat(format(packageVersion('Seurat')))")), stdout = TRUE))
new_version <- tail(new_version, 1)
if (length(new_version) == 0 || is.na(package_version(new_version, strict = FALSE)) ||
    package_version(new_version) < min_version) {
  stop("Seurat installation failed; see the messages above.", call. = FALSE)
}
cat(sprintf("\nSeurat %s installed.\n", new_version))
if (in_conda) {
  cat("Note: conda still lists the old r-seurat. Updating that package with conda\n",
      "will replace this copy; run this script again afterwards if needed.\n", sep = "")
}
