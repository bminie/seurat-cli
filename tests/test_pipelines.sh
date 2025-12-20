#!/bin/bash
# =============================================================================
# Test Seurat CLI Workflow Pipelines
# =============================================================================
#
# This script tests the Snakemake and Nextflow pipelines using demo data.
# It automatically manages the conda environment, creating it if needed.
#
# Usage:
#   ./tests/test_pipelines.sh                    # Test both pipelines
#   ./tests/test_pipelines.sh snakemake          # Test Snakemake only
#   ./tests/test_pipelines.sh nextflow           # Test Nextflow only
#   ./tests/test_pipelines.sh --dry-run          # Dry run only (no execution)
#   ./tests/test_pipelines.sh --check            # Check dependencies only
#   ./tests/test_pipelines.sh snakemake --dry-run # Dry run Snakemake only
#
# Options:
#   --dry-run    Only run dry runs (validate configuration without execution)
#   --check      Only check if dependencies are installed
#   --help       Show this help message
#
# =============================================================================

set -e

# Colors for output
RED='\033[0;31m'
GREEN='\033[0;32m'
YELLOW='\033[1;33m'
BLUE='\033[0;34m'
BOLD='\033[1m'
NC='\033[0m' # No Color

# Get script directory
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PROJECT_DIR="$(dirname "$SCRIPT_DIR")"

# Test output directories
SNAKEMAKE_OUTPUT="${PROJECT_DIR}/test_output/snakemake_test"
NEXTFLOW_OUTPUT="${PROJECT_DIR}/test_output/nextflow_test"

# Environment name
CONDA_ENV_NAME="seurat-cli"

# Default options
TEST_TARGET="both"
DRY_RUN_ONLY=false
CHECK_ONLY=false
USE_CONDA_RUN=false

# -----------------------------------------------------------------------------
# Help Function
# -----------------------------------------------------------------------------

show_help() {
    echo ""
    echo "Seurat CLI Pipeline Test Script"
    echo ""
    echo "USAGE:"
    echo "    ./tests/test_pipelines.sh [TARGET] [OPTIONS]"
    echo ""
    echo "TARGETS:"
    echo "    snakemake    Test Snakemake pipeline only"
    echo "    nextflow     Test Nextflow pipeline only"
    echo "    both         Test both pipelines (default)"
    echo ""
    echo "OPTIONS:"
    echo "    --dry-run    Only run dry runs (validate configuration without execution)"
    echo "    --check      Only check if dependencies are installed"
    echo "    --help       Show this help message"
    echo ""
    echo "EXAMPLES:"
    echo "    ./tests/test_pipelines.sh                      # Test both pipelines"
    echo "    ./tests/test_pipelines.sh snakemake            # Test Snakemake only"
    echo "    ./tests/test_pipelines.sh nextflow             # Test Nextflow only"
    echo "    ./tests/test_pipelines.sh --dry-run            # Dry run both pipelines"
    echo "    ./tests/test_pipelines.sh snakemake --dry-run  # Dry run Snakemake only"
    echo "    ./tests/test_pipelines.sh --check              # Check dependencies only"
    echo ""
    echo "ENVIRONMENT:"
    echo "    This script automatically manages the conda environment '${CONDA_ENV_NAME}'."
    echo "    If the environment doesn't exist, it will be created."
    echo "    If dependencies are missing, they will be installed."
    echo ""
    echo "REQUIREMENTS:"
    echo "    - Conda or Mamba must be installed"
    echo "    - environment.yml must exist in the project root"
    echo ""
}

# -----------------------------------------------------------------------------
# Parse Arguments
# -----------------------------------------------------------------------------

parse_args() {
    while [[ $# -gt 0 ]]; do
        case "$1" in
            snakemake|nextflow|both)
                TEST_TARGET="$1"
                shift
                ;;
            --dry-run)
                DRY_RUN_ONLY=true
                shift
                ;;
            --check)
                CHECK_ONLY=true
                shift
                ;;
            --help|-h)
                show_help
                exit 0
                ;;
            *)
                echo -e "${RED}Unknown option: $1${NC}"
                echo "Use --help for usage information"
                exit 1
                ;;
        esac
    done
}

# -----------------------------------------------------------------------------
# Conda/Mamba Environment Functions
# -----------------------------------------------------------------------------

# Detect package manager (mamba preferred over conda, Miniforge preferred)
detect_pkg_manager() {
    # Check for Miniforge mamba first (typically has newer version)
    local miniforge_mamba="$HOME/miniforge3/bin/mamba"
    if [[ -x "$miniforge_mamba" ]]; then
        echo "$miniforge_mamba"
        return
    fi

    # Check for mamba in PATH
    if command -v mamba &> /dev/null; then
        echo "mamba"
        return
    fi

    # Check for Miniforge conda
    local miniforge_conda="$HOME/miniforge3/bin/conda"
    if [[ -x "$miniforge_conda" ]]; then
        echo "$miniforge_conda"
        return
    fi

    # Fall back to conda in PATH
    if command -v conda &> /dev/null; then
        echo "conda"
        return
    fi

    echo ""
}

# Check if conda environment exists
env_exists() {
    local env_name="$1"
    local pkg_manager="${2:-conda}"
    # Match env name followed by whitespace (handles indentation) or at end of path
    "$pkg_manager" env list 2>/dev/null | grep -qE "(^|[[:space:]])${env_name}[[:space:]]" || \
    "$pkg_manager" env list 2>/dev/null | grep -q "/${env_name}$"
}

# Check if we're already in the target environment
in_target_env() {
    [[ "$CONDA_DEFAULT_ENV" == "$CONDA_ENV_NAME" ]]
}

# Get conda base path for sourcing
get_conda_base() {
    if [[ -n "$CONDA_EXE" ]]; then
        dirname "$(dirname "$CONDA_EXE")"
    elif command -v conda &> /dev/null; then
        conda info --base 2>/dev/null
    else
        echo ""
    fi
}

# Create the conda environment
create_environment() {
    local pkg_manager="$1"
    local pkg_name=$(basename "$pkg_manager")

    echo -e "${YELLOW}Creating conda environment '${CONDA_ENV_NAME}'...${NC}"

    if [[ "$pkg_name" == "mamba" ]]; then
        echo "Using mamba - this typically takes 2-5 minutes."
    else
        echo -e "${YELLOW}Using conda - this may take 30+ minutes or fail.${NC}"
        echo "Consider installing mamba for faster environment creation:"
        echo -e "  ${BLUE}conda install -n base -c conda-forge mamba${NC}"
    fi
    echo ""

    cd "$PROJECT_DIR"

    "$pkg_manager" env create -f environment.yml

    if [[ $? -eq 0 ]]; then
        echo ""
        echo -e "${GREEN}Environment '${CONDA_ENV_NAME}' created successfully!${NC}"
        return 0
    else
        echo -e "${RED}Failed to create environment${NC}"
        return 1
    fi
}

# Run a command in the conda environment
run_in_env() {
    local pkg_manager="$1"
    shift
    local cmd="$@"

    if in_target_env; then
        # Already in the environment, run directly
        eval "$cmd"
    else
        # Use conda run to execute in the environment
        "$pkg_manager" run -n "$CONDA_ENV_NAME" bash -c "$cmd"
    fi
}

# Check if a command exists in the environment
check_command_in_env() {
    local pkg_manager="$1"
    local cmd="$2"

    if in_target_env; then
        command -v "$cmd" &> /dev/null
    else
        "$pkg_manager" run -n "$CONDA_ENV_NAME" which "$cmd" &> /dev/null 2>&1
    fi
}

# -----------------------------------------------------------------------------
# Environment Setup and Validation
# -----------------------------------------------------------------------------

setup_environment() {
    echo -e "${BOLD}Checking Environment${NC}"
    echo "-------------------------------------------------------------"

    # Check for conda/mamba
    local pkg_manager=$(detect_pkg_manager)
    local pkg_name=$(basename "$pkg_manager")

    if [[ -z "$pkg_manager" ]]; then
        echo -e "${RED}Error: Neither conda nor mamba is installed.${NC}"
        echo ""
        echo "Please install Miniforge (recommended) or Miniconda:"
        echo -e "  ${BLUE}https://github.com/conda-forge/miniforge#miniforge3${NC}"
        echo ""
        return 1
    fi

    # Warn if using conda instead of mamba
    if [[ "$pkg_name" == "conda" ]]; then
        echo -e "  Package manager.... ${YELLOW}${pkg_name}${NC}"
        echo ""
        echo -e "  ${YELLOW}Warning: Using conda instead of mamba.${NC}"
        echo -e "  Environment creation may take 30+ minutes or fail."
        echo -e "  We strongly recommend installing mamba for faster dependency solving:"
        echo -e "    ${BLUE}conda install -n base -c conda-forge mamba${NC}"
        echo ""
    else
        echo -e "  Package manager.... ${GREEN}${pkg_name}${NC}"
    fi

    # Check if environment exists
    echo -n "  Environment........ "
    if env_exists "$CONDA_ENV_NAME" "$pkg_manager"; then
        echo -e "${GREEN}exists${NC}"
    else
        echo -e "${YELLOW}not found${NC}"
        echo ""

        # Check if environment.yml exists
        if [[ ! -f "${PROJECT_DIR}/environment.yml" ]]; then
            echo -e "${RED}Error: environment.yml not found in ${PROJECT_DIR}${NC}"
            return 1
        fi

        # Ask to create environment
        read -p "  Would you like to create the '${CONDA_ENV_NAME}' environment now? [Y/n] " -n 1 -r
        echo ""

        if [[ ! $REPLY =~ ^[Nn]$ ]]; then
            echo ""
            create_environment "$pkg_manager" || return 1
            echo ""
        else
            echo ""
            echo "To create the environment manually, run:"
            echo -e "  ${BLUE}make env${NC}"
            echo ""
            return 1
        fi
    fi

    # Check if we're in the environment or need to use conda run
    echo -n "  Active environment. "
    if in_target_env; then
        echo -e "${GREEN}${CONDA_ENV_NAME} (active)${NC}"
        USE_CONDA_RUN=false
    else
        echo -e "${YELLOW}${CONDA_DEFAULT_ENV:-none}${NC} (will use '${pkg_name} run')"
        USE_CONDA_RUN=true
    fi

    echo ""

    # Store package manager for later use
    PKG_MANAGER="$pkg_manager"
    return 0
}

# -----------------------------------------------------------------------------
# Dependency Checking Functions
# -----------------------------------------------------------------------------

check_dependencies() {
    echo -e "${BOLD}Checking Dependencies${NC}"
    echo "-------------------------------------------------------------"

    local missing_deps=()
    local pkg_manager="$PKG_MANAGER"

    # Check R
    echo -n "  R.................. "
    if check_command_in_env "$pkg_manager" "Rscript"; then
        local r_version
        if in_target_env; then
            r_version=$(Rscript --version 2>&1 | head -1)
        else
            r_version=$("$pkg_manager" run -n "$CONDA_ENV_NAME" Rscript --version 2>&1 | head -1)
        fi
        echo -e "${GREEN}OK${NC} (${r_version})"

        # Check R packages
        echo -n "  R packages......... "
        local r_check_cmd="Rscript -e \"
            required <- c('Seurat', 'dplyr', 'ggplot2', 'optparse')
            missing <- required[!sapply(required, requireNamespace, quietly = TRUE)]
            if (length(missing) > 0) {
                cat('MISSING:', paste(missing, collapse = ', '))
                quit(status = 1)
            }
            cat('OK')
        \""

        local r_result
        if in_target_env; then
            r_result=$(eval "$r_check_cmd" 2>/dev/null)
        else
            r_result=$("$pkg_manager" run -n "$CONDA_ENV_NAME" bash -c "$r_check_cmd" 2>/dev/null)
        fi

        if [[ "$r_result" == "OK" ]]; then
            echo -e "${GREEN}OK${NC}"
        else
            echo -e "${RED}MISSING${NC}"
            echo -e "    ${YELLOW}${r_result}${NC}"
            missing_deps+=("R packages")
        fi
    else
        echo -e "${RED}NOT FOUND${NC}"
        missing_deps+=("R")
    fi

    # Check Snakemake (only if testing snakemake or both)
    if [[ "$TEST_TARGET" == "snakemake" || "$TEST_TARGET" == "both" ]]; then
        echo -n "  Snakemake.......... "
        if check_command_in_env "$pkg_manager" "snakemake"; then
            local sm_version
            if in_target_env; then
                sm_version=$(snakemake --version 2>&1)
            else
                sm_version=$("$pkg_manager" run -n "$CONDA_ENV_NAME" snakemake --version 2>&1)
            fi
            echo -e "${GREEN}OK${NC} (v${sm_version})"
        else
            echo -e "${RED}NOT FOUND${NC}"
            missing_deps+=("snakemake")
        fi
    fi

    # Check Nextflow (only if testing nextflow or both)
    if [[ "$TEST_TARGET" == "nextflow" || "$TEST_TARGET" == "both" ]]; then
        echo -n "  Nextflow........... "
        if check_command_in_env "$pkg_manager" "nextflow"; then
            local nf_version
            if in_target_env; then
                nf_version=$(nextflow -version 2>&1 | grep -oE 'version [0-9.]+' | head -1 || echo "unknown")
            else
                nf_version=$("$pkg_manager" run -n "$CONDA_ENV_NAME" nextflow -version 2>&1 | grep -oE 'version [0-9.]+' | head -1 || echo "unknown")
            fi
            echo -e "${GREEN}OK${NC} (${nf_version})"
        else
            echo -e "${RED}NOT FOUND${NC}"
            missing_deps+=("nextflow")
        fi
    fi

    echo ""

    # Handle missing dependencies
    if [[ ${#missing_deps[@]} -gt 0 ]]; then
        echo -e "${YELLOW}Missing dependencies: ${missing_deps[*]}${NC}"
        echo ""
        echo "The environment may need to be updated. Try:"
        echo -e "  ${BLUE}${pkg_manager} env update -n ${CONDA_ENV_NAME} -f environment.yml${NC}"
        echo ""
        return 1
    else
        echo -e "${GREEN}All dependencies satisfied!${NC}"
        return 0
    fi
}

# -----------------------------------------------------------------------------
# Install R Extras (packages not available in conda)
# -----------------------------------------------------------------------------

install_r_extras() {
    echo -e "${BOLD}Checking R Extras (GitHub packages & demo data)${NC}"
    echo "-------------------------------------------------------------"

    local pkg_manager="$PKG_MANAGER"

    # R script to check and install extras
    local r_install_script='
        # Function to install if missing
        install_if_missing <- function(pkg, install_func) {
            if (!requireNamespace(pkg, quietly = TRUE)) {
                cat(sprintf("  Installing %s...\n", pkg))
                tryCatch({
                    install_func()
                    cat(sprintf("  %s installed successfully\n", pkg))
                    return(TRUE)
                }, error = function(e) {
                    cat(sprintf("  Warning: Failed to install %s: %s\n", pkg, e$message))
                    return(FALSE)
                })
            } else {
                cat(sprintf("  %s: already installed\n", pkg))
                return(TRUE)
            }
        }

        # Ensure remotes is available
        if (!requireNamespace("remotes", quietly = TRUE)) {
            install.packages("remotes", repos = "https://cloud.r-project.org", quiet = TRUE)
        }

        # Install SeuratData (required for demo datasets)
        install_if_missing("SeuratData", function() {
            remotes::install_github("satijalab/seurat-data", quiet = TRUE, upgrade = "never")
        })

        # Install presto (faster marker finding)
        install_if_missing("presto", function() {
            remotes::install_github("immunogenomics/presto", quiet = TRUE, upgrade = "never")
        })

        # Install demo dataset if SeuratData is available
        if (requireNamespace("SeuratData", quietly = TRUE)) {
            library(SeuratData)
            available <- AvailableData()
            if ("pbmc3k" %in% available$Dataset) {
                installed <- InstalledData()
                if (!"pbmc3k" %in% installed$Dataset) {
                    cat("  Installing pbmc3k demo dataset...\n")
                    tryCatch({
                        InstallData("pbmc3k", quiet = TRUE)
                        cat("  pbmc3k dataset installed successfully\n")
                    }, error = function(e) {
                        cat(sprintf("  Warning: Failed to install pbmc3k: %s\n", e$message))
                    })
                } else {
                    cat("  pbmc3k dataset: already installed\n")
                }
            }
        }

        cat("  R extras check complete\n")
    '

    echo "  Checking for SeuratData, presto, and demo datasets..."

    if in_target_env; then
        Rscript -e "$r_install_script" 2>/dev/null
    else
        "$pkg_manager" run -n "$CONDA_ENV_NAME" Rscript -e "$r_install_script" 2>/dev/null
    fi

    echo ""
    return 0
}

# -----------------------------------------------------------------------------
# Snakemake Test Functions
# -----------------------------------------------------------------------------

test_snakemake_dry_run() {
    echo -e "${YELLOW}Snakemake Dry Run...${NC}"

    # Clean and create test directory
    rm -rf "$SNAKEMAKE_OUTPUT"
    mkdir -p "$SNAKEMAKE_OUTPUT"

    # Create a minimal test config
    cat > "${SNAKEMAKE_OUTPUT}/test_config.yaml" << EOF
workflow: "sctransform"
scripts_dir: "${PROJECT_DIR}/scripts"
output_dir: "${SNAKEMAKE_OUTPUT}/results"
input_path: "DEMO"
species: "human"
threads: 2

sctransform:
  n_variable_features: 3000
  resolution: 0.5
  dims_use: "1:30"
  extra_args: "--demo"

de:
  mode: "all_markers"
  test_use: "wilcox"
  logfc_threshold: 0.25
  min_pct: 0.1
  extra_args: "--demo"

visualization:
  reduction: "umap"
  extra_args: "--demo"
EOF

    echo "  Validating Snakemake configuration..."

    local sm_cmd="cd '${PROJECT_DIR}/workflows/snakemake' && snakemake --configfile '${SNAKEMAKE_OUTPUT}/test_config.yaml' -n --quiet"

    if run_in_env "$PKG_MANAGER" "$sm_cmd" 2>&1; then
        echo -e "  ${GREEN}✓ Configuration valid${NC}"
        echo ""
        echo "  Planned execution:"
        run_in_env "$PKG_MANAGER" "$sm_cmd" 2>&1 | head -20
        echo ""
        return 0
    else
        echo -e "  ${RED}✗ Configuration invalid${NC}"
        run_in_env "$PKG_MANAGER" "cd '${PROJECT_DIR}/workflows/snakemake' && snakemake --configfile '${SNAKEMAKE_OUTPUT}/test_config.yaml' -n" 2>&1
        return 1
    fi
}

test_snakemake_full() {
    echo -e "${YELLOW}Snakemake Full Test...${NC}"

    # Run dry run first if config doesn't exist
    if [[ ! -f "${SNAKEMAKE_OUTPUT}/test_config.yaml" ]]; then
        test_snakemake_dry_run || return 1
    fi

    echo "  Running Snakemake pipeline (this may take several minutes)..."

    local sm_cmd="cd '${PROJECT_DIR}/workflows/snakemake' && snakemake --configfile '${SNAKEMAKE_OUTPUT}/test_config.yaml' --cores 2"

    if run_in_env "$PKG_MANAGER" "$sm_cmd"; then
        echo -e "  ${GREEN}✓ Pipeline completed${NC}"
    else
        echo -e "  ${RED}✗ Pipeline failed${NC}"
        return 1
    fi

    # Check outputs
    echo "  Checking outputs..."
    local expected_files=(
        "results/01_sctransform/seurat_sctransform.rds"
        "results/03_differential_expression/de_results_all.csv"
        "results/visualization/visualization_summary.md"
    )

    local all_found=true
    for f in "${expected_files[@]}"; do
        if [[ -f "${SNAKEMAKE_OUTPUT}/${f}" ]]; then
            echo -e "    ${GREEN}✓${NC} ${f}"
        else
            echo -e "    ${RED}✗${NC} ${f} (not found)"
            all_found=false
        fi
    done

    if $all_found; then
        return 0
    else
        return 1
    fi
}

# -----------------------------------------------------------------------------
# Nextflow Test Functions
# -----------------------------------------------------------------------------

test_nextflow_dry_run() {
    echo -e "${YELLOW}Nextflow Dry Run...${NC}"

    echo "  Validating Nextflow configuration..."

    local nf_cmd="cd '${PROJECT_DIR}/workflows/nextflow' && nextflow run main.nf -profile test --demo --outdir '${NEXTFLOW_OUTPUT}/results' -preview"

    if run_in_env "$PKG_MANAGER" "$nf_cmd" 2>&1; then
        echo -e "  ${GREEN}✓ Configuration valid${NC}"
        return 0
    else
        echo -e "  ${RED}✗ Configuration invalid${NC}"
        return 1
    fi
}

test_nextflow_full() {
    echo -e "${YELLOW}Nextflow Full Test...${NC}"

    # Clean previous test output
    rm -rf "$NEXTFLOW_OUTPUT"
    mkdir -p "$NEXTFLOW_OUTPUT"

    echo "  Running Nextflow pipeline (this may take several minutes)..."

    local nf_cmd="cd '${PROJECT_DIR}/workflows/nextflow' && nextflow run main.nf \
        -profile test \
        --demo \
        --outdir '${NEXTFLOW_OUTPUT}/results' \
        --workflow sctransform \
        --threads 2 \
        -with-report '${NEXTFLOW_OUTPUT}/report.html' \
        -with-timeline '${NEXTFLOW_OUTPUT}/timeline.html'"

    if run_in_env "$PKG_MANAGER" "$nf_cmd"; then
        echo -e "  ${GREEN}✓ Pipeline completed${NC}"
    else
        echo -e "  ${RED}✗ Pipeline failed${NC}"
        return 1
    fi

    # Check outputs
    echo "  Checking outputs..."
    local expected_files=(
        "results/01_sctransform/seurat_sctransform.rds"
        "results/03_differential_expression/de_results_all.csv"
        "results/visualization/visualization_summary.md"
    )

    local all_found=true
    for f in "${expected_files[@]}"; do
        if [[ -f "${NEXTFLOW_OUTPUT}/${f}" ]]; then
            echo -e "    ${GREEN}✓${NC} ${f}"
        else
            echo -e "    ${RED}✗${NC} ${f} (not found)"
            all_found=false
        fi
    done

    if $all_found; then
        return 0
    else
        return 1
    fi
}

# -----------------------------------------------------------------------------
# Main
# -----------------------------------------------------------------------------

parse_args "$@"

echo ""
echo "============================================================="
echo "  Seurat CLI Pipeline Tests"
echo "============================================================="
echo ""

# Setup and validate environment
if ! setup_environment; then
    echo ""
    exit 1
fi

# Check dependencies
if ! check_dependencies; then
    echo ""
    exit 1
fi

# Install R extras (SeuratData, presto, demo datasets) - only for full tests
if [[ "$DRY_RUN_ONLY" != true && "$CHECK_ONLY" != true ]]; then
    install_r_extras
fi

# If --check flag, exit after dependency check
if [[ "$CHECK_ONLY" == true ]]; then
    echo ""
    exit 0
fi

echo ""

# Track results
SNAKEMAKE_RESULT=0
NEXTFLOW_RESULT=0

# Run tests based on target and mode
if [[ "$DRY_RUN_ONLY" == true ]]; then
    # Dry run mode
    echo -e "${BOLD}Running Dry Runs Only${NC}"
    echo "-------------------------------------------------------------"
    echo ""

    if [[ "$TEST_TARGET" == "snakemake" || "$TEST_TARGET" == "both" ]]; then
        test_snakemake_dry_run || SNAKEMAKE_RESULT=1
        echo ""
    fi

    if [[ "$TEST_TARGET" == "nextflow" || "$TEST_TARGET" == "both" ]]; then
        test_nextflow_dry_run || NEXTFLOW_RESULT=1
        echo ""
    fi
else
    # Full test mode
    echo -e "${BOLD}Running Full Tests${NC}"
    echo "-------------------------------------------------------------"
    echo ""

    if [[ "$TEST_TARGET" == "snakemake" || "$TEST_TARGET" == "both" ]]; then
        test_snakemake_full || SNAKEMAKE_RESULT=1
        echo ""
    fi

    if [[ "$TEST_TARGET" == "nextflow" || "$TEST_TARGET" == "both" ]]; then
        test_nextflow_full || NEXTFLOW_RESULT=1
        echo ""
    fi
fi

# Summary
echo "============================================================="
echo "  Test Summary"
echo "============================================================="

if [[ "$TEST_TARGET" == "snakemake" || "$TEST_TARGET" == "both" ]]; then
    if [[ $SNAKEMAKE_RESULT -eq 0 ]]; then
        echo -e "  Snakemake: ${GREEN}PASSED${NC}"
    else
        echo -e "  Snakemake: ${RED}FAILED${NC}"
    fi
fi

if [[ "$TEST_TARGET" == "nextflow" || "$TEST_TARGET" == "both" ]]; then
    if [[ $NEXTFLOW_RESULT -eq 0 ]]; then
        echo -e "  Nextflow:  ${GREEN}PASSED${NC}"
    else
        echo -e "  Nextflow:  ${RED}FAILED${NC}"
    fi
fi

echo ""

# Exit with error if any test failed
exit $((SNAKEMAKE_RESULT + NEXTFLOW_RESULT))
