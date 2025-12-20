# Makefile for Seurat CLI
# Common tasks for development and testing

.PHONY: help setup test test-quick test-all test-pipelines test-dry-run check-deps clean lint env env-r

help:
	@echo "Seurat CLI - Available Commands"
	@echo ""
	@echo "Setup:"
	@echo "  make env              - Create conda environment (full) [requires mamba - recommended]"
	@echo "  make env-r            - Create conda environment (R-only) [requires mamba - recommended]"
	@echo "  make setup            - Install R packages via setup.R [requires R >= 4.2]"
	@echo "  make setup-all        - Install all R packages including optional [requires R >= 4.2]"
	@echo ""
	@echo "  Note: mamba is strongly recommended over conda for faster environment creation."
	@echo "  Install mamba: conda install -n base -c conda-forge mamba"
	@echo ""
	@echo "Testing:"
	@echo "  make test-quick       - Run quick tests (scripts 1-2)"
	@echo "  make test             - Run all demo tests"
	@echo "  make test-1           - Run specific script test (1-6)"
	@echo "  make test-pipelines   - Test Snakemake/Nextflow pipelines (full)"
	@echo "  make test-snakemake   - Test Snakemake pipeline only"
	@echo "  make test-nextflow    - Test Nextflow pipeline only"
	@echo "  make test-dry-run     - Dry run both pipelines (validate only)"
	@echo "  make check-deps       - Check pipeline dependencies"
	@echo ""
	@echo "Other:"
	@echo "  make lint             - Check R code style"
	@echo "  make clean            - Remove test output"
	@echo ""

# Conda environments
env:
	@if command -v mamba >/dev/null 2>&1; then \
		echo "Creating environment with mamba (recommended)..."; \
		echo "This typically takes 2-5 minutes."; \
		echo ""; \
		mamba env create -f environment.yml; \
	else \
		echo ""; \
		echo "WARNING: mamba not found, using conda instead."; \
		echo "This may take 30+ minutes or fail due to complex dependencies."; \
		echo ""; \
		echo "We strongly recommend installing mamba first:"; \
		echo "  conda install -n base -c conda-forge mamba"; \
		echo ""; \
		echo "Press Ctrl+C to cancel, or wait 5 seconds to continue with conda..."; \
		sleep 5; \
		echo "Creating environment with conda..."; \
		conda env create -f environment.yml; \
	fi
	@echo ""
	@echo "Environment created! Activate with: conda activate seurat-cli"

env-r:
	@if command -v mamba >/dev/null 2>&1; then \
		echo "Creating R-only environment with mamba (recommended)..."; \
		echo "This typically takes 1-3 minutes."; \
		echo ""; \
		mamba env create -f environment-r.yml; \
	else \
		echo ""; \
		echo "WARNING: mamba not found, using conda instead."; \
		echo "This may take 15+ minutes due to complex dependencies."; \
		echo ""; \
		echo "We strongly recommend installing mamba first:"; \
		echo "  conda install -n base -c conda-forge mamba"; \
		echo ""; \
		echo "Press Ctrl+C to cancel, or wait 5 seconds to continue with conda..."; \
		sleep 5; \
		echo "Creating environment with conda..."; \
		conda env create -f environment-r.yml; \
	fi
	@echo ""
	@echo "Environment created! Activate with: conda activate seurat-cli-r"

# Manual R setup
setup:
	Rscript setup.R

setup-all:
	Rscript setup.R --all

# Testing
test-quick:
	Rscript tests/run_demo_tests.R --quick

test:
	Rscript tests/run_demo_tests.R

test-1:
	Rscript tests/run_demo_tests.R --script 1

test-2:
	Rscript tests/run_demo_tests.R --script 2

test-3:
	Rscript tests/run_demo_tests.R --script 3

test-4:
	Rscript tests/run_demo_tests.R --script 4

test-5:
	Rscript tests/run_demo_tests.R --script 5

test-6:
	Rscript tests/run_demo_tests.R --script 6

# Pipeline tests
test-pipelines:
	./tests/test_pipelines.sh

test-snakemake:
	./tests/test_pipelines.sh snakemake

test-nextflow:
	./tests/test_pipelines.sh nextflow

test-dry-run:
	./tests/test_pipelines.sh --dry-run

check-deps:
	./tests/test_pipelines.sh --check

# Linting
lint:
	@echo "Checking R syntax..."
	@for f in scripts/*.R utils/*.R; do \
		echo "Checking $$f..."; \
		Rscript -e "parse('$$f')" > /dev/null 2>&1 && echo "  ✓ OK" || echo "  ✗ FAILED"; \
	done

# Cleanup
clean:
	rm -rf test_output/
	rm -rf output/
	rm -rf results/
	rm -f *.log
	rm -f .RData .Rhistory

# Show help for a specific script
help-1:
	Rscript scripts/01_basic_analysis.R --help

help-2:
	Rscript scripts/02_sctransform.R --help

help-3:
	Rscript scripts/03_integration.R --help

help-4:
	Rscript scripts/04_differential_expression.R --help

help-5:
	Rscript scripts/05_cell_cycle.R --help

help-6:
	Rscript scripts/06_visualization.R --help
