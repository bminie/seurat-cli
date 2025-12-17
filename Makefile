# Makefile for Seurat CLI
# Common tasks for development and testing

.PHONY: help setup test test-quick test-all clean lint

help:
	@echo "Seurat CLI - Available Commands"
	@echo ""
	@echo "  make setup       - Install required R packages"
	@echo "  make setup-all   - Install all packages including optional"
	@echo "  make test-quick  - Run quick tests (scripts 1-2)"
	@echo "  make test        - Run all demo tests"
	@echo "  make test-1      - Run specific script test (1-6)"
	@echo "  make lint        - Check R code style"
	@echo "  make clean       - Remove test output"
	@echo ""

# Setup
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
