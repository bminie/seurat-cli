# Changelog

All notable changes to this project will be documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.0.0/),
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## [Unreleased]

### Added
- Two-tier validation framework (`tests/validate_outputs.R`)
  - Tier 1: Structural validation (file existence, CSV structure, cluster assignments, Seurat object integrity)
  - Tier 2: Biological validation for demo datasets (expected cluster counts, cell counts, known markers)
- New test runner options: `--skip_validation`, `--validation_only`
- `cluster_assignments.csv` output for cell cycle script (05) for consistency with other scripts

### Changed
- `tests/run_demo_tests.R` now integrates the validation framework and reports detailed check counts

### Fixed
- Fixed clustering bug in `scripts/03_integration.R` where `FindNeighbors`/`FindClusters` used incorrect graph names, causing all cells to be assigned to a single cluster instead of the expected 15 clusters. The fix explicitly sets `graph.name` parameter to ensure consistent graph naming across Seurat v5.

## [1.0.0] - 2024-12-16

### Added
- Initial release with 6 analysis scripts

#### Scripts
- `01_basic_analysis.R` - Standard scRNA-seq workflow (QC, normalization, clustering, markers)
- `02_sctransform.R` - SCTransform normalization with v1/v2 support
- `03_integration.R` - Multi-dataset integration (CCA, RPCA, Harmony, FastMNN)
- `04_differential_expression.R` - Flexible DE testing with multiple methods
- `05_cell_cycle.R` - Cell cycle scoring and regression
- `06_visualization.R` - Comprehensive visualization suite

#### Features
- Command-line interface with `optparse` for all scripts
- Auto-detection of input formats (10x MTX, H5, RDS, CSV)
- Demo mode (`--demo` flag) for testing with built-in datasets
- Comprehensive logging with timestamps
- Parameter tracking in output (`parameters.csv`)
- Session info recording (`session_info.txt`)
- Multi-threading support where applicable
- Species-specific handling (human/mouse)

#### Utility Functions (`utils/common.R`)
- Package management and installation helpers
- Data loading for multiple formats
- QC metric calculation
- Output directory management
- Logging infrastructure
- Common argument parsing

#### Documentation
- Comprehensive README with examples
- Parameter tables for each script
- Example workflows (single sample, integration, full pipeline)
- Installation instructions

---

## Version History Summary

| Version | Date | Highlights |
|---------|------|------------|
| 1.0.0 | 2024-12-16 | Initial release with 6 scripts |
