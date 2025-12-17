# Changelog

All notable changes to this project will be documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.0.0/),
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## [Unreleased]

### Added
- Nothing yet

### Changed
- Nothing yet

### Fixed
- Nothing yet

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
