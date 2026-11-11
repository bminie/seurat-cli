# Changelog

All notable changes to this project will be documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.0.0/),
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## [Unreleased]

### Added
- Docker support with Dockerfile for containerized execution
  - Includes all dependencies, demo datasets, and scripts
  - Compatible with Singularity for HPC clusters
- Two-tier validation framework (`tests/validate_outputs.R`)
  - Tier 1: Structural validation (file existence, CSV structure, cluster assignments, Seurat object integrity)
  - Tier 2: Biological validation for demo datasets (expected cluster counts, cell counts, known markers)
- New test runner options: `--skip_validation`, `--validation_only`
- `cluster_assignments.csv` output for cell cycle script (05) for consistency with other scripts
- GitHub Actions workflows for Nextflow and Snakemake CI
- `--max_plot_features` option for `06_visualization.R` (default 36) to cap how many marker-file genes are plotted

### Changed
- `tests/run_demo_tests.R` now integrates the validation framework and reports detailed check counts
- `06_visualization.R` now takes the top `--top_markers` genes per cluster from `--marker_file` (interleaved across clusters) instead of every gene in the file
- Nextflow pipeline requires Nextflow >= 24.04; resource caps use the built-in `resourceLimits` instead of a custom `check_max` function
- Snakemake demo config only passes `--demo` to entry-point steps, so downstream steps use upstream outputs; demo can now run the `basic` and `integration` workflows
- Conda environments: removed `bioconductor-batchelor` (no osx-arm64 build), install presto from bioconda, add `r-remotes`, drop the `defaults` channel, require Snakemake >= 8 and Nextflow >= 24.04
- Snakemake container docs use `--sdm apptainer` (Snakemake >= 8) instead of the removed `--use-singularity`
- `make clean` also removes Nextflow/Snakemake work directories and pipeline result folders

### Fixed
- Nextflow pipeline failed to start on Nextflow >= 26 (strict syntax parser): moved top-level `--help`/input checks into the workflow and removed the function from `nextflow.config`
- `max_memory`, `max_cpus`, and `max_time` were never applied (the old `check_max` helper always fell back to the default)
- Nextflow `test` profile resource settings were overridden by the main `process` block; profiles now take precedence
- Nextflow `docker`/`singularity` profiles pointed at host script paths; they now use `/app/scripts` inside the image
- `basic` workflow in both Snakemake and Nextflow expected `seurat_analyzed.rds`, but `01_basic_analysis.R` writes `seurat_object.rds`
- Visualization step in pipelines plotted every significant DE gene (~300 PNGs, >200 MB for the demo)
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
