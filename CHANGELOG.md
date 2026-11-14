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
- `--no_dim_plots`, `--no_feature_plots`, `--no_violin_plots`, `--no_dot_plot`, `--no_heatmap`, and `--no_label_clusters` flags for `06_visualization.R` (these plots were on with no way to turn them off)
- `load_demo_dataset()` helper in `utils/common.R`, used by every script's `--demo` mode
- Demo Tests CI workflow (`.github/workflows/demo-tests.yml`): builds `environment-r.yml`, installs the demo datasets, and runs `tests/run_demo_tests.R` (the first CI job that runs the analyses)
- Tier 1 validation checks that every PNG figure has content (catches blank images); adds `r-png` to the conda environments
- Nextflow CI runs `nextflow lint` (strict syntax)
- `--max_features` and `--max_mt` options for `02_sctransform.R` to filter cells before SCTransform (off by default, matching the SCTransform vignette)
- Seurat 6 support:
  - `07_sample_level.R`: sample-level analysis with scSLIDE (landmark density per sample, diffusion map trajectories, pseudobulk trajectory DE per cell type, annotated heatmaps). Demo uses ifnb with donor IDs (16 samples)
  - `08_spatial.R`: 10x Xenium and Atera imaging data (QC, clustering, markers, spatial plots, transcript positions for selected genes, zoomed view with cell outlines, Atera Azimuth annotations)
  - `install_seurat.R` (and `make install-seurat`): installs the required Seurat from CRAN when the installed one (e.g. conda-forge's) is older; `setup.R` runs it first
  - Pseudobulk heatmap in `04_differential_expression.R --pseudobulk` (`PlotAnnotatedHeatmap`)
  - `--umap_init` for scripts 01, 02, 03, 05 (RunUMAP `uwot.init`); `--svd_method` for CCA integration in 03; `--threads` for 05
  - `tests/make_xenium_fixture.R` writes a synthetic Xenium folder; the demo tests now cover scripts 01-08
  - Snakemake and Nextflow `sample_level` (03 → 07, or 07 on an existing RDS) and `spatial` (08) workflows; `tests/test_pipelines.sh` runs both in each pipeline

### Changed
- `tests/run_demo_tests.R` now integrates the validation framework and reports detailed check counts
- `06_visualization.R` now takes the top `--top_markers` genes per cluster from `--marker_file` (interleaved across clusters) instead of every gene in the file
- Nextflow pipeline requires Nextflow >= 24.04; resource caps use the built-in `resourceLimits` instead of a custom `check_max` function
- Snakemake demo config only passes `--demo` to entry-point steps, so downstream steps use upstream outputs; demo can now run the `basic` and `integration` workflows
- Conda environments: removed `bioconductor-batchelor` (no osx-arm64 build), install presto from bioconda, add `r-remotes`, drop the `defaults` channel, require Snakemake >= 8 and Nextflow >= 24.04
- Snakemake container docs use `--sdm apptainer` (Snakemake >= 8) instead of the removed `--use-singularity`
- `make clean` also removes Nextflow/Snakemake work directories and pipeline result folders
- Requires Seurat >= 6.0: scripts stop with install instructions on older versions. Docker image uses `rocker/r-ver:4.6.2` (Seurat 6 is newer than the R 4.5.3 package snapshot). Conda environments add the Seurat 6 build dependencies and the packages for the new features
- `--threads` now uses Seurat's native multithreading (`setThreads()`) instead of `future` multisession workers: threads share one copy of the data, so memory no longer grows with the thread count (SCTransform on pbmc3k: 2.6 GB with 1 thread, 2.8 GB with 4)
- `setup.R` installs ComplexHeatmap and destiny by default and the spatial readers with `--all`; it installs destiny 3.26 from Bioconductor 3.23 when the older release can't install (its `smoother` dependency was archived on CRAN), and sets a CRAN mirror when none is configured

### Fixed
- `03_integration.R` failed with more than one input file: merging Seurat v5 objects already splits the counts by input, and splitting again errored. Layers are now joined before splitting by sample (the demo, which starts from one object, was unaffected)
- Snakemake pipeline silently ran the `basic` workflow for an unknown `workflow` name; it now stops with the list of workflows
- Nextflow pipeline failed to start on Nextflow >= 26 (strict syntax parser): moved top-level `--help`/input checks into the workflow and removed the function from `nextflow.config`
- `max_memory`, `max_cpus`, and `max_time` were never applied (the old `check_max` helper always fell back to the default)
- Nextflow `test` profile resource settings were overridden by the main `process` block; profiles now take precedence
- Nextflow `docker`/`singularity` profiles pointed at host script paths; they now use `/app/scripts` inside the image
- `basic` workflow in both Snakemake and Nextflow expected `seurat_analyzed.rds`, but `01_basic_analysis.R` writes `seurat_object.rds`
- Visualization step in pipelines plotted every significant DE gene (~300 PNGs, >200 MB for the demo)
- `01_basic_analysis.R` saved a blank `06_pca_heatmap.png` (`DimHeatmap` drew with base graphics and returned nothing to save)
- `03_integration.R` ignored `--reference`, `--k_anchor`, and `--dims_use` during integration; they are now passed to CCA/RPCA integration
- `03_integration.R --normalization SCT` integrated as if the data were log-normalized; it now passes `normalization.method = "SCT"`
- `03_integration.R --method FastMNNIntegration` failed with an unclear error; it now checks for SeuratWrappers and batchelor and says how to install them. Unknown methods and normalizations are rejected up front
- `03_integration.R` errored when `--resolutions` did not include `--resolution`, and silently dropped `--conditions` when the count did not match the number of samples
- Non-demo runs of `03_integration.R` required SeuratData (only needed for demo data)
- Demo mode re-ran `InstallData()` on every run because it checked for a package named `pbmc3k` instead of `pbmc3k.SeuratData`
- `05_cell_cycle.R` on SCTransform input overwrote the SCT residuals with `ScaleData` before PCA, and regressed cell cycle with `ScaleData`; it now keeps the residuals and regresses by re-running `SCTransform` with `vars.to.regress`
- `05_cell_cycle.R` used the deprecated `GetAssayData(slot = )`
- Test runner misreported per-tier check counts (Tier 2 counted each check's fields, not checks)
- CI `--help` check could never fail (a non-zero exit only raised a warning); it now fails the job
- `make lint` printed failures but still exited successfully; it now exits non-zero and also checks `tests/*.R`
- `04_differential_expression.R --pseudobulk` without `--pseudobulk_group` silently ran standard DE; it now stops with a clear message, as does `--pseudobulk` outside `between_conditions` mode
- Pseudobulk DE failed to match conditions containing `_` (AggregateExpression renames them with `-`) and could not group by `orig.ident` (AggregateExpression overwrites it); it now checks for at least 2 replicates per condition
- `04_differential_expression.R` checks that `--condition_col`, `--pseudobulk_group`, and the condition values exist before running
- Marker heatmaps in scripts 02, 03, 04, and 06 dropped genes without scaled values; missing genes are now scaled (or get SCT residuals) for plotting
- Docker image could not install Seurat: it used R 4.3.2 but installed from current CRAN, where Matrix needs R >= 4.4. It now uses `rocker/r-ver:4.5.3` with the base image's date-pinned package snapshot, and the build fails if any required R package is missing (install failures only warned before)
- Nextflow `docker` profile failed loading SeuratData: tasks run as the host user, who has no home in the image; the profile now sets `HOME=/tmp`
- Nextflow `test` profile ran 2 threads in a 4 GB limit, which SCTransform exceeds (~5.8 GB) when the limit is enforced (e.g. in Docker); it now uses 1 thread (~2 GB)
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
