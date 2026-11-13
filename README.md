# Seurat CLI

[![Demo Tests](https://github.com/bminie/seurat-cli/actions/workflows/demo-tests.yml/badge.svg?branch=main)](https://github.com/bminie/seurat-cli/actions/workflows/demo-tests.yml)
[![R Syntax Check](https://github.com/bminie/seurat-cli/actions/workflows/r-check.yml/badge.svg?branch=main)](https://github.com/bminie/seurat-cli/actions/workflows/r-check.yml)
[![Nextflow](https://github.com/bminie/seurat-cli/actions/workflows/nextflow.yml/badge.svg?branch=main)](https://github.com/bminie/seurat-cli/actions/workflows/nextflow.yml)
[![Snakemake](https://github.com/bminie/seurat-cli/actions/workflows/snakemake.yml/badge.svg?branch=main)](https://github.com/bminie/seurat-cli/actions/workflows/snakemake.yml)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)

A collection of flexible, parameterized R scripts that implement the major [Seurat](https://satijalab.org/seurat/) vignettes for single-cell RNA-seq analysis.

## Overview

These scripts provide command-line interfaces to common Seurat workflows, allowing users to:
- Run standard analyses with custom parameters
- Process different input formats (10x, CSV, RDS)
- Automate batch processing
- Generate reproducible results with logged parameters

The scripts require **Seurat 6.0 or newer**. Besides the classic vignettes they
cover the Seurat 6 additions: sample-level analysis with scSLIDE
(`07_sample_level.R`) and 10x Xenium/Atera imaging data (`08_spatial.R`).

## Installation

### Prerequisites

Before installing, ensure you have one of the following:

**For Option 1 (Mamba - Recommended):**
- [Miniforge](https://github.com/conda-forge/miniforge#miniforge3) (includes mamba) - **Recommended**
- Or install mamba into existing conda: `conda install -n base -c conda-forge mamba`

> **Note:** We strongly recommend mamba over conda. Mamba uses a faster C++ dependency solver and typically completes environment creation in 2-5 minutes, while conda can take 30+ minutes or fail on complex environments like this one.

**For Option 1 (Conda - Slower):**
- [Miniconda](https://docs.conda.io/en/latest/miniconda.html) or [Anaconda](https://www.anaconda.com/download)
- If using conda, first enable the libmamba solver for faster resolution:
  ```bash
  conda install -n base -c conda-forge conda-libmamba-solver
  conda config --set solver libmamba
  ```

**For Option 2 (Manual):**
- [R](https://cran.r-project.org/) version 4.2 or higher
- (Optional) [RStudio](https://posit.co/download/rstudio-desktop/)

### Clone the Repository

```bash
git clone https://github.com/bminie/seurat-cli.git
cd seurat-cli
```

### Option 1: Mamba/Conda (Recommended)

Create an isolated environment with all dependencies pre-configured:

```bash
# Using mamba (recommended - much faster)
mamba env create -f environment.yml
conda activate seurat-cli

# Or using conda (slower, may take 30+ minutes)
conda env create -f environment.yml
conda activate seurat-cli

# R-only environment (lighter, no workflow tools)
mamba env create -f environment-r.yml   # or: conda env create -f environment-r.yml
conda activate seurat-cli-r

# Then, in the activated environment: Seurat 6, the demo datasets, and the
# packages conda can't provide
Rscript setup.R
```

> **Seurat 6 and conda:** conda-forge can lag behind CRAN; at the time of
> writing it still has Seurat 5.5. `setup.R` runs `install_seurat.R`, which
> builds the current Seurat from CRAN on top of the environment (about a
> minute; the environment includes the compilers and dependencies). Run it on
> its own with `Rscript install_seurat.R` (or `make install-seurat`), or check
> the installed version with `Rscript install_seurat.R --check`. If you later
> update `r-seurat` with conda, run it again. Once conda-forge has Seurat 6 it
> has nothing to do.

### Option 2: Manual Installation

If you prefer to use your existing R installation:

```bash
# Install R dependencies
Rscript setup.R           # Core packages
Rscript setup.R --all     # All packages including optional

# For workflow pipelines, also install:
pip install snakemake                      # Snakemake
curl -s https://get.nextflow.io | bash     # Nextflow
```

Or open in RStudio by double-clicking `seurat-cli.Rproj`.

### Option 3: Docker

Run without installing any dependencies:

```bash
# Build the image (the second tag is the name the Nextflow docker profile uses)
docker build -t seurat-cli -t bminie/seurat-cli:latest .

# On Apple Silicon, build an amd64 image for x86-64 Linux/HPC machines instead.
# Takes ~15 min: Bioconductor packages compile from source under emulation.
docker buildx build --platform linux/amd64 --load -t seurat-cli:amd64 .

# Run with demo data
docker run --rm -v $(pwd)/output:/output seurat-cli \
  Rscript scripts/01_basic_analysis.R --demo --output /output

# Run with your data
docker run --rm -v /path/to/data:/data -v $(pwd)/output:/output seurat-cli \
  Rscript scripts/01_basic_analysis.R --input /data --output /output

# Interactive session
docker run --rm -it -v $(pwd):/work seurat-cli bash

# Get help for any script
docker run --rm seurat-cli Rscript scripts/01_basic_analysis.R --help
```

**Run multiple scripts in sequence (pipeline):**
```bash
docker run --rm -v $(pwd)/output:/output seurat-cli bash -c "\
  Rscript scripts/02_sctransform.R --demo --output /output/01_sct && \
  Rscript scripts/04_differential_expression.R \
    --input /output/01_sct/seurat_sctransform.rds \
    --output /output/02_de && \
  Rscript scripts/06_visualization.R \
    --input /output/01_sct/seurat_sctransform.rds \
    --marker_file /output/02_de/de_results_significant.csv \
    --output /output/03_viz"
```

**Singularity (for HPC clusters):**
```bash
# Pull from Docker Hub (after image is published)
singularity pull docker://bminie/seurat-cli:latest

# Run a single script
singularity exec seurat-cli_latest.sif \
  Rscript /app/scripts/01_basic_analysis.R --demo --output ./output

# Run with your data
singularity exec --bind /path/to/data:/data,./output:/output seurat-cli_latest.sif \
  Rscript /app/scripts/01_basic_analysis.R --input /data --output /output
```

**Use with workflow pipelines:**

> The image is not on Docker Hub yet. Build it locally first (see above) so the
> Nextflow `docker` profile can find `bminie/seurat-cli:latest`. The `singularity`
> profile and Snakemake `--sdm apptainer` pull from Docker Hub, so they need the
> image to be published.

```bash
# Nextflow with Docker
cd workflows/nextflow
nextflow run main.nf --demo -profile docker

# Nextflow with Singularity (HPC)
nextflow run main.nf --demo -profile singularity

# Snakemake with Apptainer/Singularity
cd workflows/snakemake
snakemake --configfile config_demo.yaml --cores 4 --sdm apptainer
```

## Quick Start

```bash
# Basic PBMC-style analysis
Rscript scripts/01_basic_analysis.R --input /path/to/10x/data --output ./results

# Run with demo data
Rscript scripts/01_basic_analysis.R --demo --output ./demo_results

# SCTransform normalization
Rscript scripts/02_sctransform.R --input /path/to/data --output ./sct_results

# View all options for any script
Rscript scripts/01_basic_analysis.R --help
```

## Scripts

### 1. Basic Analysis (`01_basic_analysis.R`)

Standard scRNA-seq workflow including QC, normalization, clustering, and marker identification.

**Key Parameters:**
| Parameter | Default | Description |
|-----------|---------|-------------|
| `--input` | required | Input file/directory |
| `--output` | `output/pbmc3k` | Output directory |
| `--min_features` | 200 | Min features per cell |
| `--max_features` | 2500 | Max features per cell |
| `--max_mt` | 5 | Max % mitochondrial |
| `--n_variable_features` | 2000 | Variable features |
| `--dims_use` | `1:10` | PCs for clustering |
| `--resolution` | 0.5 | Clustering resolution |
| `--umap_init` | `spectral` | UMAP initialization (`spectral`, `pca`, `spca`, `random`, ...) |
| `--threads` | 1 | Threads for Seurat's multithreaded steps |
| `--demo` | FALSE | Use demo dataset |

**Example:**
```bash
# Custom QC thresholds
Rscript scripts/01_basic_analysis.R \
  --input /path/to/filtered_feature_bc_matrix \
  --output ./my_analysis \
  --min_features 500 \
  --max_features 5000 \
  --max_mt 10 \
  --resolution 0.8 \
  --dims_use "1:15"
```

---

### 2. SCTransform (`02_sctransform.R`)

Improved normalization using regularized negative binomial regression.

**Key Parameters:**
| Parameter | Default | Description |
|-----------|---------|-------------|
| `--input` | required | Input file/directory |
| `--n_variable_features` | 3000 | Variable features |
| `--vars_to_regress` | NULL | Variables to regress (comma-separated) |
| `--vst_flavor` | `v2` | SCTransform version (v1, v2) |
| `--max_features` | none | Max features per cell; filters cells before SCTransform when set |
| `--max_mt` | none | Max % mitochondrial; filters cells before SCTransform when set |
| `--dims_use` | `1:30` | PCs for downstream |
| `--umap_init` | `spectral` | UMAP initialization |
| `--threads` | 1 | Threads for Seurat's multithreaded steps |

**Example:**
```bash
# With mitochondrial regression
Rscript scripts/02_sctransform.R \
  --input ./data/sample.h5 \
  --output ./sct_results \
  --vars_to_regress "percent.mt" \
  --vst_flavor v2
```

---

### 3. Integration (`03_integration.R`)

Integrate multiple datasets using Seurat v5 methods.

**Key Parameters:**
| Parameter | Default | Description |
|-----------|---------|-------------|
| `--input_list` | required | File with paths or comma-separated paths |
| `--method` | `CCAIntegration` | Integration method |
| `--sample_names` | auto | Sample names |
| `--conditions` | NULL | Condition labels |
| `--normalization` | `LogNormalize` | Normalization (LogNormalize, SCT) |
| `--dims_use` | `1:30` | PCs used for integration, UMAP, and clustering |
| `--reference` | NULL | Reference sample index, 1-based (CCA/RPCA only) |
| `--k_anchor` | 5 | Neighbors used to pick anchors (CCA/RPCA only) |
| `--svd_method` | `auto` | SVD backend for CCA: `irlba`, or `rspectra`, which avoids building the cell-by-cell matrix (`auto` picks `rspectra` for large inputs) |
| `--umap_init` | `spectral` | UMAP initialization |
| `--resolutions` | NULL | Extra clustering resolutions (comma-separated); `--resolution` sets `seurat_clusters` |

**Integration Methods:**
- `CCAIntegration` - Canonical correlation analysis
- `RPCAIntegration` - Reciprocal PCA (faster)
- `HarmonyIntegration` - Harmony
- `FastMNNIntegration` - Fast MNN (needs [SeuratWrappers](https://github.com/satijalab/seurat-wrappers) and Bioconductor `batchelor`, which are not installed by default)

**Example:**
```bash
# Create sample list file
echo "/path/to/sample1
/path/to/sample2
/path/to/sample3" > samples.txt

# Run integration
Rscript scripts/03_integration.R \
  --input_list samples.txt \
  --sample_names "Ctrl,Stim1,Stim2" \
  --conditions "control,treatment,treatment" \
  --method CCAIntegration \
  --output ./integrated

# Or with demo data
Rscript scripts/03_integration.R --demo --output ./demo_integrated
```

---

### 4. Differential Expression (`04_differential_expression.R`)

Flexible DE testing with multiple frameworks.

**Key Parameters:**
| Parameter | Default | Description |
|-----------|---------|-------------|
| `--input` | required | Seurat RDS file |
| `--mode` | `all_markers` | Analysis mode |
| `--test_use` | `wilcox` | Statistical test |
| `--min_pct` | 0.1 | Min fraction expressing |
| `--logfc_threshold` | 0.25 | Log FC threshold |
| `--pseudobulk` | FALSE | Pseudobulk DESeq2 (`between_conditions` mode only) |
| `--pseudobulk_group` | NULL | Sample/replicate column; required with `--pseudobulk`, needs at least 2 samples per condition |

With `--pseudobulk`, the script also draws `pseudobulk_heatmap.png`: the top up-
and down-regulated genes across the pseudobulk samples, ordered and annotated by
condition (Seurat 6 `PlotAnnotatedHeatmap`; needs Bioconductor `ComplexHeatmap`,
which `setup.R` installs).

**Modes:**
- `all_markers` - Find markers for all clusters
- `between_clusters` - Compare two clusters
- `between_conditions` - Compare conditions

**Tests:**
- `wilcox` - Wilcoxon rank sum (fast, default)
- `t` - Student's t-test
- `MAST` - Model-based
- `DESeq2` - For pseudobulk
- `negbinom` - Negative binomial

**Example:**
```bash
# All cluster markers
Rscript scripts/04_differential_expression.R \
  --input ./results/seurat_object.rds \
  --mode all_markers \
  --test_use wilcox \
  --output ./de_results

# Condition comparison within a cluster
Rscript scripts/04_differential_expression.R \
  --input ./integrated/seurat_integrated.rds \
  --mode between_conditions \
  --condition_col "condition" \
  --condition_1 "treatment" \
  --condition_2 "control" \
  --subset_cluster "5" \
  --output ./de_cluster5
```

---

### 5. Cell Cycle (`05_cell_cycle.R`)

Cell cycle scoring and regression.

**Key Parameters:**
| Parameter | Default | Description |
|-----------|---------|-------------|
| `--input` | required | Input file or Seurat RDS |
| `--species` | `human` | Species (human, mouse) |
| `--regress_cc` | FALSE | Regress out cell cycle |
| `--regress_cc_difference` | FALSE | Regress difference only |
| `--gene_set_version` | `2019` | Gene set version |
| `--umap_init` | `spectral` | UMAP initialization |
| `--threads` | 1 | Threads for Seurat's multithreaded steps |

**Example:**
```bash
# Score cell cycle
Rscript scripts/05_cell_cycle.R \
  --input ./results/seurat_object.rds \
  --species human \
  --output ./cell_cycle

# With regression
Rscript scripts/05_cell_cycle.R \
  --input ./data/sample.h5 \
  --regress_cc \
  --output ./cc_regressed
```

---

### 6. Visualization (`06_visualization.R`)

Generate comprehensive visualizations.

**Key Parameters:**
| Parameter | Default | Description |
|-----------|---------|-------------|
| `--input` | required | Seurat RDS file |
| `--features` | NULL | Genes to plot (comma-separated) |
| `--marker_file` | NULL | CSV with marker genes |
| `--top_markers` | 10 | Top markers per cluster taken from the marker file |
| `--max_plot_features` | 36 | Cap on marker-file genes plotted (`--features` genes always kept; 0 = no limit) |
| `--group_by` | `seurat_clusters` | Grouping variable |
| `--split_by` | NULL | Split variable |
| `--reduction` | `umap` | Reduction for plots |
| `--no_dim_plots`, `--no_feature_plots`, `--no_violin_plots`, `--no_dot_plot`, `--no_heatmap` | — | Skip a plot type (all on by default) |
| `--ridge_plots`, `--scatter_plots` | FALSE | Add ridge / feature scatter plots |
| `--no_label_clusters` | — | Don't label clusters on dim plots |

**Example:**
```bash
# Basic visualization
Rscript scripts/06_visualization.R \
  --input ./results/seurat_object.rds \
  --features "CD3D,CD4,CD8A,MS4A1,CD14,FCGR3A" \
  --output ./figures

# With split by condition
Rscript scripts/06_visualization.R \
  --input ./integrated/seurat_integrated.rds \
  --marker_file ./de_results/cluster_markers_top10.csv \
  --split_by "condition" \
  --output ./figures_split
```

---

### 7. Sample-Level Analysis (`07_sample_level.R`)

Seurat 6's [scSLIDE](https://satijalab.org/seurat/articles/scSLIDE_vignette)
workflow for studies with many samples (e.g. donors): each sample is summarized
by how its cells spread over a set of landmark cells, a diffusion map orders the
samples along the main axes of variation, and pseudobulk trajectory DE finds the
genes that change along each axis, per cell type.

Steps: `PrepareSampleObject` (balanced sketch, landmarks, response-oriented PLS
embedding, cell-to-landmark neighbors) → `GenerateSampleObject` → `RunDiffusionMap`
→ `TrajDETest` → `PlotAnnotatedHeatmap`. Instead of the vignette's cloud Azimuth
annotation (which uploads the data), the cell-type embedding is an existing
reduction, such as the `integrated` reduction from `03_integration.R`.

**Key Parameters:**
| Parameter | Default | Description |
|-----------|---------|-------------|
| `--input` | required | Seurat RDS with normalized data from several samples |
| `--sample_col` | required | Metadata column with sample IDs |
| `--response` | required | Sample-level variable the PLS embedding should capture (e.g. disease) |
| `--celltype_col` | NULL | Cell types: grouping for gene selection, landmark names, DE per cell type |
| `--reduction` / `--dims` | `pca` / `1:30` | Cell-type embedding combined with PLS |
| `--modality` | `multi` | `multi` (WNN of reduction + PLS) or `single` (PLS only) |
| `--pls_method` | `cppls` | `cppls`, `plsr`, or `spls` |
| `--ncells_per_sample` / `--ncells_landmark` | 1000 / 2000 | Training cells per sample / landmark cells |
| `--traj_dcs` | `1` | Diffusion components to test for DE |
| `--de_celltypes` | `all` | Cell types for trajectory DE |
| `--latent_vars` | NULL | Sample-level covariates (log library size is always included) |
| `--min_samples` | 6 | Minimum pseudobulk samples per cell type |

The script needs at least 4 samples; with fewer than ~10, only the first one or
two diffusion components carry information (the rest are skipped). The demo uses
ifnb with donor IDs from Seurat's DE tutorial (8 donors × 2 conditions = 16
samples; the ~5 MB donor files are downloaded once and cached).

**Outputs:** `sample_embedding.csv` (samples with their diffusion components),
`traj_de_<celltype>_DC<k>.csv`, `traj_de_summary.csv`, `sample_object.rds`,
diffusion map and heatmap figures.

**Example:**
```bash
Rscript scripts/07_sample_level.R \
  --input ./integrated/seurat_integrated.rds \
  --sample_col donor --response disease --celltype_col cell_type \
  --reduction integrated --latent_vars age,sex \
  --output ./sample_level

# Demo
Rscript scripts/07_sample_level.R --demo --output ./demo_sample_level
```

---

### 8. Imaging Spatial Data (`08_spatial.R`)

10x Genomics **Xenium** and **Atera** output folders, following the Xenium part
of Seurat's [spatial vignette](https://satijalab.org/seurat/articles/spatial_vignette_2)
and the [Atera tutorial](https://satijalab.org/seurat/articles/atera_tutorial):
QC, SCTransform, clustering, markers, clusters and genes on the tissue, transcript
positions, and an optional zoomed view with cell outlines. Atera's Azimuth cell
type annotations are plotted when present. The platform is detected from the folder.

**Key Parameters:**
| Parameter | Default | Description |
|-----------|---------|-------------|
| `--input` | required | Xenium or Atera output folder |
| `--platform` | `auto` | `auto`, `xenium`, or `atera` |
| `--molecules` | NULL | Genes whose transcript positions to load and plot (only these are read) |
| `--segmentations` | `auto` | Outlines to load: `cell`, `nucleus`, `none` (auto: cell for Xenium, none for Atera) |
| `--crop` | NULL | Zoomed view `xmin,xmax,ymin,ymax` in microns (see `cell_coordinates.csv`) |
| `--normalization` | `SCT` | `SCT` or `LogNormalize` |
| `--resolution` | 0.3 | Clustering resolution |
| `--features` | top markers | Genes to plot on the tissue |
| `--min_counts` | 1 | Minimum transcripts per cell |
| `--bpcells_dir` | NULL | Atera: keep counts on disk with BPCells |

Reading Parquet transcripts (newer Xenium outputs) needs `arrow`; reading Atera
needs `blosc`, `jsonlite` and `data.table`; `--crop` needs `sf`. `setup.R --all` and the conda
environments install them. On the 10x Xenium mouse brain "tiny subset" (36,602
cells, 248 genes) the script runs in about a minute and peaks at ~6 GB.

There is no demo dataset (the smallest public Xenium example is a 3.5 GB
download); the test suite runs the script on a small synthetic Xenium folder
written by `tests/make_xenium_fixture.R`.

**Example:**
```bash
Rscript scripts/08_spatial.R \
  --input ./Xenium_V1_FF_Mouse_Brain_Coronal_Subset_CTX_HP_outs \
  --molecules Gad1,Sst,Pvalb,Gfap \
  --crop 2000,2600,1000,1500 \
  --threads 4 --output ./xenium_results
```

---

## Input Formats

The scripts support multiple input formats:

| Format | Example | Description |
|--------|---------|-------------|
| `10x_mtx` | directory with matrix.mtx | 10x Genomics matrix |
| `10x_h5` | `sample.h5` | 10x HDF5 format |
| `rds` | `seurat.rds` | R data file |
| `csv` | `counts.csv` | CSV matrix |
| `seurat_data` | `pbmc3k` | SeuratData package |

Format is auto-detected, or specify with `--format`.

## Output Structure

Each script creates a structured output directory:

```
output/
├── seurat_object.rds      # Processed Seurat object
├── cluster_assignments.csv # Cell-cluster mapping
├── cluster_markers_all.csv # DE results
├── cluster_markers_top10.csv
├── variable_features.csv
├── umap_coordinates.csv
├── parameters.csv         # Run parameters
├── session_info.txt       # R session info
├── analysis.log           # Run log
└── *.png                  # Visualizations
```

## Example Workflows

### Complete Analysis Pipeline

```bash
# 1. Initial processing with SCTransform
Rscript scripts/02_sctransform.R \
  --input /path/to/data \
  --output ./01_sctransform \
  --vars_to_regress "percent.mt"

# 2. Cell cycle scoring
Rscript scripts/05_cell_cycle.R \
  --input ./01_sctransform/seurat_sctransform.rds \
  --regress_cc_difference \
  --output ./02_cell_cycle

# 3. Find markers
Rscript scripts/04_differential_expression.R \
  --input ./02_cell_cycle/seurat_cell_cycle.rds \
  --mode all_markers \
  --output ./03_markers

# 4. Generate figures
Rscript scripts/06_visualization.R \
  --input ./02_cell_cycle/seurat_cell_cycle.rds \
  --marker_file ./03_markers/cluster_markers_top10.csv \
  --output ./04_figures
```

### Multi-Sample Integration

```bash
# Prepare sample list
cat > samples.txt << EOF
/data/sample1/filtered_feature_bc_matrix
/data/sample2/filtered_feature_bc_matrix
/data/sample3/filtered_feature_bc_matrix
EOF

# Run integration
Rscript scripts/03_integration.R \
  --input_list samples.txt \
  --sample_names "Control,Treatment_A,Treatment_B" \
  --conditions "ctrl,treat,treat" \
  --method CCAIntegration \
  --normalization SCT \
  --output ./integrated_analysis

# Differential expression between conditions
Rscript scripts/04_differential_expression.R \
  --input ./integrated_analysis/seurat_integrated.rds \
  --mode between_conditions \
  --condition_col condition \
  --condition_1 treat \
  --condition_2 ctrl \
  --pseudobulk \
  --pseudobulk_group orig.ident \
  --output ./de_pseudobulk
```

## Workflow Pipelines

For automated, reproducible analyses, use the included **Snakemake** or **Nextflow** pipelines.

**Prerequisites:** Install [Snakemake](https://snakemake.readthedocs.io/) (`pip install snakemake`) or [Nextflow](https://www.nextflow.io/) (`curl -s https://get.nextflow.io | bash`).

### Available Workflows

| Workflow | Pipeline | Description |
|----------|----------|-------------|
| `basic` | 01 → 04 → 06 | Basic analysis with log normalization |
| `sctransform` | 02 → 04 → 06 | SCTransform normalization |
| `integration` | 03 → 04 → 06 | Multi-sample integration |
| `cell_cycle` | 02 → 05 → 04 → 06 | SCTransform with cell cycle scoring |
| `full` | 02 → 05 → 04 → 06 | Complete analysis with cell cycle |

### Snakemake

```bash
cd workflows/snakemake

# Edit config.yaml with your input paths and parameters
# Then run:
snakemake --configfile config.yaml --cores 4

# Dry run to see what will execute
snakemake --configfile config.yaml --cores 4 -n
```

### Nextflow

```bash
cd workflows/nextflow

# Run with your data
nextflow run main.nf --input /path/to/data --outdir results

# Run with demo data
nextflow run main.nf --demo --outdir demo_results

# Run specific workflow
nextflow run main.nf --input /path/to/data --workflow full

# Run with Docker
nextflow run main.nf --input /path/to/data -profile docker
```

See [workflows/README.md](workflows/README.md) for detailed documentation.

## Testing

Verify your installation by running the demo tests:

```bash
# Quick test (scripts 1-2 only, ~5 min)
Rscript tests/run_demo_tests.R --quick

# Full test suite (scripts 1-8, ~5 min)
Rscript tests/run_demo_tests.R

# Test specific script (1-8; script 8 runs on a generated synthetic Xenium folder)
Rscript tests/run_demo_tests.R --script 1

# Validation options
Rscript tests/run_demo_tests.R --skip_validation    # Skip output validation
Rscript tests/run_demo_tests.R --validation_only    # Only validate existing output

# Or using make
make test-quick
make test

# Test workflow pipelines (requires Snakemake/Nextflow)
./tests/test_pipelines.sh           # Test both
./tests/test_pipelines.sh snakemake # Snakemake only
./tests/test_pipelines.sh nextflow  # Nextflow only
```

### Validation Framework

The test runner includes a two-tier validation system:

- **Tier 1 (Structural)**: Validates output files exist, CSV structure is correct, cluster assignments are valid integers, Seurat objects load properly, and every PNG figure has content (not a blank image). Runs on all data.
- **Tier 2 (Biological)**: Validates demo-specific expectations like expected cluster counts, cell counts, and known marker genes. Only runs on demo datasets.

## Installation Requirements

Use the setup script (recommended):
```bash
Rscript setup.R            # Core + recommended packages + demo datasets
Rscript setup.R --minimal  # Core packages only
Rscript setup.R --bioc     # Core + recommended + Bioconductor DE packages
Rscript setup.R --all      # Everything including optional integration packages
```

### Package Overview

| Category | Packages | Purpose |
|----------|----------|---------|
| **Core** | Seurat, dplyr, ggplot2, patchwork, optparse, Matrix, future, viridis | Required for all scripts |
| **Recommended** | SeuratData, presto, hdf5r, sctransform, glmGamPoi | Demo data, faster markers, HDF5 support, SCTransform v2 |
| **Bioconductor DE** | MAST, DESeq2, SingleCellExperiment | Alternative DE methods |
| **Seurat 6 features** | ComplexHeatmap, destiny | Annotated heatmaps (04, 07), diffusion maps (07) |
| **Optional** | harmony, batchelor, scran, scater | Integration methods |
| **Optional (`--all`)** | arrow, blosc, jsonlite, data.table, sf, spls | Xenium Parquet / Atera zarr reading and zoomed views (08), sparse PLS (07) |

### Demo Datasets

The setup script also installs demo datasets:
- **pbmc3k** (~8MB) - Used by scripts 01, 02, 04, 05, 06
- **ifnb** (~400MB) - Used by scripts 03 (integration) and 07 (sample-level; plus ~5 MB of donor IDs downloaded on first use and cached in the user cache directory, or `$SEURAT_CLI_CACHE`)

### Manual Installation

```r
# Core packages
install.packages(c("Seurat", "dplyr", "ggplot2", "patchwork",
                   "optparse", "Matrix", "future", "viridis"))

# Recommended (from GitHub)
install.packages("remotes")
remotes::install_github("satijalab/seurat-data")
remotes::install_github("immunogenomics/presto")

# Demo datasets (after SeuratData is installed)
SeuratData::InstallData("pbmc3k")
SeuratData::InstallData("ifnb")

# Optional Bioconductor packages
BiocManager::install(c("MAST", "DESeq2", "glmGamPoi"))

# Optional integration packages
install.packages("harmony")
BiocManager::install(c("batchelor", "scran", "scater"))
```

## Tips

1. **Demo Mode**: Use `--demo` flag to test scripts with built-in datasets
2. **Parallelization**: Use `--threads N` to run Seurat's multithreaded steps (normalization, scaling, SCTransform, PCA, neighbors, clustering, UMAP, CCA/RPCA integration) on N threads. The threads share one copy of the data, so memory use barely changes with N. Gains depend on dataset size: integration of ifnb (14k cells) runs ~16% faster with 4 threads; pbmc3k is too small to benefit much
3. **Reproducibility**: Results include `parameters.csv` and `session_info.txt`
4. **Custom Features**: Most scripts accept `--features` for specific genes
5. **Species**: Set `--species mouse` for mouse data (adjusts MT gene patterns)

## Project Structure

```
seurat-cli/
├── scripts/                    # Main analysis scripts
│   ├── 01_basic_analysis.R
│   ├── 02_sctransform.R
│   ├── 03_integration.R
│   ├── 04_differential_expression.R
│   ├── 05_cell_cycle.R
│   ├── 06_visualization.R
│   ├── 07_sample_level.R       # scSLIDE (Seurat 6)
│   └── 08_spatial.R            # Xenium / Atera (Seurat 6)
├── workflows/                  # Pipeline workflows
│   ├── snakemake/             # Snakemake pipeline
│   │   ├── Snakefile
│   │   └── config.yaml
│   ├── nextflow/              # Nextflow pipeline
│   │   ├── main.nf
│   │   └── nextflow.config
│   └── README.md
├── utils/
│   └── common.R               # Shared utility functions
├── tests/
│   ├── run_demo_tests.R       # R script test runner with validation
│   ├── validate_outputs.R     # Two-tier validation framework
│   ├── make_xenium_fixture.R  # Synthetic Xenium data for testing 08
│   └── test_pipelines.sh      # Pipeline test runner
├── config/                    # Configuration templates
├── output/                    # Default output directory
├── Dockerfile                 # Docker image definition
├── environment.yml            # Conda/mamba environment
├── setup.R                    # Dependency installer
├── install_seurat.R           # Installs Seurat 6 from CRAN when conda's is older
├── Makefile                   # Common commands
└── seurat-cli.Rproj           # RStudio project file
```

## Enhancements Beyond Vignettes

While these scripts faithfully implement the Seurat vignettes, they include several enhancements for robust CLI usage:

### Seurat v5 Compatibility

- **Automatic object updates**: Uses `UpdateSeuratObject()` for compatibility with older Seurat objects
- **Layer-based integration**: Script 03 properly splits RNA assay into layers before `IntegrateLayers()` as required by Seurat v5

### Error Handling

- **Empty marker results**: All scripts gracefully handle cases where `FindAllMarkers()` returns no results (e.g., due to data characteristics or stringent filtering). Instead of crashing, scripts log a warning and create empty output files.
- **Empty DE results**: Script 04 provides informative messages when no differentially expressed genes are found, suggesting potential causes (filtering parameters, cell counts, preprocessing).
- **Cell cycle gene validation**: Script 05 warns when few cell cycle genes are found in the dataset (may indicate species mismatch) and fails gracefully with helpful messages if no genes are found.
- **Metadata validation**: Script 06 validates that grouping and splitting columns exist in metadata before attempting to use them, with helpful warnings listing available columns.
- **Reduction availability**: Script 06 checks for available dimensional reductions and skips relevant plots with informative messages if none exist.

### Reproducibility

- **Parameter logging**: All scripts save `parameters.csv` with the exact parameters used for each run
- **Session info**: All scripts save `session_info.txt` with R version and package versions
- **Comprehensive logging**: Each script creates an `analysis.log` file tracking execution progress

### Execution Flexibility

- **Dual execution support**: Scripts work both via `Rscript` command line and `source()` in RStudio
- **Demo mode**: Scripts 01-07 support `--demo` for testing with built-in datasets (script 08 has no small public dataset; the tests generate a synthetic one)
- **Multi-threading**: `--threads` sets Seurat's native multithreading (`setThreads()`, Seurat >= 5.6), so threads share one copy of the data

## References

- [Seurat Documentation](https://satijalab.org/seurat/)
- [PBMC 3K Tutorial](https://satijalab.org/seurat/articles/pbmc3k_tutorial)
- [SCTransform Paper](https://genomebiology.biomedcentral.com/articles/10.1186/s13059-019-1874-1)
- [Integration Methods](https://satijalab.org/seurat/articles/integration_introduction)

## Author

Created and maintained by [@bminie](https://github.com/bminie).

## Contributing

Contributions are welcome! Please see [CONTRIBUTING.md](CONTRIBUTING.md) for guidelines.

## License

MIT License - see [LICENSE](LICENSE) for details.

Scripts are provided as educational implementations of published methods.
