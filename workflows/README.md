# Seurat CLI Workflow Pipelines

This directory contains workflow pipelines for running Seurat CLI scripts in sequence using either **Snakemake** or **Nextflow**.

## Prerequisites

### Option 1: Mamba (Recommended)

**Requires:** [Miniforge](https://github.com/conda-forge/miniforge#miniforge3) (includes mamba) or [Mamba](https://mamba.readthedocs.io/)

> **Why Mamba?** This environment has complex dependencies (R, Seurat, Snakemake, Nextflow, Bioconductor packages). Mamba's C++ solver creates the environment in 2-5 minutes, while conda can take 30+ minutes or fail.

**Install Miniforge (includes mamba):**
- **macOS/Linux:**
  ```bash
  curl -L -O "https://github.com/conda-forge/miniforge/releases/latest/download/Miniforge3-$(uname)-$(uname -m).sh"
  bash Miniforge3-$(uname)-$(uname -m).sh
  ```
- **Windows:** Download from https://github.com/conda-forge/miniforge#miniforge3

**Or add mamba to existing conda:**
```bash
conda install -n base -c conda-forge mamba
```

**Then create the environment:**

```bash
# From the project root (recommended)
mamba env create -f environment.yml
conda activate seurat-cli
```

This installs R, Seurat, Snakemake, Nextflow, and all dependencies in one command.

### Option 1b: Conda (Slower Alternative)

If you prefer conda, first enable the faster libmamba solver:

```bash
conda install -n base -c conda-forge conda-libmamba-solver
conda config --set solver libmamba
conda env create -f environment.yml
conda activate seurat-cli
```

> **Warning:** Without libmamba solver, conda may take 30+ minutes or fail to solve the environment.

### Option 2: Manual Installation

If you prefer manual installation, you'll need to install each component separately.

**Requires:**
- [R](https://cran.r-project.org/) version 4.2 or higher
- [Python](https://www.python.org/) 3.8+ (for Snakemake)
- [Java](https://adoptium.net/) 11+ (for Nextflow)

**R Dependencies:**
```bash
Rscript setup.R        # Core packages
Rscript setup.R --all  # All packages including optional
```

**Snakemake** (Python-based):
```bash
pip install snakemake
```

**Nextflow** (Java-based):
```bash
# Check Java version (needs 11+)
java -version

# Install Nextflow
curl -s https://get.nextflow.io | bash
sudo mv nextflow /usr/local/bin/
```

## Available Workflows

Both pipelines support the following workflow types:

| Workflow | Pipeline | Description |
|----------|----------|-------------|
| `basic` | 01 → 04 → 06 | Basic analysis with log normalization |
| `sctransform` | 02 → 04 → 06 | SCTransform normalization |
| `integration` | 03 → 04 → 06 | Multi-sample integration |
| `cell_cycle` | 02 → 05 → 04 → 06 | SCTransform with cell cycle scoring |
| `full` | 02 → 05 → 04 → 06 | Complete analysis with cell cycle |

**Script Legend:**
- 01: Basic Analysis
- 02: SCTransform
- 03: Integration
- 04: Differential Expression
- 05: Cell Cycle
- 06: Visualization

## Snakemake

### Quick Start

```bash
cd workflows/snakemake

# Dry run to see what will be executed
snakemake --configfile config.yaml -n

# Run with 4 cores
snakemake --configfile config.yaml --cores 4

# Run demo
snakemake --configfile config_demo.yaml --cores 4
```

### Configuration

Edit `config.yaml` to set:
- `workflow`: Which workflow to run
- `input_path`: Path to your data
- `output_dir`: Where to save results
- Script-specific parameters

### Running Partial Pipelines

```bash
# Run only up to differential expression
snakemake --configfile config.yaml --cores 4 --until differential_expression

# Run only basic analysis
snakemake --configfile config.yaml --cores 4 basic_analysis_only
```

## Nextflow

### Quick Start

```bash
cd workflows/nextflow

# Show help
nextflow run main.nf --help

# Run with demo data
nextflow run main.nf --demo --outdir demo_results

# Run with your data
nextflow run main.nf --input /path/to/data --outdir results

# Run specific workflow
nextflow run main.nf --input /path/to/data --workflow sctransform
```

### Profiles

```bash
# Run with Docker
nextflow run main.nf --demo -profile docker

# Run on SLURM cluster
nextflow run main.nf --input /path/to/data -profile slurm

# Run demo test
nextflow run main.nf -profile test
```

### Configuration

Edit `nextflow.config` or pass parameters on command line:

```bash
nextflow run main.nf \
    --input /path/to/data \
    --workflow full \
    --cc_regress_difference \
    --sct_vars_to_regress "percent.mt" \
    --outdir my_results
```

## Input Requirements

### Single Sample (basic, sctransform, cell_cycle, full)

- 10x Genomics directory (with `matrix.mtx`, `genes.tsv`, `barcodes.tsv`)
- 10x HDF5 file (`.h5`)
- Seurat RDS file (`.rds`)
- CSV count matrix

### Multi-Sample (integration)

Create a text file with paths to each sample (one per line):

```
# samples.txt
/path/to/sample1
/path/to/sample2
/path/to/sample3
```

Then use:
- Snakemake: Set `input_list: "samples.txt"` in config
- Nextflow: `--input_list samples.txt --workflow integration`

## Output Structure

Both pipelines create similar output structures:

```
results/
├── 01_sctransform/          # Or 01_basic, 01_integration
│   ├── seurat_*.rds
│   ├── cluster_markers_*.csv
│   └── *.png
├── 02_cell_cycle/           # If workflow includes cell cycle
│   ├── seurat_cell_cycle.rds
│   └── cell_cycle_scores.csv
├── 03_differential_expression/
│   ├── de_results_all.csv
│   └── de_results_significant.csv
├── visualization/
│   ├── visualization_summary.md
│   └── *.png
├── logs/                    # Snakemake only
└── pipeline_info/           # Nextflow only
    ├── timeline.html
    ├── report.html
    └── dag.svg
```

## Testing

```bash
# Test Snakemake pipeline
cd workflows/snakemake
snakemake --configfile config_demo.yaml --cores 2 -n  # Dry run
snakemake --configfile config_demo.yaml --cores 2     # Actual run

# Test Nextflow pipeline
cd workflows/nextflow
nextflow run main.nf -profile test
```
