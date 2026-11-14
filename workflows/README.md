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
The scripts need Seurat >= 6.0, which conda-forge may not have yet; then run
`Rscript setup.R` in the activated environment (it builds Seurat 6 from CRAN via
`install_seurat.R` and installs the demo datasets). `tests/test_pipelines.sh`
does this automatically.

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
| `sample_level` | 03 → 07 | scSLIDE sample-level analysis of integrated samples (or 07 alone on an existing RDS) |
| `spatial` | 08 | Xenium/Atera imaging data (input is the output folder; no demo) |

**Script Legend:**
- 01: Basic Analysis
- 02: SCTransform
- 03: Integration
- 04: Differential Expression
- 05: Cell Cycle
- 06: Visualization
- 07: Sample-Level Analysis (scSLIDE)
- 08: Imaging Spatial (Xenium/Atera)

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

# Run with Apptainer or Singularity
nextflow run main.nf --demo -profile apptainer

# Run on a SLURM cluster with Apptainer (see "Running on an HPC cluster")
nextflow run main.nf --input /path/to/data -profile slurm,apptainer

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

Create a text file with paths to each sample (one per line; blank lines and
lines starting with `#` are skipped). Use absolute paths, or paths relative to
the directory you launch the pipeline from:

```
# samples.txt
/path/to/sample1
/path/to/sample2
/path/to/sample3
```

Nextflow stages each listed sample into the integration task (as
`sample1/`, `sample2/`, ...), so samples with the same file name are fine and
containers can read them. Snakemake checks that every listed sample exists
before it starts.

Then use:
- Snakemake: Set `input_list: "samples.txt"` in config
- Nextflow: `--input_list samples.txt --workflow integration`

### Sample-Level Analysis (sample_level)

Runs integration and then `07_sample_level.R` on the integrated object. Each
sample in the list is one sample (e.g. donor) for scSLIDE, and the integration
conditions are the response, so set them (`integration.conditions` /
`--int_conditions`). At least 4 samples are needed, ideally 10 or more. By
default cell types are the integration clusters; change `sample_level.celltype_col`
(`--sl_celltype_col`) if the object has annotations.

To analyze an existing Seurat object instead, set `sample_level.input_rds`
(`--sl_input_rds`) and the column names (`sample_col`, `response`, ...). The demo
(`config_demo.yaml` with `workflow=sample_level`, or `-profile demo --workflow
sample_level`) runs script 07's own demo: ifnb with donor IDs, 16 samples.

### Imaging Spatial Data (spatial)

Set the input to a 10x Xenium or Atera output folder (`input_path` /
`--input`). Options are under `spatial` (`--sp_*`): `molecules`, `crop`,
`normalization`, `resolution`, `platform`. There is no demo dataset; the pipeline
tests use the synthetic folder from `tests/make_xenium_fixture.R`.

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
├── sample_level/            # sample_level workflow
│   ├── traj_de_summary.csv
│   ├── sample_embedding.csv
│   └── *.png
├── spatial/                 # spatial workflow
│   ├── seurat_spatial.rds
│   ├── cluster_assignments.csv
│   └── *.png
├── logs/                    # Snakemake only
└── pipeline_info/           # Nextflow only
    ├── timeline.html
    ├── report.html
    └── dag.svg
```

## Container Support

Both pipelines can run every step in the `ghcr.io/bminie/seurat-cli` image
(linux/amd64), which has all R packages, the scripts (`/app/scripts`), and the
demo datasets.

### Nextflow with Containers

```bash
# Docker (local development)
nextflow run main.nf --demo -profile docker

# Apptainer or Singularity (HPC)
nextflow run main.nf --demo -profile apptainer
nextflow run main.nf --demo -profile singularity

# A different image, or a local .sif file
nextflow run main.nf --demo -profile apptainer --container /shared/seurat-cli.sif
```

### Snakemake with Containers

```bash
# Apptainer/Singularity (Snakemake >= 8). Rules use the scripts in the image.
snakemake --configfile config.yaml --cores 4 --sdm apptainer

# Only the working directory is mounted in the container; bind data elsewhere
snakemake --configfile config.yaml --cores 4 --sdm apptainer \
    --apptainer-args "--bind /data"
```

Set `container` in the config to use a different image or a local `.sif` file.

### Building the Container Locally

```bash
# From the project root
docker build -t seurat-cli .

# Use it with the Nextflow docker profile
nextflow run main.nf --demo -profile docker --container seurat-cli

# Convert to a .sif file for Apptainer/Singularity
apptainer build seurat-cli.sif docker-daemon://seurat-cli:latest
```

## Running on an HPC Cluster

Run the pipeline from a login node (in `tmux`/`screen`, or as a long
low-resource job). It submits each step as its own cluster job, and each job
runs in the container through Apptainer or Singularity. Each step asks for its
own memory and time: 16 GB / 4 h for integration, sample-level, and spatial,
6-8 GB for the rest, and `threads` CPUs.

### Get the image onto the cluster

The pipelines pull the image on first use. Pull it once to shared storage
instead, which is faster and works when compute nodes have no internet access:

```bash
apptainer pull /shared/containers/seurat-cli.sif docker://ghcr.io/bminie/seurat-cli:latest
```

Then pass `--container /shared/containers/seurat-cli.sif` (Nextflow) or set
`container: /shared/containers/seurat-cli.sif` in the Snakemake config. To let
Nextflow pull it, set `NXF_APPTAINER_CACHEDIR` (or `NXF_SINGULARITY_CACHEDIR`)
to shared storage so runs reuse one copy.

### Nextflow (SLURM or PBS)

```bash
cd workflows/nextflow
nextflow run main.nf -profile slurm,apptainer \
    --input_list samples.txt --workflow integration \
    --container /shared/containers/seurat-cli.sif \
    --queue short --cluster_options="--account=mylab"
```

- `slurm` or `pbs`, plus `apptainer` or `singularity`
- `--queue`: partition/queue (default: the scheduler's default)
- `--cluster_options`: extra submit options. Write it with `=`; otherwise
  Nextflow reads a value starting with `--` as a separate flag.
- `--max_memory`, `--max_cpus`, `--max_time`: caps for any one job (default
  16 GB, 8, 24 h). Raise them for large datasets, e.g. `--max_memory 64.GB`.
- `--queue_size`: most jobs queued or running at once (default 50)

Nextflow mounts input files into the container itself.

### Snakemake (SLURM)

Needs `snakemake-executor-plugin-slurm` (in `environment.yml`):

```bash
cd workflows/snakemake
snakemake --configfile my_config.yaml --workflow-profile profiles/slurm \
    --apptainer-args "--bind /data,/scratch"
```

`profiles/slurm/config.yaml` runs every rule as a SLURM job in the container;
copy it to set your partition, account, binds, and image cache folder. Memory
and time per rule come from `resources` in the config:

```yaml
container: /shared/containers/seurat-cli.sif
threads: 8
resources:
  default: {mem_mb: 16000, runtime: 240}
  integration: {mem_mb: 64000, runtime: 480}
```

Bind every folder outside the working directory that holds input data or
paths from the sample list (`--apptainer-args "--bind ..."`). If `sacct`
reports that accounting is disabled on your cluster, add
`--slurm-status-command squeue`.

For other schedulers (LSF, SGE, ...), set the matching Nextflow executor in a
config file passed with `-c` (`process.executor = 'lsf'`), or use the matching
Snakemake executor plugin.

## Testing

```bash
# Test Snakemake pipeline
cd workflows/snakemake
snakemake --configfile config_demo.yaml --cores 2 -n  # Dry run
snakemake --configfile config_demo.yaml --cores 2     # Actual run

# Test Nextflow pipeline
cd workflows/nextflow
nextflow run main.nf -profile test

# Test with containers
nextflow run main.nf -profile test,docker
```
