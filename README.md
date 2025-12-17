# Seurat CLI

[![R Syntax Check](https://github.com/bminie/seurat-cli/actions/workflows/r-check.yml/badge.svg)](https://github.com/bminie/seurat-cli/actions/workflows/r-check.yml)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)

A collection of flexible, parameterized R scripts that implement the major [Seurat](https://satijalab.org/seurat/) vignettes for single-cell RNA-seq analysis.

## Overview

These scripts provide command-line interfaces to common Seurat workflows, allowing users to:
- Run standard analyses with custom parameters
- Process different input formats (10x, CSV, RDS)
- Automate batch processing
- Generate reproducible results with logged parameters

## Installation

```bash
# Clone the repository
git clone https://github.com/bminie/seurat-cli.git
cd seurat-cli

# Install dependencies
Rscript setup.R           # Core packages
Rscript setup.R --all     # All packages including optional
```

Or open in RStudio by double-clicking `seurat-cli.Rproj`.

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
| `--dims_use` | `1:30` | PCs for downstream |

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

**Integration Methods:**
- `CCAIntegration` - Canonical correlation analysis
- `RPCAIntegration` - Reciprocal PCA (faster)
- `HarmonyIntegration` - Harmony
- `FastMNNIntegration` - Fast MNN

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
| `--pseudobulk` | FALSE | Use pseudobulk DE |

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
| `--group_by` | `seurat_clusters` | Grouping variable |
| `--split_by` | NULL | Split variable |
| `--reduction` | `umap` | Reduction for plots |

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

## Testing

Verify your installation by running the demo tests:

```bash
# Quick test (scripts 1-2 only, ~5 min)
Rscript tests/run_demo_tests.R --quick

# Full test suite (~15-30 min)
Rscript tests/run_demo_tests.R

# Test specific script
Rscript tests/run_demo_tests.R --script 1

# Or using make
make test-quick
make test
```

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
| **Optional** | harmony, batchelor, scran, scater | Integration methods |

### Demo Datasets

The setup script also installs demo datasets:
- **pbmc3k** (~8MB) - Used by scripts 01, 02, 04, 05, 06
- **ifnb** (~400MB) - Used by script 03 (integration)

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
2. **Parallelization**: Use `--threads N` for multi-threaded operations
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
│   └── 06_visualization.R
├── utils/
│   └── common.R               # Shared utility functions
├── tests/
│   └── run_demo_tests.R       # Test runner
├── config/                    # Configuration templates
├── output/                    # Default output directory
├── setup.R                    # Dependency installer
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
- **Demo mode**: All scripts support `--demo` flag for testing with built-in datasets
- **Multi-threading**: Scripts support `--threads` parameter for parallel execution

## References

- [Seurat Documentation](https://satijalab.org/seurat/)
- [PBMC 3K Tutorial](https://satijalab.org/seurat/articles/pbmc3k_tutorial)
- [SCTransform Paper](https://genomebiology.biomedcentral.com/articles/10.1186/s13059-019-1874-1)
- [Integration Methods](https://satijalab.org/seurat/articles/integration_introduction)

## Contributing

Contributions are welcome! Please see [CONTRIBUTING.md](CONTRIBUTING.md) for guidelines.

## License

MIT License - see [LICENSE](LICENSE) for details.

Scripts are provided as educational implementations of published methods.
