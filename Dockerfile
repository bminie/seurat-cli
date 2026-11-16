# =============================================================================
# Seurat CLI Docker Image
# =============================================================================
#
# Build:
#   docker build -t seurat-cli .
#
# Run:
#   docker run --rm -v $(pwd)/data:/data -v $(pwd)/output:/output seurat-cli \
#     Rscript scripts/01_basic_analysis.R --input /data/sample --output /output
#
# Run with demo:
#   docker run --rm -v $(pwd)/output:/output seurat-cli \
#     Rscript scripts/01_basic_analysis.R --demo --output /output
#
# Interactive:
#   docker run --rm -it -v $(pwd):/work seurat-cli bash
#
# Published image (linux/amd64, built by .github/workflows/docker-publish.yml):
#   docker pull ghcr.io/bminie/seurat-cli:latest
#   docker pull bminie7/seurat-cli:latest           (Docker Hub)
#   apptainer pull seurat-cli.sif docker://ghcr.io/bminie/seurat-cli:latest
#
# =============================================================================

# rocker installs R packages from Posit Package Manager: the newest R release
# follows the "latest" snapshot, and each older release is pinned to the
# snapshot from its last day as the current release. Seurat 6 needs a snapshot
# from November 2026 or later, so this must be R 4.6.2 or newer. The installs
# below use the image's default repos (no repos= override) so packages stay
# compatible with R.
FROM rocker/r-ver:4.6.2

LABEL maintainer="bminie"
LABEL description="Seurat CLI - Command-line tools for single-cell RNA-seq and imaging spatial (Xenium/Atera) analysis"
LABEL org.opencontainers.image.source="https://github.com/bminie/seurat-cli"

# Set environment variables
ENV DEBIAN_FRONTEND=noninteractive
ENV TZ=UTC

# Install system dependencies
RUN apt-get update && apt-get install -y --no-install-recommends \
    # Build tools
    build-essential \
    libcurl4-openssl-dev \
    libssl-dev \
    libxml2-dev \
    libfontconfig1-dev \
    libfreetype-dev \
    libpng-dev \
    libtiff-dev \
    libjpeg-dev \
    # HDF5 for hdf5r
    libhdf5-dev \
    # GLPK for igraph (used by Seurat)
    libglpk-dev \
    # FFTW for sctransform
    libfftw3-dev \
    # GSL for some Bioconductor packages
    libgsl-dev \
    # GDAL for spatial packages
    libgdal-dev \
    libgeos-dev \
    libproj-dev \
    # udunits2 for the units package (needed by sf)
    libudunits2-dev \
    # blosc for reading Atera zarr archives (R blosc package)
    libblosc-dev \
    pkg-config \
    # Other utilities
    git \
    procps \
    && apt-get clean \
    && rm -rf /var/lib/apt/lists/*

# Install packages in parallel
RUN echo 'options(Ncpus = max(1L, parallel::detectCores()))' >> "${R_HOME}/etc/Rprofile.site"

# Install core R packages (CRAN snapshot from the base image)
RUN R -e "install.packages(c( \
    'remotes', \
    'BiocManager', \
    'optparse', \
    'dplyr', \
    'ggplot2', \
    'patchwork', \
    'Matrix', \
    'future', \
    'viridis', \
    'hdf5r', \
    'R.utils', \
    'png' \
))"

# Install Seurat and sctransform
RUN R -e "install.packages(c('Seurat', 'sctransform'))"

# Install glmGamPoi for faster SCTransform v2
RUN R -e "BiocManager::install('glmGamPoi', ask=FALSE, update=FALSE)"

# Install SeuratData and demo datasets
RUN R -e "remotes::install_github('satijalab/seurat-data', quiet=TRUE, upgrade='never')" \
    && R -e "SeuratData::InstallData('pbmc3k')" \
    && R -e "SeuratData::InstallData('ifnb')"

# Install presto for fast marker detection
RUN R -e "remotes::install_github('immunogenomics/presto', quiet=TRUE, upgrade='never')"

# Install Bioconductor packages for DE analysis
RUN R -e "BiocManager::install(c('MAST', 'DESeq2', 'SingleCellExperiment'), ask=FALSE, update=FALSE)"

# Install optional integration packages
RUN R -e "install.packages('harmony')" \
    && R -e "BiocManager::install(c('batchelor', 'scran', 'scater'), ask=FALSE, update=FALSE)"

# Seurat 6 features: annotated heatmaps (04, 07), diffusion maps (07),
# sparse PLS (07), Xenium Parquet / Atera zarr readers and zoomed views (08)
RUN R -e "BiocManager::install(c('ComplexHeatmap', 'destiny'), ask=FALSE, update=FALSE)" \
    && R -e "install.packages(c('arrow', 'blosc', 'jsonlite', 'data.table', 'sf', 'spls'))"

# Fail the build if any package is missing. install.packages() only warns when
# a package fails to install, so without this check the image could build
# without Seurat.
RUN Rscript -e " \
    pkgs <- c('optparse', 'dplyr', 'ggplot2', 'patchwork', 'Matrix', 'future', \
              'viridis', 'hdf5r', 'png', 'Seurat', 'SeuratObject', 'sctransform', \
              'glmGamPoi', 'SeuratData', 'pbmc3k.SeuratData', 'ifnb.SeuratData', \
              'presto', 'MAST', 'DESeq2', 'SingleCellExperiment', 'harmony', \
              'batchelor', 'scran', 'scater', 'ComplexHeatmap', 'destiny', \
              'arrow', 'blosc', 'jsonlite', 'data.table', 'sf', 'spls'); \
    missing <- pkgs[!vapply(pkgs, requireNamespace, logical(1), quietly = TRUE)]; \
    if (length(missing) > 0) { \
      message('Missing R packages: ', paste(missing, collapse = ', ')); quit(status = 1) \
    }; \
    if (packageVersion('Seurat') < '6.0.0') { \
      message('Seurat >= 6.0.0 required, got ', packageVersion('Seurat')); quit(status = 1) \
    }; \
    cat('All', length(pkgs), 'R packages available. Seurat', \
        as.character(packageVersion('Seurat')), '\n')"

# -----------------------------------------------------------------------------
# Cloud support
# -----------------------------------------------------------------------------

# AWS CLI v2: Nextflow's AWS Batch tasks copy their files to and from S3 with
# it. The download is checked against the AWS CLI team's signing key
# (docker/aws-cli.asc, from the AWS CLI install guide; fingerprint
# FB5D B77F D5C1 18B8 0511 ADA8 A631 0ACC 4672 475C). The key expires on
# 2027-07-01: replace the file with AWS's renewed key then.
ARG AWSCLI_VERSION=2.40.3
COPY docker/aws-cli.asc /tmp/aws-cli.asc
RUN apt-get update && apt-get install -y --no-install-recommends curl unzip gnupg \
    && rm -rf /var/lib/apt/lists/* \
    && cd /tmp \
    && url="https://awscli.amazonaws.com/awscli-exe-linux-$(uname -m)-${AWSCLI_VERSION}.zip" \
    && curl -fsSL -o awscliv2.zip "$url" \
    && curl -fsSL -o awscliv2.zip.sig "$url.sig" \
    && export GNUPGHOME="$(mktemp -d)" \
    && gpg --batch --quiet --import /tmp/aws-cli.asc \
    && gpg --batch --verify awscliv2.zip.sig awscliv2.zip \
    && unzip -q awscliv2.zip \
    && ./aws/install \
    && rm -rf /tmp/aws /tmp/awscliv2.zip /tmp/awscliv2.zip.sig /tmp/aws-cli.asc "$GNUPGHOME" \
    && aws --version

# Snakemake and its storage plugins. Snakemake's cloud executors (AWS Batch,
# Google Batch, Azure Batch) run Snakemake inside this image for each job, so
# these versions must match the Snakemake that launches the workflow
# (workflows/snakemake/environment-cloud.yml pins the same ones). The
# environment goes first on PATH: the executors run "python -m snakemake" and
# "pip install", and Ubuntu's own Python doesn't allow pip installs.
# python3-dev: some dependencies (immutables) have no prebuilt package for
# this Python version, so pip compiles them.
ARG SNAKEMAKE_VERSION=9.27.1
RUN apt-get update && apt-get install -y --no-install-recommends python3 python3-venv python3-dev \
    && rm -rf /var/lib/apt/lists/* \
    && python3 -m venv /opt/snakemake \
    && /opt/snakemake/bin/pip install --no-cache-dir \
        "snakemake==${SNAKEMAKE_VERSION}" \
        "snakemake-storage-plugin-s3==0.3.6" \
        "snakemake-storage-plugin-gcs==1.1.4" \
        "snakemake-storage-plugin-azure==0.4.4" \
    && /opt/snakemake/bin/snakemake --version
ENV PATH="/opt/snakemake/bin:${PATH}"

# Tells the Snakefile where the scripts are when Snakemake runs in the image
ENV SEURAT_CLI_SCRIPTS=/app/scripts

# Create working directory
WORKDIR /app

# Copy project files
COPY scripts/ /app/scripts/
COPY utils/ /app/utils/
COPY tests/ /app/tests/

# Make scripts accessible from PATH
ENV PATH="/app/scripts:${PATH}"

# Files the scripts would otherwise download on first use (script 07's ifnb
# donor IDs), so the demos also run on compute nodes without internet access.
# World-readable, since HPC and the Nextflow docker profile run tasks as the
# host user.
ENV SEURAT_CLI_CACHE=/opt/seurat-cli/cache
RUN Rscript -e "source('/app/utils/common.R'); invisible(fetch_ifnb_donor_files())" \
    && chmod -R a+rX /opt/seurat-cli

# Don't load R packages from the user's home library. Apptainer/Singularity
# mount the host home directory, where packages built for the host's R would
# otherwise be loaded ahead of the image's. The directory doesn't exist, so R
# skips it.
ENV R_LIBS_USER=/opt/seurat-cli/no-user-library

# Set default command
CMD ["Rscript", "--help"]
