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
# Singularity (pull from Docker Hub after pushing):
#   singularity pull docker://bminie/seurat-cli:latest
#
# =============================================================================

# R 4.5.3 matches the conda environments. rocker pins R packages to a dated
# Posit Package Manager snapshot for this R version, so the installs below use
# the image's default repos (no repos= override) to stay compatible with R.
FROM rocker/r-ver:4.5.3

LABEL maintainer="bminie"
LABEL description="Seurat CLI - Command-line tools for single-cell RNA-seq analysis"
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
    libfreetype6-dev \
    libpng-dev \
    libtiff5-dev \
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

# Fail the build if any package is missing. install.packages() only warns when
# a package fails to install, so without this check the image could build
# without Seurat.
RUN Rscript -e " \
    pkgs <- c('optparse', 'dplyr', 'ggplot2', 'patchwork', 'Matrix', 'future', \
              'viridis', 'hdf5r', 'png', 'Seurat', 'SeuratObject', 'sctransform', \
              'glmGamPoi', 'SeuratData', 'pbmc3k.SeuratData', 'ifnb.SeuratData', \
              'presto', 'MAST', 'DESeq2', 'SingleCellExperiment', 'harmony', \
              'batchelor', 'scran', 'scater'); \
    missing <- pkgs[!vapply(pkgs, requireNamespace, logical(1), quietly = TRUE)]; \
    if (length(missing) > 0) { \
      message('Missing R packages: ', paste(missing, collapse = ', ')); quit(status = 1) \
    }; \
    cat('All', length(pkgs), 'R packages available. Seurat', \
        as.character(packageVersion('Seurat')), '\n')"

# Create working directory
WORKDIR /app

# Copy project files
COPY scripts/ /app/scripts/
COPY utils/ /app/utils/
COPY tests/ /app/tests/

# Make scripts accessible from PATH
ENV PATH="/app/scripts:${PATH}"

# Set default command
CMD ["Rscript", "--help"]
