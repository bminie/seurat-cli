#!/usr/bin/env Rscript
# =============================================================================
# Write a small synthetic Xenium output folder for testing 08_spatial.R
#
# 750 cells on a grid, split into three vertical tissue domains (A, B, C).
# Each domain expresses its own five marker genes (GeneA1-5, GeneB1-5,
# GeneC1-5) on top of 15 background genes, so clustering should recover the
# domains. Files follow the Xenium layout read by Seurat::LoadXenium():
#   cell_feature_matrix/          counts (MEX), incl. 2 negative control probes
#   cells.csv.gz                  cell centroids
#   transcripts.csv.gz            transcript positions and quality scores
#   cell_boundaries.csv.gz        hexagonal cell outlines
#   nucleus_boundaries.csv.gz     smaller nucleus outlines
# Counts in the matrix are the transcripts with quality >= 20, as in Xenium.
#
# Usage: Rscript tests/make_xenium_fixture.R <output_dir>
# =============================================================================

out_dir <- commandArgs(trailingOnly = TRUE)[1]
if (is.na(out_dir)) stop("Usage: Rscript tests/make_xenium_fixture.R <output_dir>")
suppressPackageStartupMessages(library(Matrix))
set.seed(1)

# Cells on a 30 x 25 grid, 10 microns apart
grid <- expand.grid(col = 0:29, row = 0:24)
n_cells <- nrow(grid)
cells <- data.frame(
  cell_id = sprintf("c%04d-1", seq_len(n_cells)),
  x_centroid = grid$col * 10 + 5 + stats::runif(n_cells, -1, 1),
  y_centroid = grid$row * 10 + 5 + stats::runif(n_cells, -1, 1)
)
domain <- cut(cells$x_centroid, breaks = c(-Inf, 100, 200, Inf), labels = c("A", "B", "C"))

markers <- paste0("Gene", rep(c("A", "B", "C"), each = 5), 1:5)
genes <- c(markers, sprintf("Bg%02d", 1:15))
controls <- c("NegControlProbe_00001", "NegControlProbe_00002")

# Expected transcripts per cell: markers high in their own domain
rate <- matrix(0.6, nrow = length(genes), ncol = n_cells, dimnames = list(genes, NULL))
for (d in c("A", "B", "C")) {
  rate[grepl(paste0("^Gene", d), genes), domain == d] <- 8
}
rate[grepl("^Bg", genes), ] <- 2

# Transcripts: one row per molecule, scattered around its cell's centroid,
# with ~10% below the quality threshold
counts_all <- matrix(stats::rpois(length(rate), rate), nrow = nrow(rate))
ctrl_counts <- matrix(stats::rbinom(2 * n_cells, 1, 0.02), nrow = 2)
feature_names <- c(genes, controls)
counts_all <- rbind(counts_all, ctrl_counts)
idx <- which(counts_all > 0, arr.ind = TRUE)
reps <- counts_all[idx]
mol_gene <- rep(idx[, 1], reps)
mol_cell <- rep(idx[, 2], reps)
n_mols <- length(mol_gene)
angle <- stats::runif(n_mols, 0, 2 * pi)
radius <- 4 * sqrt(stats::runif(n_mols))
transcripts <- data.frame(
  transcript_id = seq_len(n_mols),
  cell_id = cells$cell_id[mol_cell],
  overlaps_nucleus = as.integer(radius < 2),
  feature_name = feature_names[mol_gene],
  x_location = cells$x_centroid[mol_cell] + radius * cos(angle),
  y_location = cells$y_centroid[mol_cell] + radius * sin(angle),
  z_location = stats::runif(n_mols, 10, 20),
  qv = round(ifelse(stats::runif(n_mols) < 0.1, stats::runif(n_mols, 5, 19), stats::runif(n_mols, 20, 40)), 2)
)

# Count matrix from the transcripts that pass quality filtering
pass <- transcripts$qv >= 20
counts <- Matrix::sparseMatrix(
  i = match(transcripts$feature_name[pass], feature_names),
  j = match(transcripts$cell_id[pass], cells$cell_id),
  x = 1, dims = c(length(feature_names), n_cells))
cells$transcript_counts <- as.integer(Matrix::colSums(counts[seq_along(genes), ]))
cells$total_counts <- as.integer(Matrix::colSums(counts))
cells$cell_area <- 65
cells$nucleus_area <- 10

hexagon <- function(r) {
  a <- seq(0, 2 * pi, length.out = 7)
  data.frame(dx = r * cos(a), dy = r * sin(a))
}
outlines <- function(r) {
  h <- hexagon(r)
  data.frame(
    cell_id = rep(cells$cell_id, each = nrow(h)),
    vertex_x = rep(cells$x_centroid, each = nrow(h)) + rep(h$dx, n_cells),
    vertex_y = rep(cells$y_centroid, each = nrow(h)) + rep(h$dy, n_cells)
  )
}

# Write the folder
dir.create(file.path(out_dir, "cell_feature_matrix"), recursive = TRUE, showWarnings = FALSE)
write_gz <- function(df, name, ...) {
  con <- gzfile(file.path(out_dir, name), "w")
  on.exit(close(con))
  utils::write.csv(df, con, row.names = FALSE, ...)
}
write_tsv_gz <- function(df, path) {
  con <- gzfile(path, "w")
  on.exit(close(con))
  utils::write.table(df, con, sep = "\t", quote = FALSE, row.names = FALSE, col.names = FALSE)
}

mtx_path <- file.path(out_dir, "cell_feature_matrix", "matrix.mtx")
invisible(Matrix::writeMM(methods::as(counts, "CsparseMatrix"), mtx_path))
gzip_file <- function(path) {
  con_in <- file(path, "rb")
  con_out <- gzfile(paste0(path, ".gz"), "wb")
  writeBin(readBin(con_in, "raw", file.info(path)$size), con_out)
  close(con_in)
  close(con_out)
  unlink(path)
}
gzip_file(mtx_path)
write_tsv_gz(data.frame(id = sprintf("ENSFIX%05d", seq_along(feature_names)),
                        name = feature_names,
                        type = c(rep("Gene Expression", length(genes)),
                                 rep("Negative Control Probe", length(controls)))),
             file.path(out_dir, "cell_feature_matrix", "features.tsv.gz"))
write_tsv_gz(data.frame(cells$cell_id), file.path(out_dir, "cell_feature_matrix", "barcodes.tsv.gz"))

write_gz(cells, "cells.csv.gz")
write_gz(transcripts, "transcripts.csv.gz")
write_gz(outlines(5), "cell_boundaries.csv.gz")
write_gz(outlines(2), "nucleus_boundaries.csv.gz")

cat(sprintf("Wrote synthetic Xenium data to %s: %d cells, %d genes, %d transcripts\n",
            out_dir, n_cells, length(genes), n_mols))
