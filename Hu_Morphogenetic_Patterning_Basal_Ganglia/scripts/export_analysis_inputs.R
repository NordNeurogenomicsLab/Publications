# Export raw RNA UMI counts and curated cell selections from the manuscript objects.
# RNA counts, not SCT-corrected counts or integrated residuals, are the count input.
# This preserves the historical quantification when later Cell Ranger runs differ.
# Exported matrices contain retained MGE cells, not the cells excluded before analysis.
source("scripts/analysis_common.R")
source_dir <- Sys.getenv("JIA_SOURCE_OBJECT_DIR", unset = file.path(data_dir, "objects"))
x <- readRDS(file.path(source_dir, "MGE_1538_E11_and_E13_with_new_E13.Rds"))
m <- readRDS(file.path(source_dir, "MZ_reintegrated_on_52_PN_genes_08_15_2025.rds"))
p <- readRDS(file.path(source_dir, "VZ_SVZ_reintegrated_on_221_progenitor_genes_08_15_2025.rds"))
stopifnot(all(colnames(m) %in% colnames(x)), all(colnames(p) %in% colnames(x)))
ann <- data.frame(cell_id = colnames(x), barcode = sub("_[0-9]+$", "", colnames(x)),
  dataset = x$dataset, age = x$age, reference_cluster = as.character(x$seurat_clusters),
  historical_compartment = as.character(x$Dev_classes),
  mz_cluster = as.character(m$seurat_clusters[match(colnames(x), colnames(m))]),
  progenitor_cluster = as.character(p$seurat_clusters[match(colnames(x), colnames(p))]))
write.table(ann, file.path(data_dir, "cell_annotations.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
manifest <- list()
for (ds in unique(ann$dataset)) {
  cells <- ann$cell_id[ann$dataset == ds]
  mat <- GetAssayData(x, assay = "RNA", slot = "counts")[, cells, drop = FALSE]
  stopifnot(!anyDuplicated(ann$barcode[ann$dataset == ds]))
  colnames(mat) <- ann$barcode[match(cells, ann$cell_id)]
  dest <- file.path(data_dir, "counts", ds)
  dir.create(dest, recursive = TRUE, showWarnings = FALSE)
  Matrix::writeMM(mat, file.path(dest, "matrix.mtx"))
  write.table(data.frame(rownames(mat), rownames(mat), "Gene Expression"),
    file.path(dest, "features.tsv"), sep = "\t", quote = FALSE, row.names = FALSE, col.names = FALSE)
  writeLines(colnames(mat), file.path(dest, "barcodes.tsv"))
  # Confirm that the portable files preserve every UMI and every cell/gene label.
  roundtrip <- read_count_matrix(dest)
  stopifnot(identical(dimnames(mat), dimnames(roundtrip)), Matrix::nnzero(mat - roundtrip) == 0)
  manifest[[ds]] <- data.frame(dataset = ds, age = unique(ann$age[ann$dataset == ds]),
    path = file.path("counts", ds), format = "10x", selection = "curated")
}
write.table(do.call(rbind, manifest), file.path(data_dir, "count_inputs.tsv"),
  sep = "\t", quote = FALSE, row.names = FALSE)
files <- list.files(file.path(data_dir, "counts"), recursive = TRUE, full.names = TRUE)
write.table(data.frame(path = files, md5 = unname(tools::md5sum(files))),
  file.path(data_dir, "count_files.md5.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
record_session("export_inputs")
