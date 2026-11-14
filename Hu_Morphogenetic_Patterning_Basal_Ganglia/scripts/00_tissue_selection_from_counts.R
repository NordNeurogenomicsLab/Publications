# Optional pre-exclusion analysis for 1538TE10H11 and 1538TE10H13.
# Use original pre-exclusion count matrices, not matrices exported from retained
# cells. The classifier must see the contaminating populations to learn exclusion.
# Source: VZ H11 vs H13 notebook, E13 cell rescue and SVM sections (2023).
# Poor LMO performance motivated using GEX cells and expression-based tissue
# filtering rather than treating unreliable lipid barcodes as ground truth.
source("scripts/analysis_common.R")
source("scripts/tissue_identity_classification.R")
datasets <- c("1538TE10H11", "1538TE10H13")
checkpoint <- Sys.getenv("JIA_TISSUE_CHECKPOINT", unset = "")
if (nzchar(checkpoint)) {
  pre <- readRDS(checkpoint)
  validObject(pre)
  write.table(data.frame(path = checkpoint, md5 = unname(tools::md5sum(checkpoint))),
    file.path(out_dir, "methods", "tissue_checkpoint.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
} else {
inputs <- read.delim(file.path(data_dir, "count_inputs_preselection.tsv"), stringsAsFactors = FALSE)
stopifnot(setequal(inputs$dataset, datasets), !anyDuplicated(inputs$dataset))
objects <- list()
for (i in seq_along(datasets)) {
  row <- inputs[inputs$dataset == datasets[i], ]
  counts <- read_count_matrix(file.path(data_dir, row$path), row$format)
  if (is.list(counts)) counts <- counts[["Gene Expression"]]
  o <- CreateSeuratObject(counts, min.cells = 0, min.features = 0)
  o[["percent.mt"]] <- PercentageFeatureSet(o, pattern = "^mt-")
  o <- subset(o, subset = nCount_RNA > 1000 & percent.mt < 20)
  o <- RenameCells(o, new.names = paste0(colnames(o), "_", i))
  o$dataset <- datasets[i]; o$age <- if (i == 1) "E11" else "E13"
  objects[[i]] <- o
}
# The preliminary library models used glmGamPoi (v1 SCT) before integration.
pre <- sct_integrate(objects, flavor = "v1", fit_method = "glmGamPoi")
pre <- embed_cluster(pre, "integrated", min_dist = 0.3)
}
result <- classify_tissue_identity(pre, file.path(out_dir, "tissue_classification"))
saveRDS(result$MGE, file.path(out_dir, "objects", "MGE_older_libraries_tissue_filtered.rds"))
# Classifier refitting is a new fit: report agreement with the accepted manuscript
# selection, preserving mismatches for inspection rather than replacing its cells.
ref <- reference()
accepted <- ref$cell_id[ref$dataset %in% datasets]
all_cells <- union(colnames(result$classified), accepted)
write.table(data.frame(cell_id = all_cells, manuscript_MGE = all_cells %in% accepted,
  classifier_MGE = all_cells %in% colnames(result$MGE)),
  file.path(out_dir, "tables", "tissue_selection_comparison.tsv"),
  sep = "\t", quote = FALSE, row.names = FALSE)
record_session("00_tissue_selection")
