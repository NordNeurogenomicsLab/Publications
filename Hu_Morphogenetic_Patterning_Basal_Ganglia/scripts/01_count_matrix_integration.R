# Count import, quality control, tissue selection, and global SCT integration.
# 1538TE10H11 and 1538TE10H13 form the curated MGE reference.
# JH4327_01_TE10H13 is integrated with this reference as a separate unit.
# This two-unit design retains developmental-age variation within the reference
# while aligning shared cell states in JH4327_01_TE10H13.
# Source: Jia_1538_Aug_2023 datasets.Rmd, integration section.
source("scripts/analysis_common.R")
ref <- reference()
inputs <- read.delim(file.path(data_dir, "count_inputs.tsv"), stringsAsFactors = FALSE)
expected <- c("1538TE10H11", "1538TE10H13", "JH4327_01_TE10H13")
stopifnot(setequal(inputs$dataset, expected), !anyDuplicated(inputs$dataset))
objects <- list()
qc <- list()
for (i in seq_len(nrow(inputs))) {
  row <- inputs[i, ]; path <- file.path(data_dir, row$path)
  counts <- read_count_matrix(path, row$format)
  if (is.list(counts)) counts <- counts[["Gene Expression"]]
  o <- CreateSeuratObject(counts, assay = "RNA", min.cells = 0, min.features = 0)
  o[["percent.mt"]] <- PercentageFeatureSet(o, pattern = "^mt-")
  before <- ncol(o)
  # Actual historical filtering code used >1000 UMIs and <20% mitochondrial UMIs.
  # These gentle thresholds retain RNA complexity while excluding cells dominated
  # by mitochondrial transcripts; they were chosen for consistency across libraries.
  # The computed 5-MAD/200-feature/15%-mitochondrial flags were diagnostics, not the
  # active subset condition. Do not apply them as extra filters.
  if (row$selection == "unfiltered") o <- subset(o, subset = nCount_RNA > 1000 & percent.mt < 20)
  wanted <- ref[ref$dataset == row$dataset, , drop = FALSE]
  missing <- setdiff(wanted$barcode, colnames(o))
  if (length(missing)) stop(row$dataset, ": ", length(missing),
    " retained barcodes absent. Use original counts/raw matrix; do not silently drop cells.")
  # The annotation table records the accepted MGE cells after POA/septum exclusion.
  # This makes the exact selection reproducible independently of classifier refitting.
  # tissue_identity_classification.R supplies the documented marker/SVM procedure.
  o <- subset(o, cells = wanted$barcode)
  o <- RenameCells(o, new.names = wanted$cell_id[match(colnames(o), wanted$barcode)])
  o$dataset <- row$dataset; o$age <- row$age
  objects[[row$dataset]] <- o
  qc[[row$dataset]] <- data.frame(dataset = row$dataset, input_cells = before, retained_cells = ncol(o))
}
write.csv(do.call(rbind, qc), file.path(out_dir, "tables", "count_import_qc.csv"), row.names = FALSE)
old <- merge(objects[["1538TE10H11"]], objects[["1538TE10H13"]])
# The stored global object's SCT models have the v1 flavor (NULL) and poisson fit.
# Here poisson is the v1 fitting-routine name; SCT still estimates and regularizes
# negative-binomial dispersion rather than assuming Poisson noise throughout.
# Explicitly retain that historical model instead of adopting the Seurat 5 v2 default.
# The standalone exploratory glmGamPoi normalization of the new library was followed
# by a new default SCTransform on both integration units; the latter defines integration.
x <- sct_integrate(list(old, objects[["JH4327_01_TE10H13"]]), flavor = "v1")
x <- embed_cluster(x, "integrated", min_dist = 0.3)
x <- attach_reference(x, ref)
write.csv(as.data.frame(table(reference = x$manuscript_reference_cluster, calculated = x$seurat_clusters)),
  file.path(out_dir, "tables", "global_cluster_correspondence.csv"), row.names = FALSE)
saveRDS(x, file.path(out_dir, "objects", "MGE_1538_E11_and_E13_with_new_E13.Rds"))
record_session("01_global_integration")
