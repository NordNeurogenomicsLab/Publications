# Shared file handling and explicit analysis parameters. Run scripts from the project root.
extra_library <- Sys.getenv("JIA_R_LIBRARY", unset = "")
if (nzchar(extra_library)) .libPaths(c(extra_library, .libPaths()))
suppressPackageStartupMessages(library(Seurat))
suppressPackageStartupMessages(library(ggplot2))
options(Seurat.object.assay.version = "v3")
# Integration exports multi-GB objects internally; the original notebook allowed
# 30 GB of future globals. A sequential plan avoids extra full-object worker copies.
options(future.globals.maxSize = 30000 * 1024^2)
future::plan("sequential")
data_dir <- Sys.getenv("JIA_DATA_DIR", unset = "data")
out_dir <- Sys.getenv("JIA_OUTPUT_DIR", unset = "output")
analysis_mode <- Sys.getenv("JIA_ANALYSIS_MODE", unset = "counts")
stopifnot(analysis_mode %in% c("counts", "saved"))
for (d in c("objects", "tables", "figures", "methods"))
  dir.create(file.path(out_dir, d), recursive = TRUE, showWarnings = FALSE)
source("scripts/gene_sets.R")
reference <- function() {
  d <- read.delim(file.path(data_dir, "cell_annotations.tsv"), check.names = FALSE)
  stopifnot(!anyDuplicated(d$cell_id), all(c("cell_id", "barcode", "dataset", "age",
    "reference_cluster", "historical_compartment", "mz_cluster", "progenitor_cluster") %in% names(d)))
  rownames(d) <- d$cell_id
  d
}
read_object <- function(name) {
  p <- if (analysis_mode == "saved") file.path(data_dir, "objects", name) else file.path(out_dir, "objects", name)
  if (!file.exists(p)) stop("Missing object: ", p)
  o <- readRDS(p)
  validObject(o)
  o
}
record_session <- function(step) writeLines(capture.output(sessionInfo()),
  file.path(out_dir, "methods", paste0(step, "_sessionInfo.txt")))
save_plot <- function(p, name, width = 6, height = 4.5) {
  ggsave(file.path(out_dir, "figures", paste0(name, ".png")), p,
    width = width, height = height, dpi = 300, bg = "white")
}
# Anchor search uses 30 PCs; downstream graph/UMAP use 25. These are distinct steps.
# Restricting the feature matrix precedes SCT in the cell-state analyses.
sct_integrate <- function(objects, flavor = "v2", fit_method = NULL) {
  if (flavor == "v2" || identical(fit_method, "glmGamPoi")) {
    if (!requireNamespace("glmGamPoi", quietly = TRUE))
      stop("glmGamPoi is required for this model. Install it before fitting; do not use the native fallback estimator.")
  }
  objects <- lapply(objects, function(o) {
    DefaultAssay(o) <- "RNA"
    args <- list(object = o, assay = "RNA", vst.flavor = flavor, ncells = 5000,
      variable.features.n = 3000, vars.to.regress = NULL,
      do.correct.umi = TRUE, do.scale = FALSE, do.center = TRUE,
      return.only.var.genes = TRUE, seed.use = 1448145, verbose = FALSE)
    if (!is.null(fit_method)) args$method <- fit_method
    do.call(SCTransform, args)
  })
  features <- SelectIntegrationFeatures(objects, nfeatures = 3000)
  # nfeatures=3000 is a request: a restricted panel can only supply its available genes.
  objects <- PrepSCTIntegration(objects, anchor.features = features, verbose = FALSE)
  anchors <- FindIntegrationAnchors(objects, normalization.method = "SCT",
    anchor.features = features, dims = 1:30, verbose = FALSE)
  IntegrateData(anchors, normalization.method = "SCT", dims = 1:30,
    k.weight = 100, verbose = FALSE)
}
embed_cluster <- function(o, assay, min_dist = 0.3) {
  DefaultAssay(o) <- assay
  o <- RunPCA(o, npcs = 25, seed.use = 42, verbose = FALSE)
  o <- RunUMAP(o, dims = 1:25, n.neighbors = 30, min.dist = min_dist,
    umap.method = "uwot", metric = "cosine", seed.use = 42, verbose = FALSE)
  o <- FindNeighbors(o, dims = 1:25, k.param = 20, verbose = FALSE)
  FindClusters(o, resolution = 0.8, algorithm = 1, random.seed = 0, verbose = FALSE)
}
attach_reference <- function(o, ref) {
  stopifnot(all(colnames(o) %in% ref$cell_id))
  for (field in c("reference_cluster", "historical_compartment", "mz_cluster", "progenitor_cluster"))
    o[[paste0("manuscript_", field)]] <- ref[colnames(o), field]
  o
}
read_count_matrix <- function(path, format = "10x") {
  if (format == "h5") return(Read10X_h5(path))
  # Portable exports use plain Matrix Market files; Cell Ranger's newer layout
  # uses gzip files. ReadMtx handles the former without an OS-specific gzip tool.
  if (file.exists(file.path(path, "matrix.mtx")) && file.exists(file.path(path, "features.tsv")))
    return(ReadMtx(mtx = file.path(path, "matrix.mtx"), cells = file.path(path, "barcodes.tsv"),
      features = file.path(path, "features.tsv"), feature.column = 2))
  Read10X(path)
}
