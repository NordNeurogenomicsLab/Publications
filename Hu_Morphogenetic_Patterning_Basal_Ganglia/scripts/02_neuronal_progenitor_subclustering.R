# Cell-state analysis using biologically selected feature panels.
# Source: August 2025 custom-gene-set sections in the methods notebook.
# Restrict features BEFORE splitting by dataset and SCT reintegration. This directs
# the embedding toward neuronal identity or progenitor maturation rather than
# library depth and unrelated highly variable expression programs.
# Preserve the full-gene RNA/SCT assays: the panel defines the embedding, not the
# universe of genes available for expression figures and differential expression.
source("scripts/analysis_common.R")
x <- read_object("MGE_1538_E11_and_E13_with_new_E13.Rds")
ref <- reference()
specs <- list(
  MZ = list(clusters = c("0", "1", "3", "8", "13", "16", "6", "4"), genes = proj_neu_genes,
    filename = "MZ_reintegrated_on_52_PN_genes_08_15_2025.rds"),
  VZ_SVZ = list(clusters = c("2", "15", "7", "14", "5", "12", "11", "9"), genes = progenitor_genes,
    filename = "VZ_SVZ_reintegrated_on_221_progenitor_genes_08_15_2025.rds"))
requested_states <- strsplit(Sys.getenv("JIA_CELL_STATES", unset = "MZ,VZ_SVZ"), ",")[[1]]
stopifnot(all(requested_states %in% names(specs)))
for (name in requested_states) {
  sp <- specs[[name]]
  # Numeric clusters refer to the manuscript annotation, not a new run's arbitrary IDs.
  cells <- ref$cell_id[ref$reference_cluster %in% sp$clusters]
  full <- subset(x, cells = cells)
  genes <- intersect(sp$genes, rownames(full[["RNA"]]))
  writeLines(genes, file.path(out_dir, "methods", paste0(name, "_used_panel_genes.txt")))
  writeLines(setdiff(sp$genes, genes), file.path(out_dir, "methods", paste0(name, "_absent_panel_genes.txt")))
  # Start from RNA counts to avoid inheriting stale full-data SCT/integration models.
  restricted <- CreateSeuratObject(GetAssayData(full, assay = "RNA", slot = "counts")[genes, , drop = FALSE],
    meta.data = full[[]], min.cells = 0, min.features = 0)
  fitted <- sct_integrate(SplitObject(restricted, split.by = "dataset"), flavor = "v2")
  fitted <- embed_cluster(fitted, "integrated", min_dist = 0.01)
  # min.dist=0.01 emphasizes local subtype structure; resolution=0.8 controls
  # Louvain granularity. No cell-cycle regression was used for these panel analyses.
  stopifnot(setequal(colnames(full), colnames(fitted)))
  fitted <- fitted[, colnames(full)]
  full[["panel_integrated"]] <- fitted[["integrated"]]
  full[["pca"]] <- fitted[["pca"]]
  full[["umap"]] <- fitted[["umap"]]
  DefaultAssay(full[["pca"]]) <- "panel_integrated"
  DefaultAssay(full[["umap"]]) <- "panel_integrated"
  full$seurat_clusters <- fitted$seurat_clusters
  Idents(full) <- "seurat_clusters"
  # The saved MZ checkpoint includes the later full-gene v2 marker model. The
  # progenitor checkpoint was saved earlier and retains inherited global v1 SCT
  # models; its separate v2 cluster-marker fit belongs to script 03, not this file.
  if (name == "MZ") {
    DefaultAssay(full) <- "RNA"
    full <- SCTransform(full, vst.flavor = "v2", seed.use = 1448145, verbose = FALSE)
    full <- PrepSCTFindMarkers(full, verbose = FALSE)
  }
  if (name == "VZ_SVZ") {
    source("scripts/cell_cycle_annotation.R")
    full <- annotate_cell_cycle(full)
  }
  full <- attach_reference(full, ref)
  saveRDS(full, file.path(out_dir, "objects", sp$filename))
  write.csv(full[[]], file.path(out_dir, "tables", paste0(name, "_cell_metadata.csv")))
  writeLines(capture.output(fitted@commands), file.path(out_dir, "methods", paste0(name, "_commands.txt")))
}
record_session("02_cell_state_analysis")
