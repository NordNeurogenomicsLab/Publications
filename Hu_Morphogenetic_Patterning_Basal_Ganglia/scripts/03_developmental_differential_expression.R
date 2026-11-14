# Developmental expression changes and cluster markers using full-gene SCT data.
# The E11.5 group is ident.1: positive avg_log2FC means higher at E11.5.
# min.pct=0.25 limits testing to genes detected in >=25% of either group;
# logfc.threshold=0.25 is the effect-size prefilter. Wilcoxon tests are cell-level
# comparisons matching the source analysis, not replicate-level pseudobulk tests.
# Seurat supplies Bonferroni-adjusted p_val_adj; plots retain raw p_val as in the
# manuscript code. Neither the integration assay nor the 52-gene panel is used
# as the genome-wide differential-expression testing matrix.
source("scripts/analysis_common.R")
ref <- reference()
targets <- strsplit(Sys.getenv("JIA_DE_TARGETS", unset = "global,MZ,VZ_SVZ"), ",")[[1]]
stopifnot(all(targets %in% c("global", "MZ", "VZ_SVZ")))
if ("global" %in% targets) {
x <- read_object("MGE_1538_E11_and_E13_with_new_E13.Rds")
DefaultAssay(x) <- "SCT"
x <- PrepSCTFindMarkers(x, verbose = FALSE)
compartments <- list(VZ = c("2", "15", "5", "12", "7", "14"),
  SVZ = c("9", "11", "1"), MZ = c("3", "8", "4", "0", "6", "16", "13", "10"))
for (zone in names(compartments)) {
  a <- ref[ref$reference_cluster %in% compartments[[zone]], ]
  cells1 <- intersect(a$cell_id[a$age == "E11"], colnames(x))
  cells2 <- intersect(a$cell_id[a$age == "E13"], colnames(x))
  d <- FindMarkers(x, ident.1 = cells1, ident.2 = cells2, assay = "SCT", slot = "data",
    test.use = "wilcox", min.pct = 0.25, logfc.threshold = 0.25, only.pos = FALSE)
  d$gene_name <- rownames(d)
  write.csv(d, file.path(out_dir, "tables", paste0(zone, "_markers.csv")), row.names = FALSE)
  write.table(data.frame(cell_id = c(cells1, cells2), age = c(rep("E11", length(cells1)), rep("E13", length(cells2)))),
    file.path(out_dir, "methods", paste0(zone, "_DE_cells.tsv")), sep = "\t", quote = FALSE, row.names = FALSE)
}
rm(x); gc()
}
for (zone in intersect(c("MZ", "VZ_SVZ"), targets)) {
  filename <- if (zone == "MZ") "MZ_reintegrated_on_52_PN_genes_08_15_2025.rds" else "VZ_SVZ_reintegrated_on_221_progenitor_genes_08_15_2025.rds"
  o <- read_object(filename)
  DefaultAssay(o) <- "SCT"
  o <- PrepSCTFindMarkers(o, verbose = FALSE)
  marker_o <- o
  if (zone == "VZ_SVZ") {
    # The source refits a full-gene SCT v2 model for progenitor cluster markers,
    # then reloads the inherited-model checkpoint for the later age comparisons.
    DefaultAssay(marker_o) <- "RNA"
    marker_o <- SCTransform(marker_o, vst.flavor = "v2", seed.use = 1448145, verbose = FALSE)
    marker_o <- PrepSCTFindMarkers(marker_o, verbose = FALSE)
    saveRDS(marker_o, file.path(out_dir, "objects", "VZ_SVZ_cluster_marker_model.rds"))
  }
  Idents(marker_o) <- "seurat_clusters"
  markers <- FindAllMarkers(marker_o, assay = "SCT", test.use = "wilcox", only.pos = TRUE,
    min.pct = 0.25, logfc.threshold = 0.25)
  write.csv(markers, file.path(out_dir, "tables", paste0(zone, "_cluster_markers.csv")), row.names = FALSE)
  rm(marker_o); gc()
  Idents(o) <- "seurat_clusters"
  results <- list()
  for (cluster in levels(Idents(o))) {
    sub <- subset(o, idents = cluster)
    n <- table(factor(sub$age, levels = c("E11", "E13")))
    # The original VZ/SVZ analysis skipped cluster 9 with only two cells at one age.
    # Express that rule by group size so a rerun's arbitrary numeric labels cannot
    # accidentally exclude a different cluster or attempt a two-cell comparison.
    if (any(n < 3)) next
    Idents(sub) <- "age"
    d <- FindMarkers(sub, ident.1 = "E11", ident.2 = "E13", assay = "SCT",
      test.use = "wilcox", min.pct = 0.25, logfc.threshold = 0.25, recorrect_umi = FALSE)
    d$gene <- rownames(d); d$cluster <- cluster
    results[[cluster]] <- d
  }
  write.csv(do.call(rbind, results), file.path(out_dir, "tables", paste0(zone, "_within_cluster_age_markers.csv")), row.names = FALSE)
}
record_session("03_differential_expression")
