# Expression and composition of MGE developmental populations (Figures 3 and 4).
# Use SCT log-normalized expression for feature plots, preserving the full-gene
# assay after panel-based embedding. VZ expression and broad VZ DE use the same
# six-cluster compartment defined in the original developmental annotation code.
# The source also requests individual DE in four VZ clusters; that focused request
# does not redefine the whole VZ population for the expression panels.
source("scripts/analysis_common.R")
suppressPackageStartupMessages(library(patchwork))
ref <- reference()
x <- read_object("MGE_1538_E11_and_E13_with_new_E13.Rds")
x$plot_age <- factor(x$age, levels = c("E11", "E13"), labels = c("E11.5", "E13.5"))
save_plot(DimPlot(x, group.by = "plot_age", cols = c("E11.5" = "#F8766D", "E13.5" = "#00BFC4")), "Fig3A_age")
DefaultAssay(x) <- "SCT"
for (gene in c("Hes1", "Ascl1", "Dcx"))
  save_plot(FeaturePlot(x, features = gene, cols = c("lightgrey", "#7561A8")), paste0("Fig3A_", gene))
# Hes1, Ascl1, Dcx orient progenitor, neurogenic intermediate, and postmitotic states.
vz_cells <- ref$cell_id[ref$reference_cluster %in% c("2", "15", "5", "12", "7", "14")]
vz <- subset(x, cells = vz_cells)
for (gene in c("Id4", "Mest", "Dlx2", "Nfib"))
  save_plot(FeaturePlot(vz, features = gene, split.by = "plot_age", keep.scale = "feature", cols = c("lightgrey", "#7561A8")),
    paste0("Fig3_split_", gene), 10, 4)
m <- read_object("MZ_reintegrated_on_52_PN_genes_08_15_2025.rds")
m$plot_age <- factor(m$age, levels = c("E11", "E13"), labels = c("E11.5", "E13.5"))
save_plot(DimPlot(m, group.by = "plot_age", cols = c("E11.5" = "#00BFC4", "E13.5" = "#F8766D")), "Fig4D_age")
save_plot(DimPlot(m, label = TRUE), "Fig4D_clusters")
DefaultAssay(m) <- "SCT"
for (gene in c("Lhx8", "Mest", "Npy", "Mafb"))
  save_plot(FeaturePlot(m, features = gene, split.by = "plot_age", keep.scale = "feature", cols = c("lightgrey", "#7561A8")),
    paste0("Fig4_split_", gene), 10, 4)
# Supplementary neuronal identity panels use the same full-gene SCT assay and
# neuronal embedding. Markers are drawn from the anatomical reference (Table S4).
for (gene in c("Nr2f1", "Id4", "Ntrk1", "Isl1", "Gbx1", "Cdh8", "Th"))
  save_plot(FeaturePlot(m, features = gene, cols = c("lightgrey", "#7561A8")), paste0("FigS7_", gene))
# Denominator is all MZ cells at a given age. These are descriptive cell fractions;
# with one E11.5 library, cells are not independent biological sample replicates.
d <- as.data.frame(table(cluster = m$seurat_clusters, age = m$age))
d$percent <- 100 * d$Freq / ave(d$Freq, d$age, FUN = sum)
write.csv(d, file.path(out_dir, "tables", "MZ_cluster_age_proportions.csv"), row.names = FALSE)
save_plot(ggplot(d, aes(cluster, percent, fill = age)) + geom_col() +
  scale_fill_manual(values = c(E11 = "#00BFC4", E13 = "#F8766D"), labels = c(E11 = "E11.5", E13 = "E13.5")) + theme_classic() +
  labs(y = "Percent of MZ cells within age", x = "Seurat cluster"), "FigS7A_cluster_proportions", 9, 4)
# Eight anatomical names in the manuscript require a separate author-confirmed
# annotation key. Retain numeric identities until that key is supplied; expression
# phenotypes do not establish spatial identity by themselves.
for (zone in c("MZ", "VZ_SVZ")) {
  path <- file.path(out_dir, "tables", paste0(zone, "_cluster_markers.csv"))
  if (!file.exists(path)) next
  markers <- read.csv(path)
  o <- if (zone == "MZ") m else readRDS(file.path(out_dir, "objects", "VZ_SVZ_cluster_marker_model.rds"))
  DefaultAssay(o) <- "SCT"
  # Top positive cluster markers describe transcriptional identity after clustering.
  # Restrict to genes with stored scaled values for DoHeatmap; do not impute genes.
  genes <- unique(unlist(lapply(split(markers, markers$cluster), function(z) head(z$gene, 20))))
  genes <- intersect(genes, rownames(GetAssayData(o, assay = "SCT", slot = "scale.data")))
  if (length(genes)) save_plot(DoHeatmap(o, features = genes) + theme(axis.text.y = element_text(size = 6)),
    paste0(zone, "_marker_heatmap"), 14, 20)
}
record_session("04_expression_figures")
