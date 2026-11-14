# Developmental differential-expression figures (Figure 3B and Figure 4A).
# The statistical contrast is E11.5/E13.5. Reverse its sign only for plotting so
# E13.5-enriched genes appear on the right; preserve the original avg_log2FC column.
# Raw Wilcoxon p-values match the original plot. Adjusted p-values are retained in
# the result tables and must be used when interpreting corrected significance.
source("scripts/analysis_common.R")
suppressPackageStartupMessages(library(ggrepel))
for (zone in c("VZ", "MZ")) {
  path <- if (analysis_mode == "saved") file.path("tables", paste0(zone, "_markers.csv")) else
    file.path(out_dir, "tables", paste0(zone, "_markers.csv"))
  d <- read.csv(path)
  d$plot_log2FC <- -d$avg_log2FC
  d$direction <- ifelse(d$plot_log2FC > 0, "High at E13.5", "High at E11.5")
  # Labels come from this contrast's own table. The exploratory helper closed over
  # VZ_markers even when plotting another compartment; do not inherit that bug.
  labels <- head(d[order(d$p_val, d$gene_name), ], 50)
  p <- ggplot(d, aes(plot_log2FC, -log10(pmax(p_val, .Machine$double.xmin)), color = direction)) +
    geom_point(size = 1, alpha = 0.5) +
    geom_text_repel(data = labels, aes(label = gene_name), size = 3, max.overlaps = Inf, seed = 17) +
    scale_color_manual(values = c("High at E11.5" = "#F8766D", "High at E13.5" = "#7F8CCB")) +
    theme_light() + labs(title = paste("E11.5 vs E13.5", zone, "markers"),
      x = "Log2 fold change (E13.5 / E11.5)", y = "-Log10 P", color = NULL) +
    theme(legend.position = "bottom")
  save_plot(p, paste0(if (zone == "VZ") "Fig3B" else "Fig4A", "_volcano"), 7, 7)
  write.csv(d, file.path(out_dir, "tables", paste0(zone, "_volcano_data.csv")), row.names = FALSE)
}
record_session("05_volcano_plots")
