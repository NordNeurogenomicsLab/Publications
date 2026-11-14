# Cell-cycle annotation from the custom progenitor analysis (not regression).
# Source uses Seurat cc.genes converted to mouse title-case symbols and a modified
# phase rule: G1 when BOTH S and G2M module scores are <0.1. Keeping this threshold
# matters because the default Seurat phase assignment uses a different boundary.
annotate_cell_cycle <- function(o) {
  mouse_symbol <- function(g) paste0(toupper(substr(g, 1, 1)), tolower(substring(g, 2)))
  genes <- list(mouse_symbol(cc.genes$s.genes), mouse_symbol(cc.genes$g2m.genes))
  writeLines(genes[[1]], file.path(out_dir, "methods", "S_phase_genes.txt"))
  writeLines(genes[[2]], file.path(out_dir, "methods", "G2M_phase_genes.txt"))
  DefaultAssay(o) <- "SCT"
  # Module scores subtract expression-matched control genes. The control count
  # equals the smaller phase-gene set, matching the source function.
  o <- AddModuleScore(o, features = genes, assay = "SCT", name = "Cell.Cycle",
    ctrl = min(lengths(genes)), nbin = 24, seed = 1)
  scores <- as.matrix(o[[]][, c("Cell.Cycle1", "Cell.Cycle2")])
  o$S.Score <- scores[, 1]; o$G2M.Score <- scores[, 2]
  o$Phase <- apply(scores, 1, function(s) {
    if (all(s < 0.1)) "G1" else if (sum(s == max(s)) > 1) "Undecided" else c("S", "G2M")[which.max(s)]
  })
  o
}
