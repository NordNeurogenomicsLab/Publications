# Semantic palettes sampled from the final manuscript RGB render and
# original saved objects.  Values in original plot data are never
# modified by this styling layer.
MS <- list(
  cell = c(
    VZ = "#F499C1",
    SVZ = "#F26D70",
    ImN = "#9695C9",
    L6N = "#58B947",

    L5N = "#752676"
  ),
  category = c(
    Early = "#FDD700",
    Intermediate = "#87CDEA",
    Late = "#F2634A",
    Constitutive = "#989999",

    Multimodal = "#19AE69"
  ),
  assay = c(ATAC = "#545454", AC = "#58B947", ME3 = "#FF3D2E"),

  histology = c(WT = "#488CCB", Mut = "#7BD3F7", KO = "#7BD3F7"),

  histology_ext7_j = c(WT = "#B19ACB", KO = "#9273B4", Mut = "#9273B4"),

  deletion = c(WT = "#941C52", het = "#009192", KO = "#488CCB"),

  rna = c(WT = "#2178B4", Mut = "#E21E26"),

  vista = c(VZ = "#EC4696", `non-VZ` = "#A5A4D2"),

  inferno = c("#000004", "#57106E", "#BB3755", "#F98D0A", "#FCFFA4"),
  go = c("#FFFFFF", "#FCBBA1", "#FB6A4A", "#DE2D26"),
  feature = c("#D3D3D3", "#514F9C"),
  gene_heat = c("#8A2F91", "#000000", "#F8EC1D")
)
theme_set(
  theme_classic(base_size = 10,
    base_family = "Arial") +
    theme(
      plot.title = element_text(size = 11, face = "bold"),

      axis.text = element_text(color = "black"),

      strip.background = element_rect(fill = "grey92", color = NA),

      strip.text = element_text(face = "bold"),

      legend.key.size = grid::unit(3, "mm")
    )
)
# One manuscript palette for cluster points and the Figure 3 header strips.
cluster_colors <- setNames(c(
  "#1DBDC2", "#0FB9E7", "#49A1D9", "#838EC7",
  "#A681BA", "#C376B1", "#E06CA9", "#F069A0",
  "#7CAF41", "#24B34B", "#2EB567", "#2EB998",
  "#F3766E", "#E48624", "#CB972B", "#A9A335"
), clusters)
violin_colors <- setNames(
  (scales::hue_pal())(length(clusters)),
  clusters
)
umap_plot <- function(d, feature = NULL) {
  xy <- grep(
    "umap.*[12]$|UMAP[12]$", names(d),
    ignore.case = TRUE,
    value = TRUE
  )
  stopifnot(
    length(xy) >= 2
  )
  d$Cluster <- factor(d$Cluster, levels = clusters)
  if (is.null(feature)) {
    centers <- d %>%
      group_by(Cluster) %>%
      summarise(
        X = median(.data[[xy[1]]]),
        Y = median(.data[[xy[2]]]),
        .groups = "drop"
      )
    p <- ggplot(d, aes(.data[[xy[1]]], .data[[xy[2]]], color = Cluster)) +
      geom_point(size = 0.35) +
      scale_color_manual(values = cluster_colors) +
      ggrepel::geom_text_repel(
        data = centers,
        aes(X, Y, label = gsub("_", "", Cluster)),

        inherit.aes = FALSE,
        size = 3,
        color = "black",
        seed = 1,
        max.overlaps = Inf,
        box.padding = 0.15,
        min.segment.length = 0,
        segment.size = 0.2
      ) +
      theme(legend.position = "none")
  } else {
    # Draw larger, higher-signal points last without removing observations.
    d <- d[order(d[[feature]]), , drop = FALSE]
    p <- ggplot(d, aes(.data[[xy[1]]], .data[[xy[2]]], color = .data[[feature]])) +
    geom_point(size = 0.40) +
    scale_color_gradient(low = MS$feature[1],
      high = MS$feature[2]) +
    labs(
      title = if (feature == "Eomes")
        "Tbr2" else sub("_ATAC_TFIDF$", "", feature), color = NULL
    )
  }
  module_feature <- !is.null(feature) && grepl("^A[0-9]+$", feature)
  if (module_feature) {
    number <- as.integer(sub("A", "", feature))
    category <- if (number <= 3) "Early" else if (number <= 11) "Intermediate" else if (number <= 15) "Late" else "Constitutive"
    p <- p + annotate("segment", x = -Inf, xend = Inf, y = Inf, yend = Inf,
      color = unname(MS$category[category]), linewidth = 2)
  }
  p + coord_equal(clip = if (module_feature) "off" else "on") + labs(x = "UMAP 1",
    y = "UMAP 2")
}
export_heatmap <- function(
  mat,
  name,
  width = 7,
  height = 5,
  cluster_rows = FALSE,
  cluster_cols = FALSE,

  zscore = FALSE,
  limit = NULL,
  note = "Recreated matrix; see notebook for scaling and order"
) {
  mat <- as.matrix(mat)
  storage.mode(mat) <- "double"
  stopifnot(
    nrow(mat) >
      0, ncol(mat) > 0
  )
  if (zscore) {
    mat <- t(scale(t(mat)))
    mat[!is.finite(mat)] <- 0
  }
  if (!is.null(limit)) {
    mat[mat < -limit] <- -limit
    mat[mat > limit] <- limit
  }
  colors <- if (grepl("HOMER", name))
    viridisLite::viridis(101) else if (grepl("GO", name))
    colorRampPalette(MS$go)(101) else if (grepl("Fig1/|WGCNA/", name))
    colorRampPalette(MS$inferno)(101) else if (grepl("Quantized|Ext6/", name))
    rev(RColorBrewer::brewer.pal(11, "RdYlBu")) else colorRampPalette(c("#2166AC", "#FFFFBF", "#D73027"))(101)

  if (grepl("Fig2/", name) && grepl("GO|HOMER", name)) {
    assay <- if (grepl("2b_", name)) "AC" else if (grepl("2c_", name)) "ME3" else "ATAC"
    labels <- manuscript_module_labels(colnames(mat), assay)
    categories <- manuscript_module_categories(colnames(mat), assay)
    names(categories) <- colnames(mat)
    annotation <- ComplexHeatmap::HeatmapAnnotation(
      Category = categories,
      col = list(Category = MS$category),
      show_annotation_name = FALSE,
      show_legend = FALSE,
      simple_anno_size = grid::unit(2, "mm")
    )
    is_motif <- grepl("HOMER", name)
    display_rows <- rownames(mat)
    if (!is_motif) {
      replacements <- c("positive regulation" = "pos. reg.", "negative regulation" = "neg. reg.",
        "regulation" = "reg.", "differentiation" = "diff.", "development" = "dev.")
      for (term in names(replacements)) {
        display_rows <- gsub(term, replacements[[term]], display_rows, fixed = TRUE)
      }
      display_rows <- gsub("fibroblast growth factor receptor signaling pathway", "FGFR signaling", display_rows, fixed = TRUE)
      display_rows <- gsub("Wnt signaling pathway, planar cell polarity pathway", "Wnt signaling (PCP)", display_rows, fixed = TRUE)
      display_rows <- gsub("transcription involved in cell fate commitment", "transcription (cell fate commitment)", display_rows, fixed = TRUE)
    }
    maximum <- if (is_motif) 30 else 5
    palette <- if (is_motif) viridisLite::viridis(101) else colorRampPalette(MS$go)(101)
    heatmap <- ComplexHeatmap::Heatmap(
      mat,
      name = if (is_motif) "-ln(P)" else "Fold enrichment",
      col = circlize::colorRamp2(seq(0, maximum, length.out = 101), palette),
      cluster_rows = FALSE,
      cluster_columns = FALSE,
      show_row_dend = FALSE,
      show_column_dend = FALSE,
      row_names_side = "right",
      row_labels = display_rows,
      row_names_gp = grid::gpar(fontsize = if (is_motif) 9 else 8),
      row_names_max_width = grid::unit(width * 0.42, "in"),
      width = grid::unit(width * if (is_motif) 0.56 else 0.40, "in"),
      column_labels = labels,
      column_names_gp = grid::gpar(fontsize = 9, fontface = "bold"),
      column_names_rot = 90,
      bottom_annotation = annotation,
      heatmap_legend_param = list(
        title_gp = grid::gpar(fontsize = 8),
        labels_gp = grid::gpar(fontsize = 8),
        at = if (is_motif) c(0, 15, 30) else c(0, 1, 2, 3, 4, 5)
      ),
      rect_gp = grid::gpar(col = NA)
    )
    export_complex_heatmap(heatmap, name, width, height, note)
    return(invisible(mat))
  }
  p <- pheatmap::pheatmap(
    mat,
    cluster_rows = cluster_rows,
    cluster_cols = cluster_cols,
    silent = TRUE,

    border_color = NA,
    color = colors,
    fontsize = 8,
    show_rownames = nrow(mat) <=
      80, show_colnames = ncol(mat) <= 50
  )
  export_panel(p$gtable, name, width, height, note)
  invisible(mat)
}
volcano <- function(d,
  x = "Signed_log2FC",
  y = "neg_log10_P") {
  d$Class <- ifelse(
    d$FDR < 0.1, ifelse(d[[x]] > 0, "Up, FDR < 0.10", "Down, FDR < 0.10"),
    ifelse(
      d[["p.value"]] < 0.05, ifelse(d[[x]] > 0, "Up, P < 0.05", "Down, P < 0.05"),
      "Not significant"
    )
  )
  cols <- c(
    `Up, FDR < 0.10` = "#B51E24", `Down, FDR < 0.10` = "#4A78B7", `Up, P < 0.05` = "#F4A3A5",
    `Down, P < 0.05` = "#B1CCE3", `Not significant` = "#BFBFBF"
  )
  ggplot(d, aes(.data[[x]], .data[[y]], color = Class)) +
    geom_point(size = 0.45,
      alpha = 0.7) +
    geom_hline(
      yintercept = -log10(0.05),
      linetype = 2,
      color = "grey45",
      linewidth = 0.3
    ) +
    scale_color_manual(values = cols) +
    labs(
      x = "log2 fold change (concentration difference)",
      y = "-log10(P)",

      color = NULL
    ) +
    theme(panel.grid.major = element_line(color = "grey94",
      linewidth = 0.25))
}
plot_rpkm <- function(d, name) {
  stopifnot(
    all(
      abs(
        d$RPKM - as.double(d$Raw_count) *
          1e+09/(as.double(d$Gene_length_bp) *
          as.double(d$Library_sum_gene_counts))
      ) <
        1e-08
    )
  )
  d$Genotype <- factor(d$Genotype, levels = c("WT", "Mut"))
  p <- ggplot(d, aes(Genotype, RPKM, color = Genotype)) +
    geom_boxplot(width = 0.6,
      outlier.shape = NA,
      linewidth = 0.4) +
    geom_point(
      position = position_jitter(width = 0.08, seed = 1),

      size = 1.4
    ) +
    facet_wrap(~Gene + FlashTag_population,
      scales = "free_y",
      ncol = 4) +
    scale_color_manual(values = MS$rna) +
    labs(x = NULL,
      y = "RPKM",
      color = NULL) +
    theme(legend.position = "none")
  export_panel(
    p, name, 10, max(
      3, ceiling(length(unique(paste(d$Gene, d$FlashTag_population)))/4) *
        2.4
    ),
    note = "Manuscript WT blue/cKO red; reconstructed RPKM with original RNA model statistics"
  )
}


# Gene annotations are selected from the manuscript artwork, not from a
# new significance ranking. All coordinates and statistics remain supplied.
rna_volcano <- function(d, population) {
  labels <- read.csv(file.path(REPO, "provenance", "manuscript_gene_labels.csv"))
  selected <- labels$Gene[labels$Population == population]
  d$Direction <- ifelse(d$logFC > 0, "Up", "Down")
  d$Class <- ifelse(
    d$plot_group == "Not significant",
    "Not significant",
    paste(d$plot_group, d$Direction)
  )
  palette <- c(
    `FDR<0.05 Down` = "#488CCB",
    `FDR<0.05 Up` = "#D84E4A",
    `P<0.05 Down` = "#A6D3F1",
    `P<0.05 Up` = "#F3A6A1",
    `Not significant` = "#969696"
  )
  d <- d[order(d$plot_group == "FDR<0.05"), ]
  labeled <- d[d$gene_name %in% selected, ]
  stopifnot(nrow(labeled) == length(selected))
  top <- if (population == "VZ") 36 else 39
  counts <- data.frame(
    x = c(-3.8, 0.4),
    y = top * 0.87,
    label = c(
      sprintf("FDR < 0.05\n%d genes down\n%d genes up",
              sum(d$FDR < 0.05 & d$logFC < 0),
              sum(d$FDR < 0.05 & d$logFC > 0)),
      sprintf("P < 0.05\n%d genes down\n%d genes up",
              sum(d$PValue < 0.05 & d$logFC < 0),
              sum(d$PValue < 0.05 & d$logFC > 0))
    )
  )
  ggplot(d, aes(plot_logFC_clip4, neg_log10_PValue, color = Class)) +
    geom_point(size = 0.42, alpha = 0.85) +
    geom_hline(
      yintercept = -log10(0.05),
      linetype = 2,
      color = "grey45",
      linewidth = 0.25
    ) +
    ggrepel::geom_text_repel(
      data = labeled,
      aes(label = gene_name),
      color = unname(ifelse(labeled$Class == "Not significant", "black", palette[labeled$Class])),
      seed = 5,
      size = 2.1,
      fontface = "bold",
      box.padding = 0.12,
      point.padding = 0.05,
      min.segment.length = 0,
      segment.size = 0.15,
      max.overlaps = Inf,
      max.time = 3,
      force = 1.4,
      show.legend = FALSE
    ) +
    geom_text(
      data = counts,
      aes(x, y, label = label),
      inherit.aes = FALSE,
      hjust = 0,
      vjust = 1,
      size = 2.2,
      lineheight = 1.1
    ) +
    scale_color_manual(values = palette) +
    scale_x_continuous(breaks = seq(-4, 4, 2)) +
    coord_cartesian(xlim = c(-4.1, 4.1), ylim = c(0, top)) +
    labs(
      title = if (population == "VZ") "FlashTag-positive cells" else "FlashTag-negative cells",
      x = "log2 fold change",
      y = "-log10(P)"
    ) +
    theme(
      legend.position = "none",
      plot.title = element_text(size = 8, face = "bold"),
      axis.text = element_text(size = 7),
      axis.title = element_text(size = 8),
      plot.margin = margin(3, 3, 3, 3)
    )
}

# Resolve source module colors to manuscript names without changing columns.
manuscript_module_labels <- function(modules, assay) {
  dictionary <- read.csv(file.path(
    REPO, "data", "upstream", "Module_renaming_dictionary.csv"
  ), fileEncoding = "UTF-8-BOM")
  stopifnot("Data_type" %in% names(dictionary))
  dictionary <- dictionary[dictionary$Data_type == assay, ]
  result <- dictionary$Module_figure_names_short[
    match(modules, dictionary$Fine_clusters)
  ]
  result[is.na(result)] <- modules[is.na(result)]
  sub("^TMK", "TM", result)
}

manuscript_module_categories <- function(modules, assay) {
  labels <- manuscript_module_labels(modules, assay)
  number <- as.integer(gsub("[^0-9]", "", labels))
  if (assay == "ATAC") {
    stage <- ifelse(number <= 3, "Early", ifelse(number <= 11, "Intermediate",
                     ifelse(number <= 15, "Late", "Constitutive")))
  } else {
    stage <- rep("Intermediate", length(modules))
    if (assay == "AC") {
      stage[modules %in% c("turquoise", "Weak_early_VZ", "Weak_late_VZ", "brown")] <- "Early"
      stage[modules %in% c("yellow", "red", "Strong_late_1", "Strong_late_2")] <- "Late"
      stage[modules %in% c("purple", "Constitutive_strong_1", "Constitutive_strong_2")] <- "Constitutive"
    } else {
      stage[modules %in% c("green", "red", "greenyellow")] <- "Early"
      stage[modules %in% c("blue", "turquoise", "Weak_late")] <- "Late"
      stage[modules %in% c("Constitutive_medium", "Constitutive_strong")] <- "Constitutive"
    }
  }
  stage
}

manuscript_go_bars <- function(d, color) {
  d <- d[order(d$max_FE, decreasing = TRUE), ]
  d$Term <- factor(d$name, levels = rev(d$name))
  d$Display_label <- d$name
  substitutions <- c(
    "positive regulation" = "pos. reg.",
    "negative regulation" = "neg. reg.",
    "regulation" = "reg.",
    "differentiation" = "diff.",
    "development" = "dev.",
    "organization" = "org."
  )
  for (word in names(substitutions)) {
    d$Display_label <- gsub(word, substitutions[[word]], d$Display_label, fixed = TRUE)
  }
  ggplot(d, aes(max_FE, Term)) +
    geom_col(fill = color, width = 0.86) +
    geom_text(
      aes(x = 0.015 * max(max_FE), label = Display_label),
      hjust = 0,
      color = "white",
      size = 2.3,
      fontface = "bold"
    ) +
    scale_x_continuous(expand = expansion(mult = c(0, 0.04))) +
    labs(x = "Fold enrichment", y = NULL) +
    theme(
      axis.text.y = element_blank(),
      axis.ticks.y = element_blank(),
      axis.text.x = element_text(size = 7),
      axis.title.x = element_text(size = 8),
      plot.margin = margin(2, 2, 2, 2)
    )
}

# Cluster-header colors sampled from the original Figure 3d RGB rendering.
MS$heatmap_cluster <- setNames(c(
  "#1DBDC2", "#0FB9E7", "#49A1D9", "#838EC7",
  "#A681BA", "#C376B1", "#E06CA9", "#F069A0",
  "#7CAF41", "#24B34B", "#2EB567", "#2EB998",
  "#F3766E", "#E48624", "#CB972B", "#A9A335"
), clusters)
MS$histology_ext7_j <- c(WT = "#6A51A1", KO = "#A16FAF", Mut = "#A16FAF")

export_complex_heatmap <- function(heatmap, name, width, height, note = "Saved numerical values; manuscript display labels") {
  path <- file.path(REPO, "panels", name)
  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
  grDevices::cairo_pdf(paste0(path, ".pdf"), width = width, height = height)
  tryCatch(ComplexHeatmap::draw(heatmap), finally = grDevices::dev.off())
  grDevices::png(
    paste0(path, ".png"),
    width = width,
    height = height,
    units = "in",
    res = 180,
    type = "cairo",
    bg = "white"
  )
  tryCatch(ComplexHeatmap::draw(heatmap), finally = grDevices::dev.off())
  panel_registry[[name]] <<- data.frame(panel = name, note = note, width = width, height = height)
  invisible(heatmap)
}

footprint_plot <- function(d, comparison) {
  protection <- grep("^Protection_Score_", names(d), value = TRUE)
  total_counts <- grep("^TC_", names(d), value = TRUE)
  stopifnot(length(protection) == 2, length(total_counts) == 2)
  activity_vz <- d[[protection[1]]] + d[[total_counts[1]]]
  activity_other <- d[[protection[2]]] + d[[total_counts[2]]]
  # Native HINT activity is a difference. The manuscript displays a log2 ratio.
  stopifnot(max(abs(activity_other - activity_vz - d$TF_Activity)) < 1e-7)
  valid <- is.finite(activity_vz) & is.finite(activity_other) & activity_vz > 0 & activity_other > 0
  log2_fc <- rep(NA_real_, nrow(d))
  log2_fc[valid] <- log2(activity_other[valid] / activity_vz[valid])
  if ("log2_Fold_change" %in% names(d)) {
    stopifnot(all(is.na(d$log2_Fold_change) == !valid),
      max(abs(d$log2_Fold_change[valid] - log2_fc[valid])) < 1e-12)
  }
  d$log2_Fold_change <- log2_fc
  d$Symbol <- sub("_MOUSE.*", "", d$Motif)
  comparator <- switch(comparison, VZ_E12_vs_SVZ_E12 = "SVZ",
    VZ_E12_vs_ImN_E13 = "ImN", VZ_E12_vs_L6N_P2 = "L6N", VZ_E12_vs_L5N_P2 = "L5N")
  # Undefined ratios are retained in source data, without a fabricated pseudocount.
  plotted <- d[valid, , drop = FALSE]
  significant <- plotted$P_values < 0.05 & abs(plotted$log2_Fold_change) > 0.1
  plotted$Display_group <- ifelse(significant,
    ifelse(plotted$log2_Fold_change > 0, comparator, "VZ"), "Not significant")
  manuscript_labels <- list(
    SVZ = c("USF2", "E2F5", "E2F4", "E2F2", "SP2", "NRF1", "SP1", "SP3", "KLF3", "TFDP1", "E2F1", "SP4", "TFAP2D", "SOX9", "SP5", "SOX10", "OLIG2", "TWIST1", "ATOH1", "MEIS1", "PBX1", "NEUROD1", "IKZF1"),
    ImN = c("POU4F2", "PRRX2", "LHX2", "SOX2", "SOX10", "MAX", "ATOH1", "AHR", "TWIST1", "NEUROD2", "NEUROD1", "MYC", "BHLHA15", "ARNT", "NEUROG2"),
    L6N = c("POU4F2", "VSX2", "PAX3", "LHX2", "LHX3", "NOBOX", "ALX1", "ZBTB33", "NRF1", "CTCF", "AHR", "MAX", "RFX1", "CTCFL", "ARNT", "MYCN"),
    L5N = c("MBD2", "SP2", "E2F2", "SP1", "POU4F2", "SP3", "KLF3", "VSX2", "E2F4", "E2F5", "LHX2", "NOBOX", "PBX1", "LHX3", "E2F1", "CEBPZ", "SP5", "TFAP2D", "MEF2A", "MEF2C", "MEF2D")
  )
  # Display-only selection transcribed from the manuscript; P values are untouched.
  labeled <- plotted[plotted$Symbol %in% manuscript_labels[[comparator]], , drop = FALSE]
  labeled <- labeled[order(labeled$P_values), , drop = FALSE]
  labeled <- labeled[!duplicated(labeled$Symbol), , drop = FALSE]
  manuscript_range <- if (comparator == "SVZ") c(-1, 2) else c(-2, 1.5)
  # Extend the frame when necessary instead of clipping a finite observation.
  display_range <- c(min(manuscript_range[1], min(plotted$log2_Fold_change) - 0.05),
    max(manuscript_range[2], max(plotted$log2_Fold_change) + 0.05))
  palette <- c(VZ = MS$cell[["VZ"]], setNames(MS$cell[[comparator]], comparator),
    `Not significant` = "#BFBFBF")
  p <- ggplot(plotted, aes(log2_Fold_change,
      -log10(pmax(P_values, .Machine$double.xmin)), color = Display_group)) +
    geom_point(size = 0.65) +
    geom_hline(yintercept = -log10(0.05), linewidth = 0.25, color = "grey55") +
    geom_vline(xintercept = c(-0.1, 0.1), linewidth = 0.25, color = "grey55") +
    ggrepel::geom_text_repel(data = labeled, aes(label = Symbol), seed = 2,
      size = 2.2, color = "#24364B", fontface = "bold", max.overlaps = Inf, box.padding = 0.12,
      min.segment.length = 0, segment.size = 0.15, max.time = 3) +
    scale_color_manual(values = palette) +
    coord_cartesian(xlim = display_range) +
    labs(title = paste("VZ vs", comparator), x = "log2 fold change", y = "-log10(P)") +
    theme(legend.position = "none", plot.title = element_text(size = 10),
      axis.title = element_text(size = 9), axis.text = element_text(size = 8))
  point_layer <- ggplot_build(p)$data[[1]]
  stopifnot(nrow(point_layer) == sum(valid))
  p
}

# Fill colors sampled from the original Extended Data Figure 5 rendering.
MS$qc <- c(nCount_RNA = "#F3756D", nFeature_RNA = "#7BAF41",
           nCount_peaks = "#19BDC2", nFeature_peaks = "#A780BA")
MS$marker <- c(
  Hes1 = "#F3756D", Sox2 = "#DA8E27", Pax6 = "#ACA234", Hmga2 = "#63B346",
  Eomes = "#29B45C", Neurog2 = "#2ABAA3", Insm1 = "#00BADE", Tbr1 = "#479DD7",
  Neurod2 = "#9C85BE", Bhlhe22 = "#C473AF", Slc17a6 = "#F066A6"
)
MS$pca_class <- c(RG = "#33A02C", IPC = "#377EB8", CPN = "#E41A1C")

vista_panel <- function(d, tests, acetylation = FALSE) {
  order <- if (acetylation) c("Ac16", "Ac3", "Ac5", "Ac6", "Ac7", "Ac9", "Ac14") else c("A1", "A6", "A9", "A12", "A17")
  stage <- if (acetylation) c("Early", "Early", "Intermediate", "Intermediate", "Intermediate", "Late", "Constitutive") else c("Early", "Intermediate", "Intermediate", "Late", "Constitutive")
  stopifnot(setequal(d$Module, order))
  stopifnot(all(abs(d$pct - 100 * d$n / d$class_total) < 1e-10))
  d$Position <- match(d$Module, order)
  d$VISTA_class <- factor(d$VISTA_class, levels = c("VZ", "non-VZ"))
  maximum <- if (acetylation) 45 else 90
  labels <- data.frame(Position = seq_along(order), Module = order)
  significant <- tests[!is.na(tests$BH_FDR) & tests$BH_FDR < 0.05, ]
  significant$Position <- match(significant$Module, order)
  significant$Height <- vapply(significant$Module, function(module) max(d$pct[d$Module == module]), numeric(1)) + maximum * 0.025
  totals <- vapply(c("VZ", "non-VZ"), function(group) {
    unique(d$class_total[d$VISTA_class == group])
  }, numeric(1))
  ggplot(d, aes(Position, pct, fill = VISTA_class)) +
    geom_col(position = position_dodge(width = 0.82), width = 0.78) +
    geom_rect(
      data = labels,
      aes(xmin = Position - 0.40, xmax = Position + 0.40,
          ymin = -maximum * 0.24, ymax = -maximum * 0.035),
      inherit.aes = FALSE,
      fill = unname(MS$category[stage]),
      color = NA
    ) +
    geom_text(
      data = labels, aes(Position, -maximum * 0.135, label = Module),
      inherit.aes = FALSE, angle = 90, size = 5, fontface = "bold"
    ) +
    geom_segment(
      data = significant,
      aes(x = Position - 0.38, xend = Position + 0.38, y = Height, yend = Height),
      inherit.aes = FALSE, linewidth = 0.35
    ) +
    geom_text(
      data = significant, aes(Position, Height + maximum * 0.04, label = "*"),
      inherit.aes = FALSE, size = 5, fontface = "bold"
    ) +
    geom_hline(yintercept = 0, linewidth = 0.45) +
    annotate("segment", x = 0.5, xend = 0.5, y = 0, yend = maximum, linewidth = 0.45) +
    scale_fill_manual(
      values = MS$vista,
      labels = c(sprintf("VZ enh. (n=%d)", totals[1]), sprintf("non-VZ enh. (n=%d)", totals[2]))
    ) +
    scale_x_continuous(breaks = NULL, expand = expansion(mult = 0)) +
    scale_y_continuous(breaks = if (acetylation) seq(0, 40, 10) else seq(0, 80, 20), expand = expansion(mult = 0)) +
    coord_cartesian(xlim = c(0.5, length(order) + 0.5), ylim = c(-maximum * 0.25, maximum * 1.04), clip = "off") +
    labs(x = NULL, y = "Proportion within class (%)", fill = NULL) +
    theme_classic(base_size = 11, base_family = "Arial") +
    theme(
      axis.line = element_blank(), axis.ticks.x = element_blank(),
      axis.text.y = element_text(size = 10, color = "black"),
      axis.title.y = element_text(size = 11, face = "bold", margin = margin(r = 2)),
      legend.position = "inside", legend.position.inside = c(0.99, 0.96),
      legend.justification = c(1, 1), legend.text = element_text(size = 9, face = "bold"),
      legend.key.size = grid::unit(2.5, "mm"), legend.spacing.y = grid::unit(0, "mm"),
      legend.margin = margin(0, 0, 0, 0), plot.margin = margin(5, 3, 0, 1)
    )
}

# Extended Data 1f/g use the module colors sampled from the final artwork's
# labeled legend swatches. Counts and module identities remain source values.
plot_ext1_differential_counts <- function(sheet) {
  d <- read_table("Extended_Fig1_KC", sheet)
  palette <- read.csv(file.path(REPO, "provenance/Extended_Fig1_module_colors.csv"))
  is_ac <- sheet == "1f_DiffBind_counts"
  groups <- if (is_ac) list(
    c(1:3, 16, 17), 4:8, 9:12, c(13:15, 20, 29), c(18, 19, 21:28)
  ) else list(
    c(1:3, 16, 17), 4:9, c(10:13, 19), c(14, 15, 20), c(18, 21:26, 28:33, 35)
  )
  prefix <- if (is_ac) "Ac" else "TMK"
  order <- paste0(prefix, unlist(groups))
  order <- order[order %in% unique(d$Module)]
  stopifnot(setequal(order, unique(d$Module)), all(order %in% palette$Module))
  colors <- setNames(palette$Color, palette$Module)
  d$Module <- factor(d$Module, levels = order)
  comparison <- unique(d$Comparison)
  d$Comparison <- factor(d$Comparison, levels = comparison)
  labels <- if (is_ac) c("VZ E10.5 vs VZ E12.5", "VZ E10.5 vs L5N P2") else
    c("VZ E11.5 vs VZ E13.5", "VZ E12.5 vs L6N P2")
  # The short labels preserve the actual source contrast identities.
  stopifnot(length(comparison) == 2)
  p <- ggplot(d, aes(Comparison, Significant_interval_count, fill = Module)) +
    geom_col(width = 0.9) +
    scale_fill_manual(values = colors, breaks = order) +
    scale_x_discrete(labels = setNames(labels, comparison)) +
    scale_y_continuous(breaks = if (is_ac) seq(0, 12500, 2500) else seq(0, 5000, 1000),
      expand = expansion(mult = c(0, 0.04))) +
    labs(title = paste0(if (is_ac) "H3K27ac" else "H3K27me3",
        " significant DiffBind intervals\n(counts) by module"),
      subtitle = "Threshold: FDR < 0.05", x = "Comparison",
      y = "Number of intervals", fill = "Module") +
    guides(fill = guide_legend(ncol = if (is_ac) 5 else 6, byrow = TRUE)) +
    theme_minimal(base_size = 12, base_family = "Arial") +
    theme(plot.title = element_text(size = 14, face = "bold", hjust = 0.5),
      plot.subtitle = element_text(size = 10), panel.grid.minor = element_blank(),
      axis.text.x = element_text(size = 9), legend.position = "bottom",
      legend.text = element_text(size = 8), legend.title = element_text(size = 9),
      legend.key.height = grid::unit(2.5, "mm"), legend.key.width = grid::unit(4, "mm"))
  built <- ggplot_build(p)$data[[1]]
  expected <- tapply(d$Significant_interval_count, d$Comparison, sum)
  plotted <- tapply(built$ymax - built$ymin, built$x, sum)
  stopifnot(max(abs(as.numeric(expected) - as.numeric(plotted))) < 1e-8)
  export_panel(p, paste0("Ext1/", sheet), 6.5, 6,
    note = "Original legend colors and category stacking; original absolute counts retained; source contrast ages preserved")
}
