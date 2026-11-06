# Saved numerical layers only. No peak calling, differential tests, or LinkPeaks.
prepare_extended5_linked_scores <- function(cache, link_directory) {
  mat <- cache$SCT
  metadata <- cache$metadata
  scores <- list()
  audit <- list()
  for (assay in c("AC", "ME3")) {
    links <- read.csv(file.path(link_directory, paste0(assay, ".linked_peaks.csv")))
    key <- read_table("Fig1_KC", if (assay == "AC") "d3_Peak_annotations" else "e3_Peak_annotations")
    key$Feature <- gsub("_", "-", key$ID, fixed = TRUE)
    stopifnot(all(links$query_region %in% key$Feature))
    # Nic's Figure 3 procedure: positive links, lowest P value per query peak,
    # then unique genes per module. Apply it to the saved AC/ME3 links.
    selected <- links %>% filter(zscore >= 0) %>%
      group_by(query_region) %>%
      slice_min(pvalue, n = 1, with_ties = FALSE) %>% ungroup()
    selected$Module <- key$Fine_clusters_short[match(selected$query_region, key$Feature)]
    modules <- paste0(if (assay == "AC") "Ac" else "TMK", 1:15)
    sets <- lapply(modules, function(m) unique(selected$linked_gene[selected$Module == m]))
    names(sets) <- modules
    # Equivalent to AddModuleScore(..., nbin=24, ctrl=100, seed=1),
    # resetting the seed for each single-module call as in Nic's notebook.
    mean_expression <- Matrix::rowMeans(mat)
    mean_expression <- mean_expression[order(mean_expression)]
    set.seed(1)
    bins <- ggplot2::cut_number(mean_expression + rnorm(length(mean_expression)) / 1e30,
      n = 24, labels = FALSE, right = FALSE)
    names(bins) <- names(mean_expression)
    rng_after_bins <- .Random.seed
    raw <- matrix(NA_real_, 15, ncol(mat), dimnames = list(modules, colnames(mat)))
    controls <- list()
    for (module in modules) {
      features <- intersect(sets[[module]], rownames(mat))
      stopifnot(length(features) > 0)
      assign(".Random.seed", rng_after_bins, envir = .GlobalEnv)
      control <- unique(unlist(lapply(features, function(gene) {
        names(sample(bins[which(bins == bins[gene])], size = 100, replace = FALSE))
      }), use.names = FALSE))
      raw[module, ] <- Matrix::colMeans(mat[features, , drop = FALSE]) -
        Matrix::colMeans(mat[control, , drop = FALSE])
      controls[[module]] <- control
      audit[[paste(assay, module)]] <- data.frame(Assay = assay, Module = module,
        Linked_genes = length(sets[[module]]), SCT_genes = length(features), Controls = length(control))
    }
    scaled <- t(scale(t(raw)))
    scaled[scaled > 10] <- 10
    scores[[assay]] <- list(raw = raw, scaled = scaled,
      genes = sets, controls = controls, selected_links = selected)
    cat(assay, "linked-gene scores exported\n"); flush.console()
  }
  scores$metadata <- metadata
  scores$audit <- bind_rows(audit)
  scores
}

# The supplemental marker display selection is explicit because Nic's supplied
# final-Figure-3 notebook does not specify the supplemental top-marker count.
select_extended5_markers <- function(panel, available, top_n = 10) {
  d <- read_table("Extended_Fig5_KC", paste0(panel, "_Marker_results"))
  d <- d[d$p_val_adj < 0.05 & d$avg_log2FC > 0, ]
  d$cluster <- factor(d$cluster, levels = clusters)
  d <- d[order(d$cluster, -d$avg_log2FC, d$gene), ]
  d <- d %>% group_by(cluster) %>% slice_head(n = top_n) %>% ungroup()
  d$Display_order <- seq_len(nrow(d))
  d$Available_in_saved_layer <- d$gene %in% available
  d
}

# Heatmaps use the stored cell order within each Cell_Type_3 block. They do
# not recluster cells or recalculate marker P values.
plot_extended5_cell_heatmap <- function(mat, metadata, name,
  expression = TRUE, row_labels = FALSE) {
  metadata$Cell_Type_3 <- factor(metadata$Cell_Type_3, levels = clusters)
  cells <- rownames(metadata)[order(metadata$Cell_Type_3, seq_len(nrow(metadata)))]
  stopifnot(setequal(colnames(mat), cells), all(is.finite(mat)))
  mat <- mat[, cells, drop = FALSE]
  cluster <- metadata[cells, "Cell_Type_3"]
  cluster_key <- read.csv(file.path(REPO, "provenance/Extended5_cluster_colors.csv"))
  colors <- setNames(cluster_key$Color, cluster_key$Cluster)[clusters]
  header <- ComplexHeatmap::HeatmapAnnotation(
    Label = ComplexHeatmap::anno_block(labels = gsub("_", "", clusters),
      labels_rot = 90, labels_gp = grid::gpar(fontsize = 8),
      gp = grid::gpar(fill = "white", col = NA), height = grid::unit(8, "mm")),
    Cluster = ComplexHeatmap::anno_block(gp = grid::gpar(fill = colors, col = NA),
      height = grid::unit(1.5, "mm")),
    gap = grid::unit(0.5, "mm"), show_annotation_name = FALSE)
  col <- if (expression && !row_labels) {
    palette <- read.csv(file.path(REPO, "provenance/Extended5_RNA_marker_palette.csv"))
    circlize::colorRamp2(palette$Value, palette$Color, space = "RGB")
  } else if (expression) circlize::colorRamp2(c(-2.5, 0, 2.5),
      c("#FF00FF", "#000000", "#FFFF00"), space = "RGB") else
    { palette <- read.csv(file.path(REPO, "provenance/Extended5_marker_palette.csv"))
      circlize::colorRamp2(palette$Value, palette$Color, space = "RGB") }
  h <- ComplexHeatmap::Heatmap(mat, name = "Expression", col = col,
    cluster_rows = FALSE, cluster_columns = FALSE, cluster_column_slices = FALSE,
    column_split = factor(cluster, levels = clusters), column_gap = grid::unit(0.25, "mm"),
    column_title = NULL, top_annotation = header,
    show_column_names = FALSE, show_row_names = row_labels,
    row_names_side = "left", row_names_gp = grid::gpar(fontsize = 9),
    use_raster = TRUE, raster_quality = 3, raster_resize_mat = FALSE, raster_by_magick = FALSE,
    heatmap_legend_param = list(at = if (expression) c(-2, -1, 0, 1, 2) else c(-1, 0, 1, 2),
      title_gp = grid::gpar(fontsize = 8), labels_gp = grid::gpar(fontsize = 8)))
  export_complex_heatmap(h, name, if (row_labels) 6.8 else 6.2, if (row_labels) 2.7 else 3.0)
  save_object(list(matrix = mat, Cell_ID = cells, cluster = cluster,
    row_order = rownames(mat), display_min = if (expression) -2.5 else -1.5, display_max = 2.5), paste0("multiome/", basename(name)))
}

plot_extended_atac_mean_heatmap <- function(sheet, groups, name) {
  d <- read_table(if (grepl("^5", sheet)) "Extended_Fig5_KC" else "Extended_Fig6_KC", sheet)
  # Match manuscript equal-height module strips. Every stored row is retained.
  pieces <- lapply(groups, function(module) {
    x <- if (module == "All") d else d[d$Fine_clusters_short == module, ]
    m <- as.matrix(x[, clusters]);rownames(m) <- x$ID
    z <- t(scale(t(m)));z[!is.finite(z)] <- 0
    # Multi-module strips retain their saved source ordering. Only the
    # single TM3-linked block uses a small display-only dendrogram.
    hc <- if (length(groups) == 1 && nrow(z) <= 3500) hclust(dist(z), method = "complete") else NULL
    ord <- if (!is.null(hc)) {
      # Small display-only row clustering. No WGCNA or embedding is refitted.
      hc$order
    } else seq_len(nrow(z))
    list(matrix = if (is.null(hc)) z[ord, , drop = FALSE] else z,
      Row_tree = if (is.null(hc)) FALSE else as.dendrogram(hc),
      Display_order = rownames(z)[ord], Module = module)
  })
  palette_key <- read.csv(file.path(REPO, "provenance/Extended_ATAC_mean_palette.csv"))
  palette <- circlize::colorRamp2(palette_key$Value, palette_key$Color, space = "RGB")
  plots <- lapply(pieces, function(x) {
    strip_range <- if (length(groups) > 1 && x$Module == "A10") c(-2, 4) else c(-4, 4)
    strip_palette <- if (length(groups) > 1) circlize::colorRamp2(
      seq(strip_range[1], strip_range[2], length.out = nrow(palette_key)),
      palette_key$Color, space = "RGB") else palette
    ComplexHeatmap::Heatmap(x$matrix, name = paste0(x$Module, " z"), col = strip_palette,
      cluster_rows = x$Row_tree, row_dend_reorder = FALSE, show_row_dend = x$Module == "All", cluster_columns = FALSE,
      column_title = if (x$Module == "All") "Multiome ATAC accessibility\nacross pREs in TM3" else NULL,
      column_title_gp = grid::gpar(fontsize = 11, fontface = "bold"),
      show_row_names = FALSE, show_column_names = identical(x$Module, tail(groups, 1)),
      column_labels = gsub("_", "", clusters), column_names_rot = 90,
      column_names_gp = grid::gpar(fontsize = 9), row_title = if (x$Module == "All") NULL else x$Module,
      row_title_gp = grid::gpar(fontsize = 11),
      show_heatmap_legend = TRUE, heatmap_legend_param = list(at = if (length(groups) > 1) seq(strip_range[1], strip_range[2], by = 2) else c(-4, -2, 0, 2, 4), title = "Row z score"), use_raster = TRUE, raster_quality = 3,
      raster_resize_mat = FALSE, raster_by_magick = FALSE)
  })
  # Draw independent heatmaps in equal-height viewports, as assembled in the artwork.
  stem <- file.path(REPO, "panels", name)
  draw <- function() {
    grid::grid.newpage()
    n <- length(plots)
    for (i in seq_along(plots)) {
      grid::pushViewport(grid::viewport(x = 0.49, y = 1 - (i - 0.5) / n,
        width = 0.94, height = 0.98 / n))
      ComplexHeatmap::draw(plots[[i]], newpage = FALSE,
        padding = grid::unit(c(1, 1, 1, 1), "mm"))
      grid::upViewport()
    }
  }
  grDevices::cairo_pdf(paste0(stem, ".pdf"), width = if (length(groups) > 1) 6.2 else 4.3, height = if (length(groups) > 1) 6.8 else 5.3)
  draw();grDevices::dev.off()
  grDevices::png(paste0(stem, ".png"), width = if (length(groups) > 1) 6.2 else 4.3, height = if (length(groups) > 1) 6.8 else 5.3,
    units = "in", res = 180, type = "cairo")
  draw();grDevices::dev.off()
  save_object(list(blocks = pieces, cluster_order = clusters,
    scaling = "row z score across 16 stored cluster means",
    row_order = "single-block display uses complete-linkage Euclidean clustering when <=3500 rows; multi-strip displays retain saved input row order",
    display_limits = c(-4, 4)), paste0("multiome/", basename(name)))
  invisible(pieces)
}
