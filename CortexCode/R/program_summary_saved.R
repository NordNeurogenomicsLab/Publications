# Original program pooling code from the final author notebook.
# Loads saved cluster assignments; does not refit models.
`%||%` <- function(a, b) if (!is.null(a)) a else b

# -----------------------------
# 0) Load bundles
# -----------------------------
out_dir <- file.path(REPO, "objects/integration")
bundles_file <- file.path(out_dir, "umap_reclustering_bundles.RData")
stopifnot(file.exists(bundles_file))
base::load(bundles_file)  # loads: bundles

# 3-way bundle
b3 <- bundles$bundle_threeway_hi_res %||% bundles$bundle_threeway_high_res
stopifnot(!is.null(b3), !is.null(b3$umap_df))
df_umap_3way <- b3$umap_df

# ATAC–AC bundle
bAA <- bundles$bundle_atac_ac_hi_res %||% bundles$bundle_atac_ac_high_res %||%
  bundles$bundle_atac_ac %||% bundles$bundle_atac_ac_hi %||%
  bundles$bundle_atac_ac_high
stopifnot(!is.null(bAA), !is.null(bAA$umap_df))
df_umap_atac_ac <- bAA$umap_df

# ATAC–ME3 LOW-RES bundle
bAM_low <- bundles$bundle_atac_me3_low_res %||%
  bundles$bundle_atac_me3_low %||%
  bundles$bundle_atac_me3_low_resolution %||%
  bundles$bundle_atac_me3_lo_res %||%
  bundles$bundle_atac_me3_lr
stopifnot(!is.null(bAM_low), !is.null(bAM_low$umap_df))
df_umap_atac_me3_low <- bAM_low$umap_df

# -----------------------------
# 1) Parameters
# -----------------------------
stages <- c("VZ", "SVZ", "ImN", "L6N", "L5N")
max_level <- 10
total_connections <- 157984L

modality_colors <- c(
  ATAC = "#545554",
  AC   = "#4eb748",
  ME3  = "#ed3b25"
)

# -----------------------------
# 2) Detect cluster columns
# -----------------------------
detect_cluster_col <- function(df, candidates) {
  hit <- candidates[candidates %in% names(df)]
  if (length(hit) == 0) NA_character_ else hit[[1]]
}

cluster_col_3way <- detect_cluster_col(
  df_umap_3way,
  c("three_way_cluster", "threeway_cluster", "cluster", "cluster_id")
)

cluster_col_aa <- detect_cluster_col(
  df_umap_atac_ac,
  c("two_way_cluster", "two_way_cluster_id", "twoWay_cluster", "cluster", "cluster_id")
)

cluster_col_am <- detect_cluster_col(
  df_umap_atac_me3_low,
  c("two_way_cluster", "two_way_cluster_id", "twoWay_cluster", "cluster", "cluster_id")
)

stopifnot(!is.na(cluster_col_3way), !is.na(cluster_col_aa), !is.na(cluster_col_am))

# -----------------------------
# 3) Robust vector parsing → long format
# -----------------------------
parse_vec_string <- function(x, stages) {
  if (is.list(x)) {
    mat <- base::do.call(rbind, lapply(x, function(v) as.numeric(v)))
    colnames(mat) <- stages
    return(mat)
  }

  if (is.numeric(x)) {
    mat <- matrix(as.numeric(x), ncol = length(stages), byrow = TRUE)
    colnames(mat) <- stages
    return(mat)
  }

  xs <- as.character(x)
  xs <- gsub("\\s+", "", xs)
  xs <- gsub(";", ",", xs, fixed = TRUE)
  xs <- gsub("-", ",", xs, fixed = TRUE)
  pieces <- strsplit(xs, ",", fixed = TRUE)
  nums <- unlist(lapply(pieces, as.numeric), use.names = FALSE)
  mat <- matrix(nums, ncol = length(stages), byrow = TRUE)
  colnames(mat) <- stages
  mat
}

vec_to_long_from_umap <- function(df_in, mod, stages) {
  col <- paste0(mod, "_vec")
  stopifnot(col %in% names(df_in))

  mat <- parse_vec_string(df_in[[col]], stages)

  base::as.data.frame(mat) %>%
    dplyr::mutate(row_id = dplyr::row_number()) %>%
    tidyr::pivot_longer(
      cols = -row_id,
      names_to = "Stage",
      values_to = "quantized_value"
    ) %>%
    dplyr::mutate(
      Modality = mod,
      Stage = factor(Stage, levels = stages)
    )
}

# -----------------------------
# 4) Program summarizers
# -----------------------------
summarize_program_3way <- function(df_umap, cluster_col, clusters, program_name,
                                   stages, max_level, total_connections) {

  sel <- as.character(clusters)

  df_sel <- df_umap %>%
    dplyr::mutate(cluster = as.character(.data[[cluster_col]])) %>%
    dplyr::filter(cluster %in% sel)

  n_loci <- nrow(df_sel)
  stopifnot(n_loci > 0)

  pct <- 100 * (n_loci / as.numeric(total_connections))

  long_df <- dplyr::bind_rows(
    vec_to_long_from_umap(df_sel, "ATAC", stages),
    vec_to_long_from_umap(df_sel, "AC",   stages),
    vec_to_long_from_umap(df_sel, "ME3",  stages)
  ) %>%
    dplyr::mutate(Modality = factor(Modality, levels = c("ATAC", "AC", "ME3")))

  long_df %>%
    dplyr::group_by(Modality, Stage) %>%
    dplyr::summarise(
      mean_quantized = mean(quantized_value, na.rm = TRUE),
      sd_quantized   = stats::sd(quantized_value, na.rm = TRUE),
      .groups = "drop"
    ) %>%
    dplyr::mutate(
      Program = paste0(program_name, "\n(n=", n_loci, " loci; ", sprintf("%.2f", pct), "% of all)"),
      ymin = pmax(0, mean_quantized - sd_quantized),
      ymax = pmin(max_level, mean_quantized + sd_quantized)
    )
}

summarize_program_pair <- function(df_umap, cluster_col, clusters, program_name,
                                   stages, max_level, total_connections,
                                   mod1, mod2) {

  sel <- as.character(clusters)

  df_sel <- df_umap %>%
    dplyr::mutate(cluster = as.character(.data[[cluster_col]])) %>%
    dplyr::filter(cluster %in% sel)

  n_loci <- nrow(df_sel)
  stopifnot(n_loci > 0)

  pct <- 100 * (n_loci / as.numeric(total_connections))

  long_df <- dplyr::bind_rows(
    vec_to_long_from_umap(df_sel, mod1, stages),
    vec_to_long_from_umap(df_sel, mod2, stages)
  ) %>%
    dplyr::mutate(Modality = factor(Modality, levels = c("ATAC", "AC", "ME3")))

  long_df %>%
    dplyr::group_by(Modality, Stage) %>%
    dplyr::summarise(
      mean_quantized = mean(quantized_value, na.rm = TRUE),
      sd_quantized   = stats::sd(quantized_value, na.rm = TRUE),
      .groups = "drop"
    ) %>%
    dplyr::mutate(
      Program = paste0(program_name, "\n(n=", n_loci, " loci; ", sprintf("%.2f", pct), "% of all)"),
      ymin = pmax(0, mean_quantized - sd_quantized),
      ymax = pmin(max_level, mean_quantized + sd_quantized)
    )
}

# -----------------------------
# 5) Program definitions
# NOTE:
# Cluster numbers are bundle-specific. For example, cluster 4 in the
# 3-way bundle is unrelated to cluster 4 in the ATAC-ME3 bundle.
# -----------------------------
program_defs_3way <- list(
  "1a. VZ/SVZ decommissioning with repression" = c(6, 16, 28),
  "1d. Activation and derepression of SVZ pREs" = c(11, 26, 36),
  "1f. Repression dominant" = c(4),

  "2a. Activation with derepression" = c(8, 21, 34, 39),
  "2b. Derepression-led activation" = c(20, 29, 35),
  "3a. Progressive activation with derepression" = c(1, 3, 12, 18, 27),
  "3b. Repression of embryonic cell fates" = c(7, 17, 24, 25),
  "4a. Constitutively active" = c(10, 19, 37, 38),
  "4b. Constitutively active / bivalent" = c(5, 13, 32, 33)
)

program_defs_atac_ac <- list(
  "1b. H3K27me3-independent decommissioning" = c(1, 3, 20, 25, 26, 28, 29, 32, 33, 37, 39),
  "1e. Transient activation of SVZ pREs (H3K27me3 independent)" = c(2, 5, 10, 35),
  "2c. Progressive activation (H3K27me3-independent)" = c(9, 11, 12, 13, 16, 22, 30, 44),
  "3c. Progressive activation (H3K27me3-independent)" = c(14, 15, 24, 34, 40, 42, 43),
  "4c. ATAC-H3K27ac constitutive" = c(17, 18, 21, 36, 38, 41)
)

program_defs_atac_me3_low <- list(
  "1c. Progressive repression" = c(3, 4),
  "2d. ATAC/H3K27me3 derepression" = c(6),
  "4d. Progressive derepression (non-cortical brain pREs?)" = c(2)
)

# -----------------------------
# 6) Build summary table
# -----------------------------
prog_3way_df <- dplyr::bind_rows(lapply(names(program_defs_3way), function(nm) {
  summarize_program_3way(
    df_umap_3way,
    cluster_col_3way,
    program_defs_3way[[nm]],
    nm,
    stages = stages,
    max_level = max_level,
    total_connections = total_connections
  )
}))

prog_atac_ac_df <- dplyr::bind_rows(lapply(names(program_defs_atac_ac), function(nm) {
  summarize_program_pair(
    df_umap_atac_ac,
    cluster_col_aa,
    program_defs_atac_ac[[nm]],
    nm,
    stages = stages,
    max_level = max_level,
    total_connections = total_connections,
    mod1 = "ATAC",
    mod2 = "AC"
  )
}))

prog_atac_me3_df <- dplyr::bind_rows(lapply(names(program_defs_atac_me3_low), function(nm) {
  summarize_program_pair(
    df_umap_atac_me3_low,
    cluster_col_am,
    program_defs_atac_me3_low[[nm]],
    nm,
    stages = stages,
    max_level = max_level,
    total_connections = total_connections,
    mod1 = "ATAC",
    mod2 = "ME3"
  )
}))

prog_df <- dplyr::bind_rows(prog_3way_df, prog_atac_ac_df, prog_atac_me3_df)
