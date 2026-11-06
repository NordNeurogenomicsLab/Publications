# Recorded developmental-stage aggregation and ordinal encoder.
# Source: Kdm6b_ME3_AC_and_ATAC_big_picture_Fg.6_v4.Rmd.
construct_metadata <- function(df, modality) {
  # Pick suffix by modality
  suffix <- ifelse(modality == "ATAC", "_cpm$", "_fpkm$")

  # Extract sample names
  sample_names <- grep(
    suffix, names(df),
    value = TRUE
  )
  sample_clean <- gsub(suffix, "", sample_names)
  sample_clean <- gsub("^(Atac_|Hac_|Htk_)", "", sample_clean)  # strip prefixes

  metadata <- data.frame(
    RawName = sample_names,
    Samples = sample_clean,
    stringsAsFactors = FALSE
  )

  # Age: extract E## or P##
  metadata$age <- gsub(".*(E\\d+|P\\d+).*", "\\1", metadata$Samples)

  # Region: _C = Caudal, _R = Rostral
  metadata$region <- ifelse(
    grepl("_C", metadata$Samples),
    "Caudal", ifelse(
      grepl("_R", metadata$Samples),
      "Rostral", NA
    )
  )

  # Developmental compartment
  metadata$dev <- NA
  metadata$dev[grepl("Nex", metadata$Samples)] <- "Neurons"
  metadata$dev[grepl("Vz", metadata$Samples)] <- "VZ"
  metadata$dev[grepl("Svz", metadata$Samples)] <- "SVZ"
  metadata$dev[grepl("Ntsr", metadata$Samples)] <- "L6"
  metadata$dev[grepl("Rbp", metadata$Samples)] <- "L5"

  return(metadata)
}
# Usage
build_stage_matrix <- function(df, metadata, suffix) {

  value_cols <- grep(
    paste0(suffix, "$"),
    names(df),
    value = TRUE
  )

  df %>%
    dplyr::select(ID, all_of(value_cols)) %>%
    tidyr::pivot_longer(
      cols = all_of(value_cols),
      names_to = "RawName",
      values_to = "value"
    ) %>%
    dplyr::left_join(metadata, by = "RawName") %>%
    dplyr::filter(dev %in% c("VZ", "SVZ", "Neurons", "L6", "L5")) %>%
    dplyr::mutate(
      stage = dplyr::recode(
        dev,
        VZ = "VZ",
        SVZ = "SVZ",
        Neurons = "ImN",
        L6 = "L6N",

        L5 = "L5N"
      )
    ) %>%
    dplyr::group_by(ID, stage) %>%
    dplyr::summarise(
      value = median(value, na.rm = TRUE),

      .groups = "drop"
    ) %>%
    tidyr::pivot_wider(names_from = stage,
      values_from = value) %>%
    dplyr::arrange(ID)
}
ordinalize <- function(mat,
  stages,
  zero_q = 0.05,
  n_levels = 10,
  id_col = "ID") {

  x <- as.matrix(mat[, stages, drop = FALSE])
  ids <- mat[[id_col]]

  zero_cut <- as.numeric(
    quantile(
      as.numeric(x),
      zero_q,
      na.rm = TRUE,
      type = 7
    )
  )

  nz <- as.numeric(x)[as.numeric(x) >
    zero_cut]
  q <- as.numeric(
    quantile(
      nz,
      probs = seq(0, 1, length.out = n_levels + 1),

      na.rm = TRUE,
      type = 7
    )
  )

  # guard ties
  q <- unique(q)
  if ((length(q) -
    1) < n_levels) {
    warning("Quantile breaks collapsed; falling back to rank-binning.")
    out <- matrix(
      0L, nrow = nrow(x),
      ncol = ncol(x),
      dimnames = list(ids, stages)
    )
    keep <- x > zero_cut
    rk <- rank(x[keep], ties.method = "average")
    out[keep] <- pmax(
      1L, pmin(
        n_levels, floor(
          (rk - 1)/length(rk) *
          n_levels
        ) +
          1L
      )
    )
    return(out)
  }

  encode <- function(v) {
    out <- integer(length(v))
    keep <- v > zero_cut
    out[!keep] <- 0L
    out[keep] <- as.integer(cut(v[keep],
      breaks = q,
      labels = 1:n_levels,
      include.lowest = TRUE))
    out
  }

  out <- apply(x, 2, encode)
  rownames(out) <- ids
  storage.mode(out) <- "integer"
  out
}
