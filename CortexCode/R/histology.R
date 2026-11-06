# Parsers preserve original Excel row numbers and distinguish sections
# from mice.
num <- function(x) suppressWarnings(as.numeric(x))
valid <- function(x) !is.na(x) &&
  nzchar(x)
isgeno <- function(x) valid(x) &&
  x %in% c("WT", "KO", "Mut", "het")
plot_mice <- function(d,
  name,
  ylab,
  width = 7,
  height = 4) {
  stopifnot(
    all(
      c("Mouse", "Genotype", "Region", "Value") %in%
        names(d)
    ),
    all(is.finite(d$Value))
  )
  d$Region <- factor(d$Region, levels = unique(d$Region))
  m <- d %>%
    group_by(Mouse, Genotype, Region) %>%
    summarise(
      Value = mean(Value),
      Sections = n(),
      .groups = "drop"
    )
  m$Genotype <- factor(m$Genotype, levels = c("WT", "het", "KO", "Mut"))
  cols <- if (grepl("Fig7/|Ext9/", name))
    MS$deletion else MS$histology
  if (name == "Ext7/7j_CP_micrometres")
    cols <- MS$histology_ext7_j
  one_region <- length(unique(m$Region)) == 1
  p <- ggplot(m, aes(if (one_region) Genotype else Region, Value, fill = Genotype)) +
    stat_summary(
      fun = mean,
      geom = "col",
      position = position_dodge(0.8),

      width = 0.7,
      color = "black",
      linewidth = 0.3
    ) +
    stat_summary(
      fun.data = mean_se,
      geom = "errorbar",
      position = position_dodge(0.8),

      width = 0.15,
      linewidth = 0.3
    ) +
    geom_point(
      aes(group = Genotype),
      position = position_jitterdodge(jitter.width = 0.06, dodge.width = 0.8, seed = 1),

      size = 1.5,
      color = "black"
    ) +
    scale_fill_manual(values = cols) +
    labs(x = NULL,
      y = ylab)
  if (one_region) {
    p <- p + scale_x_discrete(labels = c(WT = "WT", KO = "cKO", Mut = "cKO", het = "Het")) +
      theme(legend.position = "none")
  }
  export_panel(
    p, name, width, height, note = "Each point is a biological mouse mean; original section/hemisphere measurements retained separately"
  )
  f <- file.path(REPO, "reports", gsub("/", "_", name, fixed = TRUE))
  write.csv(
    d, paste0(f, "_measurements.csv"),
    row.names = FALSE
  )
  write.csv(
    m, paste0(f, "_mouse_means.csv"),
    row.names = FALSE
  )
  save_object(
    list(measurements = d, mouse_means = m),

    paste0("histology/", basename(name))
  )
  invisible(m)
}


# Mouse IDs and genotype labels occupy different rows in this original grid.
# Resolve genotype within each complete mouse block, including its first section.
parse_ext7_tle4_measurements <- function() {
  x <- read_raw("Extended_Fig7_KC", "Ext_7j Tle4_thickness")
  stopifnot(x[1, 1] == "1mm")
  calibration <- num(x[1, 2])
  starts <- which(!is.na(x[, 1]) & grepl("^[0-9]+-[0-9]+$", x[, 1]))
  stopifnot(length(starts) == 6, is.finite(calibration), calibration > 0)
  rows <- list()
  for (block in seq_along(starts)) {
    first <- starts[block]
    last <- if (block < length(starts)) starts[block + 1] - 1 else nrow(x)
    index <- seq.int(first, last)
    genotype <- unique(x[index, 1][vapply(x[index, 1], isgeno, logical(1))])
    stopifnot(length(genotype) == 1)
    measurement_rows <- index[
      !is.na(x[index, 2]) & grepl("^section", x[index, 2]) &
        is.finite(suppressWarnings(as.numeric(x[index, 5])))
    ]
    stopifnot(length(measurement_rows) == 4)
    for (i in measurement_rows) {
      percent <- num(x[i, 5])
      stopifnot(abs(percent - 100 * num(x[i, 3]) / num(x[i, 4])) < 1e-8)
      rows[[length(rows) + 1]] <- data.frame(
        Mouse = x[first, 1],
        Genotype = genotype,
        Section = x[i, 2],
        Excel_row = i,
        Tle4_micrometres = 1000 * num(x[i, 3]) / calibration,
        CP_micrometres = 1000 * num(x[i, 4]) / calibration,
        Tle4_percent = percent
      )
    }
    # Independently verify the average lengths against the supplied um row.
    um <- index[!is.na(x[index, 2]) & x[index, 2] == "um"]
    stopifnot(length(um) == 1)
    marker_mean <- mean(vapply(measurement_rows, function(i) num(x[i, 3]), numeric(1))) * 1000 / calibration
    cortex_mean <- mean(vapply(measurement_rows, function(i) num(x[i, 4]), numeric(1))) * 1000 / calibration
    stopifnot(abs(marker_mean - num(x[um, 3])) < 1e-8)
    stopifnot(abs(cortex_mean - num(x[um, 4])) < 1e-8)
  }
  result <- bind_rows(rows)
  mice <- unique(result[, c("Mouse", "Genotype")])
  stopifnot(nrow(result) == 24, sum(mice$Genotype == "WT") == 3, sum(mice$Genotype == "KO") == 3)
  result
}

histology_ext7_lengths <- function() {
  measured <- parse_ext7_tle4_measurements()
  for (column in c("CP_micrometres", "Tle4_micrometres")) {
    d <- measured
    d$Region <- if (column == "CP_micrometres") "Cortical plate" else "TLE4-positive layer"
    d$Value <- d[[column]]
    name <- if (column == "CP_micrometres") "Ext7/7j_CP_micrometres" else "Ext7/7j_TLE4_micrometres"
    plot_mice(d, name, "Thickness (micrometres)", 3, 3.5)
  }
  invisible(measured)
}

histology_fig5 <- function() {
  # 5e: normalized TLE4 signal in four cortical positions. Column 1
  # alternates IDs and genotypes.
  x <- read_raw("Fig5_KC", "5e_TLE4_E14")
  out <- list()
  mouse <- geno <- NA_character_
  for (i in seq_len(nrow(x))) {
    a <- x[i, 1]
    if (isgeno(a))
      geno <- a else if (valid(a))
      mouse <- a
    if (valid(x[i, 2]) &&
      grepl("^section", x[i, 2]) &&
      is.finite(num(x[i, 4])))
      for (j in 4:7) out[[length(out) +
        1]] <- data.frame(
        Mouse = mouse,
        Genotype = geno,
        Section = x[i, 2],
        Excel_row = i,

        Region = paste0("Cx", j - 3),
        Value = num(x[i, j])
      )
  }
  plot_mice(
    bind_rows(out),
    "Fig5/5e_TLE4_E14", "TLE4 signal (% LGE)"
  )
  # 5f: paired left/right measurements; metadata are carried down only
  # over numeric rows.
  x <- read_raw("Fig5_KC", "5f_Tbr1_E14")
  out <- list()
  mouse <- geno <- section <- NA_character_
  for (i in seq_len(nrow(x))) {
    if (isgeno(x[i, 1]))
      geno <- x[i, 1]
    if (valid(x[i, 2]))
      mouse <- x[i, 2]
    if (valid(x[i, 3]))
      section <- x[i, 3]
    if (is.finite(num(x[i, 4])) &&
      is.finite(num(x[i, 7])))
      for (j in c(5, 6, 8, 9)) out[[length(out) +
        1]] <- data.frame(
        Mouse = mouse,
        Genotype = geno,
        Section = section,
        Excel_row = i,

        Region = if (j %in% c(5, 8))
          "CP" else "IZ", Hemisphere = if (j < 7)
          "Left" else "Right", Value = num(x[i, j])
      )
  }
  plot_mice(
    bind_rows(out),
    "Fig5/5f_Tbr1_E14", "TBR1 signal (% pTE)", 5, 4
  )
  for (pair in list(
    c("5g_Tbr1_P10", "5", "TBR1 thickness (um)"),
    c("5h_Bcl11b_P10", "5", "BCL11B thickness (um)"),
    c("5i_Fezf2 P10", "6", "FEZF2 thickness (um)")
  )) {
    x <- read_raw("Fig5_KC", pair[1])
    out <- list()
    mouse <- geno <- NA_character_
    for (i in seq_len(nrow(x))) {
      if (isgeno(x[i, 1]))
        geno <- x[i, 1] else if (valid(x[i, 1]))
        mouse <- x[i, 1]
      if (is.finite(num(x[i, 2])) &&
        is.finite(num(x[i, as.integer(pair[2])])))
        out[[length(out) +
          1]] <- data.frame(
          Mouse = mouse,
          Genotype = geno,
          Section = x[i, 2],
          Excel_row = i,

          Region = "Dorsal cortex",
          Value = num(x[i, as.integer(pair[2])])
        )
    }
    d <- bind_rows(out)
    plot_mice(
      d, paste0("Fig5/", gsub(" ", "_", pair[1])),
      pair[3], 4, 4
    )
  }
}
histology_ext7 <- function() {
  measured <- parse_ext7_tle4_measurements()
  measured$Region <- "Dorsal cortex"
  measured$Value <- measured$Tle4_percent
  plot_mice(
    measured,
    "Ext7/7j_TLE4_thickness",
    "TLE4 thickness (% cortical thickness)",
    5,
    4
  )
  x <- read_raw("Extended_Fig7_KC", "Ext_7kTbr1_thickness")
  out <- list()
  mouse <- geno <- section <- NA_character_
  for (i in seq_len(nrow(x))) {
    if (isgeno(x[i, 1]))
      geno <- x[i, 1]
    if (valid(x[i, 2]))
      mouse <- x[i, 2]
    if (valid(x[i, 3]))
      section <- x[i, 3]
    if (is.finite(num(x[i, 4])) &&
      is.finite(num(x[i, 7])))
      for (j in c(5, 6, 8, 9)) out[[length(out) +
        1]] <- data.frame(
        Mouse = mouse,
        Genotype = geno,
        Section = section,
        Excel_row = i,

        Region = if (j %in% c(5, 8))
          "MZ" else "IZ", Hemisphere = if (j < 7)
          "Left" else "Right", Value = num(x[i, j])
      )
  }
  layer_mice <- plot_mice(
    bind_rows(out),
    "Ext7/7k_TBR1_layer_thickness", "Layer thickness (% CP)", 5, 4
  )
  for (region in c("MZ", "IZ")) {
    z <- layer_mice[layer_mice$Region == region, ]
    p <- ggplot(z, aes(Genotype, Value, fill = Genotype)) +
      stat_summary(fun = mean, geom = "col", width = 0.65, color = "black", linewidth = 0.3) +
      stat_summary(fun.data = mean_se, geom = "errorbar", width = 0.15, linewidth = 0.3) +
      geom_point(position = position_jitter(width = 0.06, seed = 1), size = 1.3) +
      scale_fill_manual(values = if (region == "MZ") MS$histology_ext7_j else MS$histology) +
      labs(title = paste("TBR1-positive", region), x = NULL, y = "Thickness (% CP)") +
      theme(legend.position = "none")
    export_panel(p, paste0("Ext7/7k_TBR1_", region), 3, 3.5)
  }
}
histology_fig7 <- function() {
  for (pair in list(
    c("Fig7_KC", "7k_l_Tbr1_416Del", "Fig7/7kl_Tbr1_416Del"),
    c("Extended_Fig9_KC", "9f_Tbr1_Nkx2-1", "Ext9/9f_Tbr1_Nkx2-1")
  )) {
    x <- read_raw(pair[1], pair[2])
    out <- list()
    geno <- NA_character_
    for (i in seq_len(nrow(x))) {
      if (isgeno(x[i, 1]))
        geno <- x[i, 1]
      if (is.finite(num(x[i, 4])) &&
        valid(x[i, 2]))
        for (j in 11:14) {
          v <- num(x[i, j])
          expected <- 100 * num(x[i, j - 6])/num(x[i, 4])
          stopifnot(
          is.finite(v),
          abs(v - expected) <
            1e-08
        )
          out[[length(out) +
          1]] <- data.frame(
          Mouse = x[i, 2],
          Genotype = geno,
          Section = x[i, 3],

          Excel_row = i,
          Region = c("Lateral CP", "Dorsal CP", "IZ", "Dorsomedial CP")[j -
            10], Value = v
        )
        }
    }
    d <- bind_rows(out)
    n <- unique(d[, c("Mouse", "Genotype")])
    stopifnot(
      sum(n$Genotype == "WT") ==
        1, sum(n$Genotype == "het") ==
        3, sum(n$Genotype == "KO") == 4
    )
    m <- plot_mice(d, pair[3], "TBR1 signal (% Nkx signal)", 7, 4)
    if (pair[1] == "Fig7_KC")
      for (region in c("Dorsal CP", "Dorsomedial CP")) {
        z <- m[m$Region == region, ]
        p <- ggplot(z, aes(Genotype, Value, fill = Genotype)) +
          stat_summary(
          fun = mean,
          geom = "col",
          width = 0.65,
          color = "black",

          linewidth = 0.4
        ) +
          stat_summary(
          fun.data = mean_se,
          geom = "errorbar",
          width = 0.15,

          linewidth = 0.4
        ) +
          geom_point(
          position = position_jitter(width = 0.06, seed = 1),

          size = 1.5
        ) +
          scale_fill_manual(values = MS$deletion) +
          labs(x = NULL,
            y = "Normalized TBR1 ISH intensity",
            title = region) +
          theme(legend.position = "none")
        export_panel(
          p, paste0(
          "Fig7/", if (region == "Dorsal CP")
            "7k_Dorsal_CP" else "7l_Dorsomedial_CP"
        ),
          3.5,
          4,
          note = "Biological mouse means and individual animals; WT n=1, het n=3, KO n=4"
        )
      }
  }
}
