# Portable entry point. Run Rscript run_all.R from the repository root.
repo_root <- function() {
  p <- normalizePath(getwd(),
    winslash = "/",
    mustWork = TRUE)
  repeat {
    if (file.exists(file.path(p, "R", "common.R")) &&
      dir.exists(file.path(p, "provenance")))
      return(p)
    q <- dirname(p)
    if (identical(p, q))
      stop("Run from the GitHub_repo directory or one of its subdirectories.")
    p <- q
  }
}
REPO <- repo_root()
source(file.path(REPO, "R/input_catalog.R"))
suppressPackageStartupMessages(
  {
    library(ggplot2)
    library(dplyr)
    library(tidyr)
    library(jsonlite)
  }
)
options(stringsAsFactors = FALSE)
set.seed(1)
clusters <- c(
  paste0("RG_", 1:8),
  paste0("IPC_", 1:4),
  paste0("CPN_", 1:4)
)
read_table <- function(book, sheet) {
  f <- file.path(
    REPO, "data", "source_tables", book, paste0(
      gsub(" ", "_", sheet, fixed = TRUE),
      ".csv"
    )
  )
  if (!file.exists(f))
    stop("Missing bundled plotting input: ", f)
  data <- read.csv(
    f,
    check.names = FALSE,
    na.strings = c("", "NA"),

    stringsAsFactors = FALSE
  )
  if (sheet %in% c("6bprime_RPKM", "6o_RPKM", "7cde_jk_RPKM")) {
    metadata <- read.csv(file.path(REPO, "data/upstream/RNA/sample_metadata.csv"))
    index <- match(data$Sample, metadata$Sample)
    stopifnot(!anyNA(index), !anyDuplicated(metadata$Sample))
    for (field in c("Genotype", "Cells", "FlashTag_population", "Litter", "Sex", "DE_population")) {
      stopifnot(identical(as.character(data[[field]]), as.character(metadata[[field]][index])))
    }
  }
  data
}
read_raw <- function(book, sheet) read.csv(
  file.path(
    REPO, "data", "raw_sheets", book, paste0(
      gsub(" ", "_", sheet, fixed = TRUE),
      ".csv"
    )
  ),
  header = FALSE,
  check.names = FALSE,
  colClasses = "character",
  na.strings = "",

  fill = TRUE
)
save_object <- function(x, name) {
  p <- file.path(REPO, "objects", paste0(name, ".rds"))
  dir.create(
    dirname(p),
    recursive = TRUE,
    showWarnings = FALSE
  )
  saveRDS(x, p, compress = "gzip")
  invisible(p)
}
panel_registry <- list()
export_panel <- function(
  p,
  name,
  width = 5,
  height = 4,
  note = "Recreated from validated numerical inputs"
) {
  f <- file.path(REPO, "panels", name)
  dir.create(
    dirname(f),
    showWarnings = FALSE,
    recursive = TRUE
  )
  ggsave(
    paste0(f, ".pdf"),
    plot = p,
    width = width,
    height = height,
    device = grDevices::cairo_pdf,

    limitsize = FALSE
  )
  ggsave(
    paste0(f, ".png"),
    plot = p,
    width = width,
    height = height,
    dpi = 180,
    bg = "white",

    limitsize = FALSE
  )
  panel_registry[[name]] <<- data.frame(panel = name,
    note = note,
    width = width,
    height = height)
  invisible(p)
}
finish_notebook <- function(name) {
  write_session_info(name, REPO)
  dir.create(
    file.path(REPO, "reports"),
    showWarnings = FALSE
  )
  if (length(panel_registry))
    write.csv(
      bind_rows(panel_registry),
      file.path(REPO, "reports", paste0(name, "_panels.csv")),
      row.names = FALSE
    )
  writeLines(
    capture.output(sessionInfo()),
    file.path(REPO, "reports", paste0(name, "_sessionInfo.txt"))
  )
}
module_order <- function(x) {
  u <- unique(as.character(x))
  u[order(
    suppressWarnings(as.numeric(gsub("[^0-9]", "", u))),
    u
  )]
}
numeric_matrix <- function(d,
  key = 1,
  cols = NULL) {
  if (is.null(cols))
    cols <- which(vapply(d, is.numeric, logical(1)))
  m <- as.matrix(d[, cols, drop = FALSE])
  rownames(m) <- make.unique(as.character(d[[key]]))
  m
}
source(file.path(REPO, "R/manuscript_style.R"))
source(file.path(REPO, "R/extended_multiome_heatmaps.R"))
