# The catalog is metadata only: these helpers never load an RDS or run analysis.
notebook_inputs <- function(notebook, root = REPO) {
  inputs <- read.csv(file.path(root, "config", "notebook_inputs.csv"),
    stringsAsFactors = FALSE, check.names = FALSE)
  inputs[inputs$Notebook == notebook,
    c("Mode", "Path", "Purpose", "Local_status"), drop = FALSE]
}

check_notebook_inputs <- function(notebook = NULL, include_optional = FALSE,
                                  root = getwd()) {
  inputs <- read.csv(file.path(root, "config", "notebook_input_files.csv"),
    stringsAsFactors = FALSE, check.names = FALSE)
  if (!is.null(notebook)) {
    inputs <- inputs[startsWith(inputs$Notebook, notebook), , drop = FALSE]
    if (!nrow(inputs)) stop("No matching notebook in the input catalog: ", notebook)
  }
  if (!include_optional) inputs <- inputs[inputs$Mode == "saved-results", , drop = FALSE]
  inputs <- inputs[!duplicated(inputs[c("Notebook", "Mode", "Path")]), , drop = FALSE]
  # Relative paths avoid Windows path-length problems in deeply nested checkouts.
  at_root <- identical(normalizePath(root, winslash = "/", mustWork = TRUE),
    normalizePath(getwd(), winslash = "/", mustWork = TRUE))
  paths <- if (at_root) inputs$Path else file.path(root, inputs$Path)
  inputs$Present <- file.exists(paths)
  sizes <- file.info(paths, extra_cols = FALSE)$size
  inputs$Size_matches <- inputs$Present &
    (is.na(inputs$Bytes) | sizes == inputs$Bytes)
  inputs$Status <- ifelse(!inputs$Present, "MISSING",
    ifelse(inputs$Size_matches, "OK", "SIZE MISMATCH"))
  inputs
}

write_session_info <- function(name, root = getwd()) {
  destination <- file.path(root, "environment", "sessions")
  dir.create(destination, recursive = TRUE, showWarnings = FALSE)
  text <- c(
    paste("Reconstruction session; notebook:", name),
    paste("Captured UTC:", format(Sys.time(), tz = "UTC", usetz = TRUE)),
    "This records the current execution environment, not the original analysis environment.",
    "", capture.output(sessionInfo())
  )
  writeLines(text, file.path(destination, paste0(name, "_sessionInfo.txt")))
  invisible(text)
}
