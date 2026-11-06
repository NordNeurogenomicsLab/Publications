# Usage: Rscript run_all.R [notebook prefix, e.g. 03] Set CORTEX_R_LIBS
# to one or more library directories separated by .Platform$path.sep.
libs <- Sys.getenv("CORTEX_R_LIBS")
if (nzchar(libs)) .libPaths(
  c(
    strsplit(libs, .Platform$path.sep, fixed = TRUE)[[1]],
    .libPaths()
  )
)
source("R/common.R")
args <- commandArgs(TRUE)
files <- sort(list.files(".",
  pattern = "^[0-9].*\\.Rmd$",
  full.names = TRUE))
if (length(args)) files <- files[grepl(
  paste0("^./", args[1]),
  files
)]
if (!length(files)) stop("No matching notebooks")
if (!dir.exists("data") || !dir.exists("objects")) stop(
  "The saved-input data release is required. See DATA_AVAILABILITY.md; the GitHub code checkout alone does not include data/ or objects/."
)
dir.create("reports", showWarnings = FALSE)
for (f in files) {
  cat("\nRUNNING ", f, "\n")
  flush.console()
  out <- file.path("reports", sub("\\.Rmd$", ".md", basename(f)))
  # Knit executes every enabled chunk, retains an inspectable Markdown
  # report, and needs no Pandoc installation. Optional HTML is a
  # separate render step.
  notebook_environment <- new.env(parent = globalenv())
  parameter_definitions <- knitr::knit_params(readLines(f, warn = FALSE))
  notebook_environment$params <- lapply(parameter_definitions, function(x) x$value)
  knitr::knit(
    f,
    output = out,
    envir = notebook_environment,

    quiet = FALSE
  )
  gc()
}
cat("\nALL SELECTED NOTEBOOKS COMPLETED\n")
