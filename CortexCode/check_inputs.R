# Run from CortexCode/: Rscript check_inputs.R [03] [--include-optional]
# Inspects paths and sizes only; no saved objects or scientific models are loaded.
source("R/input_catalog.R")
args <- commandArgs(trailingOnly = TRUE)
prefix <- args[!startsWith(args, "--")]
if (length(prefix) > 1) stop("Specify at most one notebook prefix.")
result <- check_notebook_inputs(
  notebook = if (length(prefix)) prefix else NULL,
  include_optional = "--include-optional" %in% args
)
dir.create("reports", showWarnings = FALSE)
write.csv(result, "reports/input_check.csv", row.names = FALSE)
print(as.data.frame(table(result$Status)), row.names = FALSE)
missing <- result[result$Status != "OK", c("Notebook", "Mode", "Path", "Status")]
if (nrow(missing)) {
  print(missing, row.names = FALSE)
  stop("Inputs are missing or differ in size. See INPUTS.md and reports/input_check.csv.")
}
cat("Declared inputs are present and match recorded sizes. No analysis was run.\n")
