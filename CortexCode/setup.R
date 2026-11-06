# Read-only dependency check; no package installation or upgrades.
cran <- c(
  "ggplot2", "dplyr", "tidyr", "jsonlite", "pheatmap", "knitr", "rmarkdown",
  "patchwork", "viridisLite", "RColorBrewer", "scales", "ggrepel", "circlize"
)
bioc <- c("ComplexHeatmap")
required <- c(cran, bioc)
ok <- vapply(
  required, requireNamespace, logical(1),
  quietly = TRUE
)
print(
  data.frame(
    Package = required,
    Installed = ok,
    Version = vapply(
      required, function(p) if (requireNamespace(p, quietly = TRUE)) as.character(packageVersion(p)) else NA_character_,
      character(1)
    )
  )
)
if (!all(ok)) stop(
  "Install the missing packages in your chosen R environment. ComplexHeatmap is a Bioconductor package. No packages were modified."
)
cat(
  "Default plotting dependencies are available. Optional historical regeneration additionally uses the packages listed in the original scripts.\n"
)
