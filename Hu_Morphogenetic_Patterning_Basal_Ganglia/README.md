# Hu: Morphogenetic Patterning in the Embryonic Basal Ganglia

Code for *Morphogenetic Patterning During Regional and Cell Type Specification in the Embryonic Basal Ganglia*.

## Status
Project folder prepared for Publications; scientific review remains in progress. Raw and processed GEO copies are verified; metadata remains a draft. No accession or license is invented. Confirm author ownership, contributors and license before public release.

## Recreate saved panels
Set JIA_WORK_DIR to the directory containing inputs/, figures/, and audit/. On sclab this is /mnt/disks/data1/Jia_GEO_submission/codex_review. Required saved objects and DE tables are listed in config/inputs.md.

```
Rscript scripts/recreate_saved_panels.R
Rscript scripts/recreate_VZ_documented_clusters.R
Rscript scripts/recreate_volcano.R
```

The documented VZ cluster script supersedes the first script's VZ plots based on incomplete Dev_classes metadata. Do not run alignment, normalization, integration or clustering to recreate these saved plots.

The review environment is R 4.4.0 with Seurat, SeuratObject, ggplot2, patchwork, ggrepel and dependencies. Exact review environment versions are in docs/review_environment.txt. These are review versions, not claims about the historical manuscript environment.

## Provenance and limitations
Historical Rmd files preserve original exploratory code and absolute paths. They are evidence, not a turnkey pipeline; do not execute them blindly. Exact historical Cell Ranger settings for two older libraries remain unresolved. The eight named Figure 4D groups are not stored as a recoverable mapping. Candidate marker-based assignments require review. Manual image annotations and plot styling differ.

Large data, credentials and original manuscript images are excluded. Link the final GEO accession after approval and deposition. The working GEO audit is maintained separately from this publication folder.

## Manuscript
Hu et al., bioRxiv v1 (2026), https://doi.org/10.64898/2026.04.14.718571

Historical VZ/MZ DE tables are included under tables/. No new DE tests were run. Provisional marker-based annotations are not a recovered manuscript key.
