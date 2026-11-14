# Computational verification

The structured workflow was executed with R 4.4.0, Seurat 5.3.0, glmGamPoi 1.16.0, and scAnnotatR 1.10.0. `renv.lock` records the full dependency set.

## Count inputs and analysis objects

The two recovered original Cell Ranger count matrices agree exactly with all retained cells' RNA counts in the manuscript object: zero differing gene–cell entries and zero total-UMI differences. The comparison covers 3,064 E11.5 and 1,374 older-library E13.5 cells and 32,285 genes. See `count_matrix_verification.tsv` and `config/count_matrix_checksums.tsv`.

Count import, global integration, both cell-state subclustering analyses, developmental differential expression, marker testing, and figure generation completed successfully. Serialized objects passed Seurat validity checks, cell-ID/embedding alignment checks, and the intended model-type checks:

| Object | Cells | E11.5 | E13.5 | SCT fitting routine |
| --- | ---: | ---: | ---: | --- |
| Global MGE | 7,785 | 3,064 | 4,721 | v1 poisson, two reference models |
| Neuronal MZ | 4,047 | 878 | 3,169 | v2 glmGamPoi_offset |
| Progenitor checkpoint | 3,294 | 2,185 | 1,109 | inherited v1 poisson models |
| Progenitor cluster-marker model | 3,294 | 2,185 | 1,109 | separate v2 glmGamPoi_offset |

The global and MZ means of detected genes are 4,613.958 and 4,628.570, matching the manuscript's rounded 4,614 and 4,629. The count workflow generated 26 computational plots, including marker heatmaps. The supplied-object entry point generated 24 expression, UMAP, proportion, and volcano plots. Representative final plots were visually inspected for labels, axes, faceting, and rendering.

## Tissue-selection check

The classifier module was executed on the original 5,969-cell normalized checkpoint. The current refit produced held-out AUCs of 0.988952 (POA) and 0.986231 (septum). It retained 4,398 MGE cells, differing at 40 cell labels from the manuscript's accepted 4,438-cell older-library selection. The recovered July 18 workspace contains an earlier classifier state and is not treated as the final selection model.

The main count workflow uses the explicit accepted cell annotation table, preserving the manuscript selection. A classifier refit is compared to that table and does not overwrite it. The final fitted historical models or any intervening annotation decisions are needed to explain the 40-cell difference. No thresholds were tuned to force agreement.

The preliminary count-to-normalized-checkpoint code is supplied in script 00. Classifier execution was validated using the original normalized checkpoint; its optional full preliminary normalization refit was not used as the verification entry point.

## Anatomical categories and supplementary identifiers

Figure 4D's final eight-group cell annotation still requires an author-confirmed key. Table S4 is a marker reference and Figure S7 has regional annotations, but neither supplies the complete deterministic per-cell aggregation. Numerical-cluster proportions are calculated directly and are not relabeled with an inferred key.

Ten gene-name cells in the original supplementary spreadsheets were converted to dates. The analysis CSVs preserve the original symbols; the affected cells are recorded in `supplement_gene_symbol_dates.tsv`.
