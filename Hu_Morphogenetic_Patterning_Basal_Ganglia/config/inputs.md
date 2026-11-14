# Data inputs and execution

Run commands from this project folder. Set `JIA_DATA_DIR` to a data directory and `JIA_OUTPUT_DIR` to a separate output directory. Nothing writes to the supplied manuscript objects.

## Download the analysis inputs

```bash
python3 scripts/download_public_data.py --include data --output /path/to/jia_release
export JIA_DATA_DIR=/path/to/jia_release/data
```

The [public data index](../docs/public_data.md) links count matrices, annotations, and saved objects individually. The release already includes the count manifests and curated annotations needed by scripts 00–05.

## Saved-object entry point

Place these three supplied files in `$JIA_DATA_DIR/objects/`:

- `MGE_1538_E11_and_E13_with_new_E13.Rds`: 7,785 curated MGE cells (3,064 E11.5; 4,721 E13.5).
- `MZ_reintegrated_on_52_PN_genes_08_15_2025.rds`: 4,047 MZ cells (878 E11.5; 3,169 E13.5), with full-gene expression and the neuronal-panel embedding.
- `VZ_SVZ_reintegrated_on_221_progenitor_genes_08_15_2025.rds`: 3,294 progenitor cells, with full-gene expression and the progenitor-panel embedding.

```bash
export JIA_ANALYSIS_MODE=saved
Rscript scripts/04_developmental_expression_figures.R
Rscript scripts/05_developmental_volcano_plots.R
```

The export produces raw RNA count matrices in 10x Matrix Market format, `count_inputs.tsv`, and `cell_annotations.tsv`. The latter is also distributed in `config/` and must match the count version. It contains cell ID, original 10x barcode, library, age, original global cluster, saved developmental compartment, and the two refined cluster identities. It records selection and annotation, not gene expression or private donor information.

## Count-matrix entry point

Use the matrices exported above for the historical RNA quantification. Alternatively, supply original Cell Ranger matrices and edit `count_inputs.tsv`:

| Column | Meaning |
| --- | --- |
| dataset | `1538TE10H11`, `1538TE10H13`, or `JH4327_01_TE10H13` |
| age | `E11` or `E13` (E11.5/E13.5 in the manuscript) |
| path | Count directory or H5 path relative to `JIA_DATA_DIR` |
| format | `10x` for matrix/features/barcodes; `h5` for a Cell Ranger H5 |
| selection | `curated` for exported retained-cell counts; `unfiltered` to apply QC before retained-cell selection |

Then run scripts 01–05 with `JIA_ANALYSIS_MODE=counts`. H5 input additionally needs `hdf5r`. The export route uses `Matrix` and avoids that dependency.

When starting directly from count files without running the RDS exporter, copy the distributed annotations and count manifest into the data directory, then edit count paths/selection types to match the supplied matrices:

```bash
mkdir -p "$JIA_DATA_DIR"
cp config/cell_annotations.tsv "$JIA_DATA_DIR/cell_annotations.tsv"
cp config/count_inputs.tsv "$JIA_DATA_DIR/count_inputs.tsv"
```

The January 2026 reprocessed matrices are not interchangeable with the original quantification for 1538TE10H11 and 1538TE10H13. Some retained cells are absent from the reprocessed filtered matrices and occur only in the raw matrices; the script stops on missing retained barcodes. A barcode match alone does not establish matching UMI values.

## Software

The original global object was saved with Seurat 4.1.3; the methods document lists R 4.3.3 / Seurat 5.0.3 for the later analysis. The global SCT model is v1, while the full-gene MZ marker model is v2. These flavors and random seeds are explicit in the scripts. `docs/review_environment.txt` records the installed R 4.4.0 / Seurat 5.3.0 environment used for checks, not a historical lockfile. Dependencies: Seurat, SeuratObject, Matrix, sctransform, glmGamPoi, ggplot2, patchwork, ggrepel; the tissue classifier additionally needs scAnnotatR, dplyr, and cowplot. Every analytical stage writes its actual session information.


## Pre-exclusion tissue-classifier entry point

To refit the original tissue-selection method, provide **original, pre-exclusion** matrices for 1538TE10H11 and 1538TE10H13 and a `count_inputs_preselection.tsv` following `count_inputs_preselection.example.tsv`. Run script 00 before the downstream analyses. The supplied curated count exports have already excluded contaminants and cannot be used to train the original exclusion models. Script 00 saves the classifications and their comparison to the accepted manuscript cell selection. Review any differences before replacing the supplied cell annotation table. The code specifies marker quantiles, RNA-count classifier input, train/test sampling, and overlap handling; the original 2023 notebook is included in `historical/`.

The two original pre-exclusion count matrices were recovered from the archived source folders and verified against all retained cells in the manuscript RNA assay: 3,064 E11.5 and 1,374 E13.5 cells from 1538TE10H13, 32,285 shared genes, zero differing gene–cell entries, and zero total-UMI differences. These count versions support the historical R quantification. The recovered folders contain 5,111 and 4,099 Cell Ranger-called cells respectively, before the analysis QC and tissue-selection steps.

An optional `JIA_R_LIBRARY` environment variable prepends a project-specific R library. The workflow stops if glmGamPoi is unavailable for a stage that requires it; it does not silently substitute Seurat's slower native estimator.

To test the classifiers on the exact historical normalized input, set `JIA_TISSUE_CHECKPOINT` to `1st_comparison_1538TE10H11_1538TE10HE13_5969_MGE_POA_Septum_cells_5_25_23.Rds` before running script 00. This bypasses only the preliminary QC/integration already represented by that checkpoint. The count-to-checkpoint preparation remains fully specified in the same script, and the actual checkpoint hash is recorded with the results.
