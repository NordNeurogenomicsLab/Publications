# Public single-cell analysis data

Release: **2026-11-13**, **46 scientific files**, **4,405,155,826 bytes** (4.41 GB; 4.10 GiB). Downloads are available without authentication.

[File manifest](../config/PUBLIC_DATA_FILES.tsv) · [Hosted manifest](https://storage.googleapis.com/hu-morphogenetic-patterning-306023/releases/2026-11-13/FILES.tsv) · [SHA-256 checksums](https://storage.googleapis.com/hu-morphogenetic-patterning-306023/releases/2026-11-13/SHA256SUMS)

## Download

From the manuscript code directory:

```bash
python3 scripts/download_public_data.py --output /path/to/jia_release
export JIA_DATA_DIR=/path/to/jia_release/data
export JIA_OUTPUT_DIR=/path/to/analysis_outputs
```

The downloader checks file sizes and SHA-256 checksums and skips files that already match. To download only analysis inputs, add `--include data`. For result tables and checkpoints, use `--include results`; for the original Tables S1–S4, use `--include supplementary_tables`. Multiple `--include` options may be combined. Individual files can also be downloaded from the links below. Count matrices need all three files in their library directory. See [input configuration and execution](../config/inputs.md) for the analysis commands.

## File contents

`data/objects/` contains the three supplied manuscript Seurat objects, retaining their stored cell assignments, expression assays, and embeddings. `data/counts/` contains RNA UMI counts for the accepted MGE cells; `data/original_counts/` contains the pre-exclusion Cell Ranger-called cells for 1538TE10H11 and 1538TE10H13. The annotation and count manifests use paths relative to `JIA_DATA_DIR`.

`results/count_outputs/` contains newly calculated results and checkpoints from the count-based scripts. `results/saved_outputs/` contains figure data calculated from the supplied manuscript objects. These folders are identified separately because fitted models and numerical cluster labels may differ between runs. The [scientific methods](scientific_methods.md) describe the analysis and annotation requirements.

`TableS1.xlsx`–`TableS4.xlsx` are unchanged source workbooks. The [gene-symbol correction record](supplement_gene_symbol_dates.tsv) documents symbols that Excel converted to dates; analysis CSVs retain gene symbols.

## Analysis inputs

| File | MB | Contents |
| --- | ---: | --- |
| [`data/cell_annotations.tsv`](https://storage.googleapis.com/hu-morphogenetic-patterning-306023/releases/2026-11-13/data/cell_annotations.tsv) | 0.53 | Manuscript analysis input |
| [`data/count_inputs.tsv`](https://storage.googleapis.com/hu-morphogenetic-patterning-306023/releases/2026-11-13/data/count_inputs.tsv) | 0.00 | Manuscript analysis input |
| [`data/counts/1538TE10H11/barcodes.tsv`](https://storage.googleapis.com/hu-morphogenetic-patterning-306023/releases/2026-11-13/data/counts/1538TE10H11/barcodes.tsv) | 0.06 | Manuscript analysis input |
| [`data/counts/1538TE10H11/features.tsv`](https://storage.googleapis.com/hu-morphogenetic-patterning-306023/releases/2026-11-13/data/counts/1538TE10H11/features.tsv) | 1.01 | Manuscript analysis input |
| [`data/counts/1538TE10H11/matrix.mtx`](https://storage.googleapis.com/hu-morphogenetic-patterning-306023/releases/2026-11-13/data/counts/1538TE10H11/matrix.mtx) | 162.75 | Manuscript analysis input |
| [`data/counts/1538TE10H13/barcodes.tsv`](https://storage.googleapis.com/hu-morphogenetic-patterning-306023/releases/2026-11-13/data/counts/1538TE10H13/barcodes.tsv) | 0.03 | Manuscript analysis input |
| [`data/counts/1538TE10H13/features.tsv`](https://storage.googleapis.com/hu-morphogenetic-patterning-306023/releases/2026-11-13/data/counts/1538TE10H13/features.tsv) | 1.01 | Manuscript analysis input |
| [`data/counts/1538TE10H13/matrix.mtx`](https://storage.googleapis.com/hu-morphogenetic-patterning-306023/releases/2026-11-13/data/counts/1538TE10H13/matrix.mtx) | 58.81 | Manuscript analysis input |
| [`data/counts/JH4327_01_TE10H13/barcodes.tsv`](https://storage.googleapis.com/hu-morphogenetic-patterning-306023/releases/2026-11-13/data/counts/JH4327_01_TE10H13/barcodes.tsv) | 0.06 | Manuscript analysis input |
| [`data/counts/JH4327_01_TE10H13/features.tsv`](https://storage.googleapis.com/hu-morphogenetic-patterning-306023/releases/2026-11-13/data/counts/JH4327_01_TE10H13/features.tsv) | 1.01 | Manuscript analysis input |
| [`data/counts/JH4327_01_TE10H13/matrix.mtx`](https://storage.googleapis.com/hu-morphogenetic-patterning-306023/releases/2026-11-13/data/counts/JH4327_01_TE10H13/matrix.mtx) | 220.50 | Manuscript analysis input |
| [`data/objects/MGE_1538_E11_and_E13_with_new_E13.Rds`](https://storage.googleapis.com/hu-morphogenetic-patterning-306023/releases/2026-11-13/data/objects/MGE_1538_E11_and_E13_with_new_E13.Rds) | 880.36 | Manuscript analysis input |
| [`data/objects/MZ_reintegrated_on_52_PN_genes_08_15_2025.rds`](https://storage.googleapis.com/hu-morphogenetic-patterning-306023/releases/2026-11-13/data/objects/MZ_reintegrated_on_52_PN_genes_08_15_2025.rds) | 452.55 | Manuscript analysis input |
| [`data/objects/VZ_SVZ_reintegrated_on_221_progenitor_genes_08_15_2025.rds`](https://storage.googleapis.com/hu-morphogenetic-patterning-306023/releases/2026-11-13/data/objects/VZ_SVZ_reintegrated_on_221_progenitor_genes_08_15_2025.rds) | 378.76 | Manuscript analysis input |
| [`data/original_counts/1538TE10H11/barcodes.tsv.gz`](https://storage.googleapis.com/hu-morphogenetic-patterning-306023/releases/2026-11-13/data/original_counts/1538TE10H11/barcodes.tsv.gz) | 0.03 | Pre-exclusion Cell Ranger counts |
| [`data/original_counts/1538TE10H11/features.tsv.gz`](https://storage.googleapis.com/hu-morphogenetic-patterning-306023/releases/2026-11-13/data/original_counts/1538TE10H11/features.tsv.gz) | 0.29 | Pre-exclusion Cell Ranger counts |
| [`data/original_counts/1538TE10H11/matrix.mtx.gz`](https://storage.googleapis.com/hu-morphogenetic-patterning-306023/releases/2026-11-13/data/original_counts/1538TE10H11/matrix.mtx.gz) | 81.26 | Pre-exclusion Cell Ranger counts |
| [`data/original_counts/1538TE10H13/barcodes.tsv.gz`](https://storage.googleapis.com/hu-morphogenetic-patterning-306023/releases/2026-11-13/data/original_counts/1538TE10H13/barcodes.tsv.gz) | 0.02 | Pre-exclusion Cell Ranger counts |
| [`data/original_counts/1538TE10H13/features.tsv.gz`](https://storage.googleapis.com/hu-morphogenetic-patterning-306023/releases/2026-11-13/data/original_counts/1538TE10H13/features.tsv.gz) | 0.29 | Pre-exclusion Cell Ranger counts |
| [`data/original_counts/1538TE10H13/matrix.mtx.gz`](https://storage.googleapis.com/hu-morphogenetic-patterning-306023/releases/2026-11-13/data/original_counts/1538TE10H13/matrix.mtx.gz) | 32.93 | Pre-exclusion Cell Ranger counts |
| [`data/count_inputs_preselection.tsv`](https://storage.googleapis.com/hu-morphogenetic-patterning-306023/releases/2026-11-13/data/count_inputs_preselection.tsv) | 0.00 | Pre-exclusion count manifest |

## Count-based analysis results

| File | MB | Contents |
| --- | ---: | --- |
| [`results/count_outputs/MZ_cell_metadata.csv`](https://storage.googleapis.com/hu-morphogenetic-patterning-306023/releases/2026-11-13/results/count_outputs/MZ_cell_metadata.csv) | 0.49 | Count-based analysis result |
| [`results/count_outputs/MZ_cluster_age_proportions.csv`](https://storage.googleapis.com/hu-morphogenetic-patterning-306023/releases/2026-11-13/results/count_outputs/MZ_cluster_age_proportions.csv) | 0.00 | Count-based analysis result |
| [`results/count_outputs/MZ_cluster_markers.csv`](https://storage.googleapis.com/hu-morphogenetic-patterning-306023/releases/2026-11-13/results/count_outputs/MZ_cluster_markers.csv) | 1.30 | Count-based analysis result |
| [`results/count_outputs/MZ_markers.csv`](https://storage.googleapis.com/hu-morphogenetic-patterning-306023/releases/2026-11-13/results/count_outputs/MZ_markers.csv) | 0.41 | Count-based analysis result |
| [`results/count_outputs/MZ_volcano_data.csv`](https://storage.googleapis.com/hu-morphogenetic-patterning-306023/releases/2026-11-13/results/count_outputs/MZ_volcano_data.csv) | 0.59 | Count-based analysis result |
| [`results/count_outputs/MZ_within_cluster_age_markers.csv`](https://storage.googleapis.com/hu-morphogenetic-patterning-306023/releases/2026-11-13/results/count_outputs/MZ_within_cluster_age_markers.csv) | 5.39 | Count-based analysis result |
| [`results/count_outputs/SVZ_markers.csv`](https://storage.googleapis.com/hu-morphogenetic-patterning-306023/releases/2026-11-13/results/count_outputs/SVZ_markers.csv) | 0.38 | Count-based analysis result |
| [`results/count_outputs/VZ_SVZ_cell_metadata.csv`](https://storage.googleapis.com/hu-morphogenetic-patterning-306023/releases/2026-11-13/results/count_outputs/VZ_SVZ_cell_metadata.csv) | 0.65 | Count-based analysis result |
| [`results/count_outputs/VZ_SVZ_cluster_markers.csv`](https://storage.googleapis.com/hu-morphogenetic-patterning-306023/releases/2026-11-13/results/count_outputs/VZ_SVZ_cluster_markers.csv) | 1.01 | Count-based analysis result |
| [`results/count_outputs/VZ_SVZ_within_cluster_age_markers.csv`](https://storage.googleapis.com/hu-morphogenetic-patterning-306023/releases/2026-11-13/results/count_outputs/VZ_SVZ_within_cluster_age_markers.csv) | 3.81 | Count-based analysis result |
| [`results/count_outputs/VZ_markers.csv`](https://storage.googleapis.com/hu-morphogenetic-patterning-306023/releases/2026-11-13/results/count_outputs/VZ_markers.csv) | 0.43 | Count-based analysis result |
| [`results/count_outputs/VZ_volcano_data.csv`](https://storage.googleapis.com/hu-morphogenetic-patterning-306023/releases/2026-11-13/results/count_outputs/VZ_volcano_data.csv) | 0.61 | Count-based analysis result |
| [`results/count_outputs/count_import_qc.csv`](https://storage.googleapis.com/hu-morphogenetic-patterning-306023/releases/2026-11-13/results/count_outputs/count_import_qc.csv) | 0.00 | Count-based analysis result |
| [`results/count_outputs/global_cluster_correspondence.csv`](https://storage.googleapis.com/hu-morphogenetic-patterning-306023/releases/2026-11-13/results/count_outputs/global_cluster_correspondence.csv) | 0.00 | Count-based analysis result |
| [`results/count_outputs/objects/MGE_1538_E11_and_E13_with_new_E13.Rds`](https://storage.googleapis.com/hu-morphogenetic-patterning-306023/releases/2026-11-13/results/count_outputs/objects/MGE_1538_E11_and_E13_with_new_E13.Rds) | 880.27 | Count-based analysis checkpoint |
| [`results/count_outputs/objects/MZ_reintegrated_on_52_PN_genes_08_15_2025.rds`](https://storage.googleapis.com/hu-morphogenetic-patterning-306023/releases/2026-11-13/results/count_outputs/objects/MZ_reintegrated_on_52_PN_genes_08_15_2025.rds) | 455.35 | Count-based analysis checkpoint |
| [`results/count_outputs/objects/VZ_SVZ_cluster_marker_model.rds`](https://storage.googleapis.com/hu-morphogenetic-patterning-306023/releases/2026-11-13/results/count_outputs/objects/VZ_SVZ_cluster_marker_model.rds) | 391.49 | Count-based analysis checkpoint |
| [`results/count_outputs/objects/VZ_SVZ_reintegrated_on_221_progenitor_genes_08_15_2025.rds`](https://storage.googleapis.com/hu-morphogenetic-patterning-306023/releases/2026-11-13/results/count_outputs/objects/VZ_SVZ_reintegrated_on_221_progenitor_genes_08_15_2025.rds) | 389.30 | Count-based analysis checkpoint |

## Manuscript-object figure data

| File | MB | Contents |
| --- | ---: | --- |
| [`results/saved_outputs/MZ_cluster_age_proportions.csv`](https://storage.googleapis.com/hu-morphogenetic-patterning-306023/releases/2026-11-13/results/saved_outputs/MZ_cluster_age_proportions.csv) | 0.00 | Manuscript-object figure data |
| [`results/saved_outputs/MZ_volcano_data.csv`](https://storage.googleapis.com/hu-morphogenetic-patterning-306023/releases/2026-11-13/results/saved_outputs/MZ_volcano_data.csv) | 0.40 | Manuscript-object figure data |
| [`results/saved_outputs/VZ_volcano_data.csv`](https://storage.googleapis.com/hu-morphogenetic-patterning-306023/releases/2026-11-13/results/saved_outputs/VZ_volcano_data.csv) | 0.37 | Manuscript-object figure data |

## Original supplementary tables

| File | MB | Contents |
| --- | ---: | --- |
| [`supplementary_tables/TableS1.xlsx`](https://storage.googleapis.com/hu-morphogenetic-patterning-306023/releases/2026-11-13/supplementary_tables/TableS1.xlsx) | 0.30 | Original manuscript supplementary table |
| [`supplementary_tables/TableS2.xlsx`](https://storage.googleapis.com/hu-morphogenetic-patterning-306023/releases/2026-11-13/supplementary_tables/TableS2.xlsx) | 0.32 | Original manuscript supplementary table |
| [`supplementary_tables/TableS3.xlsx`](https://storage.googleapis.com/hu-morphogenetic-patterning-306023/releases/2026-11-13/supplementary_tables/TableS3.xlsx) | 0.01 | Original manuscript supplementary table |
| [`supplementary_tables/TableS4.xlsx`](https://storage.googleapis.com/hu-morphogenetic-patterning-306023/releases/2026-11-13/supplementary_tables/TableS4.xlsx) | 0.01 | Original manuscript supplementary table |

## Release verification

All 49 hosted files (46 scientific files plus the manifest, checksum file, and release description) were verified against the staged copies by file size and CRC32C. MD5 was also checked where supplied by Cloud Storage. Every download URL was tested anonymously, and the downloaded annotation table matched its SHA-256 checksum. These checks were completed on 2026-11-13. The release uses a dated path; later revisions should use a new release directory.
