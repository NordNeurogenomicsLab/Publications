# Morphogenetic Patterning During Regional and Cell Type Specification in the Embryonic Basal Ganglia

**Jia Sheng Hu, Karol Cichewicz, Jonathan W. C. Lim, Linda J. Richards, Luis Puelles, Alex S. Nord, and John L. Rubenstein**

Correspondence: John L. Rubenstein (john.rubenstein@ucsf.edu).

[Manuscript](https://doi.org/10.64898/2026.04.14.718571) · [Author affiliations](docs/authors.md)

## Abstract

Cell type specification in the embryonic brain and spinal cord is thought to begin within molecularly defined progenitor domains that do not intermix. Our data provide an alternative model that is spatially and temporally dynamic within a basal ganglia anlage, the medial ganglionic eminence (MGE). MGE progenitor cells are progressively displaced ventrally and caudally from a rostral growth zone (the MGE/LGE sulcus). Progenitors that leave the MGE/LGE sulcus early occupy caudoventral MGE regions, while ones that leave later reside in rostrodorsal MGE regions. As they change position, their transcriptional states and cell type output change. Transcriptional analyses showed an upregulation of the Nfi TFs during the period of progenitor movement. Nfia and Nfib double mutants alter the repertoire of cortical interneuron subtypes. Overall, we present a mechanism that synchronizes regional patterning with tissue growth and links spatial and temporal specification in producing diverse neuronal subtypes.

## Single-cell RNA-seq analysis

The analysis compares E11.5 and E13.5 medial ganglionic eminence (MGE) cells labeled by the 1538 lineage following E10.5 tamoxifen administration. It includes count-matrix import, quality control, MGE cell selection, SCT normalization and integration, dimensionality reduction and clustering, developmental-stage differential expression, and gene-expression and cell-composition figures.

| Script | Scientific content |
| --- | --- |
| `scripts/00_tissue_selection_from_counts.R` | Quality control and POA/septum cell exclusion in 1538TE10H11 (E11.5) and 1538TE10H13 (E13.5) |
| `scripts/01_count_matrix_integration.R` | SCT integration of 1538TE10H11 (E11.5), 1538TE10H13 (E13.5), and JH4327_01_TE10H13 (E13.5) |
| `scripts/02_neuronal_progenitor_subclustering.R` | MZ analysis with the 52-gene neuronal panel and VZ/SVZ analysis with the progenitor panel; preserve full-gene expression for marker analysis |
| `scripts/03_developmental_differential_expression.R` | SCT-based E11.5 versus E13.5 tests within developmental compartments and neuronal/progenitor clusters |
| `scripts/04_developmental_expression_figures.R` | Global and compartment UMAPs, expression maps, cluster proportions, and marker heatmaps |
| `scripts/05_developmental_volcano_plots.R` | Volcano plots of E11.5 versus E13.5 differential expression |
| `scripts/tissue_identity_classification.R` | Document the marker-guided POA/septum classification procedure and its train/test split |

```bash
export JIA_DATA_DIR=/path/to/analysis_inputs
export JIA_OUTPUT_DIR=/path/to/analysis_outputs
Rscript scripts/01_count_matrix_integration.R
Rscript scripts/02_neuronal_progenitor_subclustering.R
Rscript scripts/03_developmental_differential_expression.R
Rscript scripts/04_developmental_expression_figures.R
Rscript scripts/05_developmental_volcano_plots.R
```

## R session information

```r
sessionInfo()
```

[Session details](docs/sessionInfo.txt).

## Data availability

The [public data release](docs/public_data.md) contains **46 files (4.41 GB)**: count matrices, cell annotations, manuscript Seurat objects, analysis checkpoints, result tables, and Tables S1–S4. Each file has a direct download URL and SHA-256 checksum in the [file manifest](config/PUBLIC_DATA_FILES.tsv). No sign-in is required.

| Data | Download |
| --- | --- |
| Cell annotations: library, age, compartment, and cluster assignments | [Annotations](https://storage.googleapis.com/hu-morphogenetic-patterning-306023/releases/2026-11-13/data/cell_annotations.tsv) |
| 1538TE10H11 (E11.5), accepted MGE cells | [UMI count matrix](https://storage.googleapis.com/hu-morphogenetic-patterning-306023/releases/2026-11-13/data/counts/1538TE10H11/matrix.mtx) |
| 1538TE10H13 (E13.5), accepted MGE cells | [UMI count matrix](https://storage.googleapis.com/hu-morphogenetic-patterning-306023/releases/2026-11-13/data/counts/1538TE10H13/matrix.mtx) |
| JH4327_01_TE10H13 (E13.5), accepted MGE cells | [UMI count matrix](https://storage.googleapis.com/hu-morphogenetic-patterning-306023/releases/2026-11-13/data/counts/JH4327_01_TE10H13/matrix.mtx) |
| Global MGE Seurat object: 7,785 cells | [RDS](https://storage.googleapis.com/hu-morphogenetic-patterning-306023/releases/2026-11-13/data/objects/MGE_1538_E11_and_E13_with_new_E13.Rds) |
| MZ Seurat object: 4,047 cells | [RDS](https://storage.googleapis.com/hu-morphogenetic-patterning-306023/releases/2026-11-13/data/objects/MZ_reintegrated_on_52_PN_genes_08_15_2025.rds) |
| VZ/SVZ Seurat object: 3,294 cells | [RDS](https://storage.googleapis.com/hu-morphogenetic-patterning-306023/releases/2026-11-13/data/objects/VZ_SVZ_reintegrated_on_221_progenitor_genes_08_15_2025.rds) |

Count-matrix feature and barcode files, pre-exclusion matrices, differential-expression tables, supplementary tables, and calculated checkpoints are linked in the [complete data index](docs/public_data.md). [Download instructions](docs/public_data.md#download) preserve the directory layout used by the scripts.

Raw sequencing reads will be linked through GEO/SRA once an accession is assigned.
