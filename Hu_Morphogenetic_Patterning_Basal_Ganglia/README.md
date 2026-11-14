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
| `scripts/00_tissue_selection_from_counts.R` | Pre-exclusion QC, preliminary integration, and POA/septum classifier fitting for the older libraries |
| `scripts/01_count_matrix_integration.R` | Import count matrices, apply documented QC and retained-cell annotations, combine the earlier libraries, and integrate the additional E13.5 library |
| `scripts/02_neuronal_progenitor_subclustering.R` | MZ analysis with the 52-gene neuronal panel and VZ/SVZ analysis with the progenitor panel; preserve full-gene expression for marker analysis |
| `scripts/03_developmental_differential_expression.R` | SCT-based E11.5 versus E13.5 tests within developmental compartments and neuronal/progenitor clusters |
| `scripts/04_developmental_expression_figures.R` | Global and compartment UMAPs, expression maps, cluster proportions, and marker heatmaps |
| `scripts/05_developmental_volcano_plots.R` | Developmental differential-expression volcano plots with an explicit fold-change convention |
| `scripts/export_analysis_inputs.R` | Export the original RNA count matrices and curated cell annotations from the supplied manuscript objects |
| `scripts/tissue_identity_classification.R` | Document the marker-guided POA/septum classification procedure and its train/test split |

Methods, parameters, gene lists, and reasons for the analysis choices are given in the scripts and [scientific methods](docs/scientific_methods.md). [Data inputs and execution](config/inputs.md) describe both a count-matrix workflow and entry points using saved Seurat objects. Script output directories are separate from the original data.

```bash
export JIA_DATA_DIR=/path/to/analysis_inputs
export JIA_OUTPUT_DIR=/path/to/analysis_outputs
Rscript scripts/01_count_matrix_integration.R
Rscript scripts/02_neuronal_progenitor_subclustering.R
Rscript scripts/03_developmental_differential_expression.R
Rscript scripts/04_developmental_expression_figures.R
Rscript scripts/05_developmental_volcano_plots.R
```

For figures from the supplied manuscript objects, set `JIA_ANALYSIS_MODE=saved` and run scripts 04–05. Script 03 can also analyze those objects. The saved objects retain the manuscript cell identities and embeddings; the count workflow calculates new embeddings and clusters and retains the manuscript annotations in separate columns for comparison. A fresh run is evaluated against the reference cell assignments rather than assuming numeric cluster IDs remain identical.

## Computational environment

`renv.lock` pins the R 4.4.0 environment used to execute the structured workflow, including Seurat, glmGamPoi, and scAnnotatR. It records the executable environment for these scripts; historical object versions are documented separately.

```bash
Rscript -e 'if (!requireNamespace("renv", quietly=TRUE)) install.packages("renv", repos="https://cloud.r-project.org"); renv::restore(lockfile="renv.lock", library="renv/library", prompt=FALSE)'
export JIA_R_LIBRARY="$PWD/renv/library"
```

[Computational verification](docs/computational_verification.md) records the execution checks, input identity tests, and annotation requirements.

## Data availability

Count matrices and Seurat objects are external data dependencies listed in [config/inputs.md](config/inputs.md). Large data files are not stored in Git. A GEO accession has not yet been assigned in this repository. The original analysis notebooks are retained in `historical/`; the structured scripts describe the methods used in the single-cell analyses.
