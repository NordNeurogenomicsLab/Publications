## Manuscript

**Profiling dynamic epigenetic programs in the developing mouse cortex reveals the role of ASD-associated Kdm6b in deep-layer neuron differentiation.**

Athéna R. Ypsilanti<sup>†,*</sup>, Eirene Markenscoff-Papadimitriou<sup>†,*</sup>, Karol Cichewicz<sup>†</sup>, Sean Whalen, Nicolas Seban, Ethan Fenton, Kenneth Lim, Shan Dong, Jada Baldoza, Yurong Wang, Matt Barraza, Fadya Binyameen, Diane Dickel, Len A. Pennacchio, Axel Visel, Stephan J. Sanders, Katherine S. Pollard, Alex S. Nord<sup>*</sup>, John Rubenstein<sup>*</sup>

† Contributed equally: Athéna R. Ypsilanti, Eirene Markenscoff-Papadimitriou and Karol Cichewicz.

\* Corresponding authors: Athéna R. Ypsilanti, Eirene Markenscoff-Papadimitriou, Alex S. Nord and John Rubenstein.

# CortexCode analysis and figures

These notebooks document the chromatin-module, multiome, integration, differential-chromatin and enrichment analyses and export the manuscript's available quantitative panels. Key scientific algorithms are included directly in the figure notebooks. Their analysis parameters default to `false`; a normal run loads the saved, verified results.

## Run from saved results

The GitHub checkout contains code and provenance records. Restore the separately supplied `data/` and `objects/` directories into this directory before running; a permanent public saved-input release has not yet been linked.

Run from the repository directory with the required R packages installed:

```sh
Rscript run_all.R
```

Alternatively, use `Rscript run_all.R 03` to knit Figure 3 alone. `run_all.R` produces Markdown reports without Pandoc. With Pandoc installed, render HTML with `rmarkdown::render("03_Figure_3_multiome.Rmd")`. `setup.R` checks dependencies; it does not install packages automatically. `CORTEX_R_LIBS` can specify additional library directories separated by the platform path separator.

## Analysis notebooks

| Notebook | Scientific code included |
|---|---|
| `00_inputs_and_provenance.Rmd` | Input schemas, file identities and plotting-table checkpoints |
| `01_WGCNA_checkpoints.Rmd` | Count preparation, histone interval merging, CPM, quality masks, soft-threshold diagnostics, complete WGCNA fit settings and hierarchical grey refinement |
| `02_Figures_1_2_modules.Rmd` | GREAT submission and GO filtering; HOMER BED/background generation and command; motif selection, sequence collapse and curated matrices; Figures 1â€“2 |
| `03_Figure_3_multiome.Rmd` | Module pooling, native Seurat scoring and controls, link selection, DotPlot statistics and score scaling; Figure 3 |
| `04_Figure_4_integration.Rmd` | Expanded interval overlap construction and singletons, stage medians, ordinal encoding, trajectory assignment, complete feature/UMAP/graph-clustering routines, selected settings and manual program definitions |
| `05_Figure_5_RNA_histology.Rmd` | Explicit mouse topGO implementation and input selection; pseudobulk VST/PCA; saved bulk RNA results and mouse-level histology summaries |
| `06_Figure_6_epigenomics.Rmd` | Complete original DiffBind workflow, peak expansion, counting, normalization, DESeq2 contrasts, module-enrichment statistics and expression-normalization formula |
| `07_Figure_6_variants.Rmd` | Fu gene sets, orthology and gene nulls; complete native balanced case/control variant-downsampling class and selected invocation |
| `08_Figure_7_enhancer_histology.Rmd` | Mouse-level enhancer-reporter measurements and histology panels |
| `10_manuscript_graphics.Rmd` | Shared manuscript styling and saved vector exports |
| `11_additional_manuscript_panels.Rmd` | Additional manuscript groupings, UMAP exports and count panels |
| `12_Extended_multiome_heatmaps.Rmd` | Saved-layer extraction, AC/ME3 linked-gene scores, explicit marker selection, and supplemental heatmap ordering/scaling |

The smaller notebooks in `pipeline/` provide supplementary input/method descriptions. They are not required to find the principal algorithms: these appear in the figure notebooks above. See [ANALYSIS_METHODS.md](ANALYSIS_METHODS.md) for parameters, execution controls and source identities.

## Deliberate analysis reruns

Enable a particular parameter explicitly, for example:

```r
rmarkdown::render(
  "01_WGCNA_checkpoints.Rmd",
  params = list(run_wgcna = TRUE)
)
```

This example documents the invocation; routine rendering leaves `run_wgcna = FALSE`. Other notebooks expose separate flags for GREAT, HOMER, multiome scoring, integration, topGO, DiffBind and resampling. Rerun results are written separately from archived inputs. Saved plotting inputs are not silently replaced by a new fit. Live annotation services and updated package/database versions can change results and require comparison with the archived tables.

## Data, objects and provenance

- `data/source_tables/`: manuscript numerical plotting inputs.
- `data/upstream/`: raw WGCNA count tables, final interval catalogs, original RNA statistics, full GREAT tables, HOMER known-motif results and target/background BEDs, multiome inputs and portable DiffBind metadata.
- `objects/`: fitted WGCNA networks, exact module-score/control checkpoints, compact stored-assay cache, original integration bundles and original resampling results.
- `R/`: shared plotting and data helpers. Key analysis classes/functions are also shown directly in the relevant notebooks.
- `provenance/`: source-script identities and SHA-256 hashes, input manifests, panel coverage, palettes and recorded HOMER commands. Historical HOMER mount paths are replaced by `${HOMER_ROOT}` in command records.
- `panels/` and `reports/`: generated exports and execution reports.

Large inputs, R objects and generated reports are saved in the local delivery but excluded from ordinary Git commits by `.gitignore`; they require a separate data deposit. This directory is not itself a completed public data deposit. [DATA_AVAILABILITY.md](DATA_AVAILABILITY.md) lists external inputs and unresolved historical details.

## Reproduction limits

The original Figure 5d GO output/configuration and original bulk-RNA model-fitting script have not been recovered. The supplied topGO code is an explicit mouse adaptation of the older Neurona implementation; it is not falsely labeled as the recovered Figure 5d run. Two recovered-method differences also need reconciliation: hierarchical grey splitting versus the manuscript's k-means description, and the selected HOMER module-complement background versus the GC-matched global background described in the manuscript. HOMER native âˆ’ln(P) values must also be distinguished from âˆ’log10(P).

Figure 6f uses the available original geneâ€“peak-overlap analysis (484 down / 250 up), whose counting unit differs from the unresolved final-artwork totals (362 / 281). Coverage and display approximations are recorded in `provenance/panel_coverage.csv`. These qualifications do not alter the supplied source statistics.

Figure assembly utilities, comparison slides and development archives are excluded from the public code release. Numerical notebooks export the available computational panels; they do not generate Illustrator schematics, browser screenshots or microscopy images.

## Related original analysis repositories

- [Kdm6b cKO differential expression](https://github.com/NordNeurogenomicsLab/CortexCode_Kdm6b_DE)
- [ncDNVs and gene-set enrichment](https://github.com/NordNeurogenomicsLab/CortexCode_enrichment)

The earlier [WGCNA notebook](WGCNA_network_generation_w_1kb_peak_merging.Rmd) is retained for historical reference; start with the numbered notebooks above for the recreated workflow.

Original rendered reports: [Kdm6b differential expression](https://nordneurogenomicslab.github.io/CortexCode_Kdm6b_DE/) and [ncDNV enrichment](https://nordneurogenomicslab.github.io/CortexCode_enrichment/).
