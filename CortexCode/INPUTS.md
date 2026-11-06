# Input files and execution map

All paths below are relative to `CortexCode/`, the folder containing `run_all.R`.
The repository includes selected compact plotting tables, histology grids, metadata and original WGCNA count tables. [TABLES_ON_GITHUB.csv](data/TABLES_ON_GITHUB.csv) lists every deposited table with size and SHA-256. Larger matrices and binary objects remain outside Git and require a separate input archive; a fresh clone cannot yet reproduce every panel. The file inventory marks each file's availability.

## Locate the biological inputs

| Input | Exact location | Meaning |
|---|---|---|
| Developmental ATAC intervals | `data/upstream/final_tables/Supp_Table_1_ATAC_master.csv` | Final annotated peak/module catalog, not raw sequencing reads |
| Developmental H3K27ac intervals | `data/upstream/final_tables/Supp_Table_2_AC_master.csv` | Final H3K27ac catalog |
| Developmental H3K27me3 intervals | `data/upstream/final_tables/Supp_Table_3_ME3_master.csv` | Final H3K27me3 catalog |
| Original WGCNA raw counts | `data/upstream/WGCNA/WGCNA_input_data/atac_idr_peaks_merged_uniform-mm10.tsv`, `hac_overlap_peaks_merged-mm10.tsv`, `htk_overlap_peaks_merged-mm10.tsv` | Original peak-by-sample count tables for optional WGCNA fitting; copied from the project WGCNA_input_data folder |
| Merged histone input checkpoints | `data/upstream/WGCNA/WGCNA_input_data/data_h3k27ac_merged_by_1kb.RData`, `data_h3k27me3_merged_by_1kb.RData` | Archived intermediate count objects; not replacements for the raw TSVs without checking their schemas |
| Fitted WGCNA networks | `data/upstream/WGCNA/WGCNA_objects/*.Rds` | Three saved networks loaded by notebook 01 |
| Bulk RNA featureCounts | `data/upstream/RNA/featureCounts_final_15samples.txt` | Original counts and feature lengths for 15 libraries; retained provenance input, not a new DE model fit |
| Bulk RNA sample sheet | `data/upstream/RNA/sample_metadata.csv` | Library IDs, genotype, litter, sex and FlashTag population |
| Bulk RNA DE, VZ | `data/upstream/RNA/VZ_wt_vs_mut_DE_Litter_Sex_cor.csv` | Original model result table |
| Bulk RNA DE, nonVZ | `data/upstream/RNA/nonVZ_wt_vs_mut_DE_Litter_Sex_cor.csv` | Original model result table |
| Multiome pseudobulk counts | `data/upstream/multiome/RNA_pseudobulk_counts.csv`, `intervals_pseudobulk_counts.csv` | Aggregated count matrices for optional VST/PCA |
| Multiome module definitions | `data/upstream/multiome/atac_modules.csv` | Interval-to-module mapping |
| AC / ME3 peak-to-gene links | `data/upstream/multiome/AC.linked_peaks.csv`, `ME3.linked_peaks.csv` | Saved assay-specific links |
| DiffBind samples | `data/upstream/diffbind/sample_metadata.csv` | BAM basenames and contrast metadata; BAM files themselves are external |
| Differential peak tables | `data/upstream/final_tables/*Flashtag*DBA_DESEQ2.csv` | Saved ATAC/AC/ME3 statistics by FlashTag population |
| Numerical figure tables | `data/source_tables/<figure>/*.csv` | Clean exports of manuscript source-data worksheets |
| Histology measurement grids | `data/raw_sheets/<figure>/*.csv` | Original irregular worksheet grids parsed by `R/histology.R` |

## Find every file and saved object

- [Per-notebook input declarations](config/notebook_inputs.csv): saved-results versus optional-rerun collections.
- [Expanded per-notebook file list](config/notebook_input_files.csv): individual filenames, byte sizes and purpose.
- [Complete file inventory](config/file_inventory.csv): every local data/object file, original source where recorded and declared consumers.
- [Saved R object inventory](SAVED_OBJECTS.md): every RDS/RData/RDA checkpoint, including retained/generated objects that are not directly loaded.

Collection declarations are conservative: a figure's table directory includes
more tables than any single panel reads. Generated caches are not automatically
treated as required inputs. Recorded copy hashes are inherited provenance, not a
claim that multi-gigabyte files were rehashed during this documentation update.

## Execution

After supplying the archived inputs at these paths, run from `CortexCode/`:

```sh
Rscript setup.R
Rscript check_inputs.R 03
Rscript run_all.R 03
```

Omit `03` to check or render all main notebooks. `check_inputs.R` checks names and
recorded sizes without loading R objects. Add `--include-optional` to inspect
declared optional inputs; missing external objects will then be reported. This is
an input completeness check, not a scientific rerun or a full checksum audit.
The output is `reports/input_check.csv`.

For a supplementary notebook, render from the same directory, for example:

```r
rmarkdown::render("pipeline/03_multiome_module_scores.Rmd", knit_root_dir = getwd())
```

The main notebooks contain the analysis implementations and figure exports.
`pipeline/` provides additional upstream documentation and selected checkpoint
inspection. Saved-result inputs are distinct from optional original analysis
inputs; turning a rerun flag on does not make missing raw files available.

## External inputs for optional analyses

| Setting | Needed for |
|---|---|
| `CORTEX_BAM_DIR` | Original BAMs and indices for DiffBind |
| `CORTEX_MULTIOME_OBJECT` | Full original Seurat object for extracting marker layers; default `data/upstream/multiome/Ctx_dev.clean.rds` is not bundled |
| `CORTEX_VARIANTS_TSV` | Original human variant TSV |
| `CORTEX_LIFTOVER_CHAIN` | hg38-to-mm10 chain |
| `CORTEX_TOPGO_CONFIG` | Explicit topGO configuration; original Figure 5d settings were not recovered |
| `CORTEX_HOMER_BIN` | HOMER executable if absent from PATH; mm10 installation also required |

See [DATA_AVAILABILITY.md](DATA_AVAILABILITY.md) for methodological limitations
and [environment/](environment/) for the captured reconstruction environment.
