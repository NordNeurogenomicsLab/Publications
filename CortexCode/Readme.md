## Manuscript

**Profiling dynamic epigenetic programs in the developing mouse cortex reveals the role of ASD-associated Kdm6b in deep-layer neuron differentiation.**

Athéna R. Ypsilanti<sup>†,*</sup>, Eirene Markenscoff-Papadimitriou<sup>†,*</sup>, Karol Cichewicz<sup>†</sup>, Sean Whalen, Nicolas Seban, Ethan Fenton, Kenneth Lim, Shan Dong, Jada Baldoza, Yurong Wang, Matt Barraza, Fadya Binyameen, Diane Dickel, Len A. Pennacchio, Axel Visel, Stephan J. Sanders, Katherine S. Pollard, Alex S. Nord<sup>*</sup>, John Rubenstein<sup>*</sup>

† Contributed equally: Athéna R. Ypsilanti, Eirene Markenscoff-Papadimitriou and Karol Cichewicz.

\* Corresponding authors: Athéna R. Ypsilanti, Eirene Markenscoff-Papadimitriou, Alex S. Nord and John Rubenstein.

## Inputs and environment

[Source_tables](SOURCE_TABLES.md) · [Tables on GitHub](data/TABLES_ON_GITHUB.csv) · [Input files and run order](INPUTS.md) · [Saved RDS/RData objects](SAVED_OBJECTS.md) · [R session and package versions](environment/README.md)

## Figure and analysis notebooks

| Notebook | Contents |
|---|---|
| [00_inputs_and_provenance.Rmd](00_inputs_and_provenance.Rmd) | Inputs and provenance |
| [01_WGCNA_checkpoints.Rmd](01_WGCNA_checkpoints.Rmd) | WGCNA networks and module refinement |
| [02_Figures_1_2_modules.Rmd](02_Figures_1_2_modules.Rmd) | Figures 1–2: modules, GO enrichment and motifs |
| [03_Figure_3_multiome.Rmd](03_Figure_3_multiome.Rmd) | Figure 3: multiome analysis |
| [04_Figure_4_integration.Rmd](04_Figure_4_integration.Rmd) | Figure 4: chromatin integration and trajectories |
| [05_Figure_5_RNA_histology.Rmd](05_Figure_5_RNA_histology.Rmd) | Figure 5: RNA, topGO and histology |
| [06_Figure_6_epigenomics.Rmd](06_Figure_6_epigenomics.Rmd) | Figure 6: differential chromatin analysis |
| [07_Figure_6_variants.Rmd](07_Figure_6_variants.Rmd) | Figure 6: variant enrichment |
| [08_Figure_7_enhancer_histology.Rmd](08_Figure_7_enhancer_histology.Rmd) | Figure 7: enhancer reporters and histology |
| [09_manuscript_graphics.Rmd](09_manuscript_graphics.Rmd) | Panel styling and exports |
| [10_additional_manuscript_panels.Rmd](10_additional_manuscript_panels.Rmd) | Additional computational panels |
| [11_Extended_multiome_heatmaps.Rmd](11_Extended_multiome_heatmaps.Rmd) | Extended Data multiome heatmaps |

## Upstream pipeline notebooks

| Notebook | Contents |
|---|---|
| [00_sequencing_inputs.Rmd](pipeline/00_sequencing_inputs.Rmd) | Sequencing inputs and processing records |
| [01_WGCNA_network_generation.Rmd](pipeline/01_WGCNA_network_generation.Rmd) | WGCNA network generation |
| [02_module_annotation.Rmd](pipeline/02_module_annotation.Rmd) | Module refinement and manuscript nomenclature |
| [03_multiome_module_scores.Rmd](pipeline/03_multiome_module_scores.Rmd) | Multiome accessibility and linked-gene scores |
| [04_integration_trajectories.Rmd](pipeline/04_integration_trajectories.Rmd) | Chromatin integration and developmental trajectories |
| [05_RNA_and_pseudobulk.Rmd](pipeline/05_RNA_and_pseudobulk.Rmd) | RNA results and pseudobulk PCA |
| [06_epigenomic_statistics.Rmd](pipeline/06_epigenomic_statistics.Rmd) | Differential chromatin statistics and enrichment |
| [07_variant_enrichment.Rmd](pipeline/07_variant_enrichment.Rmd) | Balanced variant enrichment |

## Input data tables

[Box data and saved objects](https://ucdavis.box.com/s/jrmg33n2dcljxadenvx47ji8c6azg02s) (308 files, 5.2 GB). Download the release and merge its `data/` and `objects/` folders into this code directory. [Paths and checksums](config/BOX_RELEASE_FILES.csv).

### Primary RNA-seq and chromatin count tables

| File name | Contents | Repository |
|---|---|---|
| <code>featureCounts_final_</code><br><code>15samples.txt</code> | Kdm6b cKO and WT bulk RNA-seq gene counts and feature lengths (15 libraries) | [Box](https://ucdavis.app.box.com/s/jrmg33n2dcljxadenvx47ji8c6azg02s/file/2575320512459 "data/upstream/RNA/featureCounts_final_15samples.txt") |
| <code>hac_overlap_peaks_merged-</code><br><code>mm10.tsv</code> | Developmental H3K27ac ChIP-seq peak-by-sample counts | [File](data/upstream/WGCNA/WGCNA_input_data/hac_overlap_peaks_merged-mm10.tsv "data/upstream/WGCNA/WGCNA_input_data/hac_overlap_peaks_merged-mm10.tsv") |
| <code>htk_overlap_peaks_merged-</code><br><code>mm10.tsv</code> | Developmental H3K27me3 ChIP-seq peak-by-sample counts | [File](data/upstream/WGCNA/WGCNA_input_data/htk_overlap_peaks_merged-mm10.tsv "data/upstream/WGCNA/WGCNA_input_data/htk_overlap_peaks_merged-mm10.tsv") |
| <code>atac_idr_peaks_merged_</code><br><code>uniform-mm10.tsv</code> | Developmental ATAC-seq peak-by-sample counts | [File](data/upstream/WGCNA/WGCNA_input_data/atac_idr_peaks_merged_uniform-mm10.tsv "data/upstream/WGCNA/WGCNA_input_data/atac_idr_peaks_merged_uniform-mm10.tsv") |

Kdm6b cKO chromatin and RNA differential results: [ATAC/ChIP-seq and RNA-seq results](#kdm6b-cko-diffbind-and-rna-seq-results).

### Metadata

| File name | Contents | Repository |
|---|---|---|
| `sample_metadata.csv` | Bulk RNA library IDs, genotype, litter, sex and FlashTag population | [File](data/upstream/RNA/sample_metadata.csv "data/upstream/RNA/sample_metadata.csv") |
| `sample_metadata.csv` | ATAC/ChIP-seq sample metadata and BAM basenames for DiffBind | [Box](https://ucdavis.app.box.com/s/jrmg33n2dcljxadenvx47ji8c6azg02s/file/2575320687685 "data/upstream/diffbind/sample_metadata.csv") |
| `multiome_cell_metadata.csv` | Multiome cell identities and annotations | [Box](https://ucdavis.app.box.com/s/jrmg33n2dcljxadenvx47ji8c6azg02s/file/2575321679668 "data/upstream/multiome/multiome_cell_metadata.csv") |

### WGCNA annotations and integrated data

| File name | Contents | Repository |
|---|---|---|
| <code>Module_renaming_</code><br><code>dictionary.csv</code> | Original module labels mapped to manuscript names | [File](data/upstream/Module_renaming_dictionary.csv "data/upstream/Module_renaming_dictionary.csv") |
| `atac_modules.csv` | ATAC interval-to-module assignments | [File](data/upstream/multiome/atac_modules.csv "data/upstream/multiome/atac_modules.csv") |
| `Supp_Table_1_ATAC_master.csv` | Final annotated ATAC intervals and modules | [Box](https://ucdavis.app.box.com/s/jrmg33n2dcljxadenvx47ji8c6azg02s/file/2575313834073 "data/upstream/final_tables/Supp_Table_1_ATAC_master.csv") |
| `Supp_Table_2_AC_master.csv` | Final annotated H3K27ac intervals and modules | [Box](https://ucdavis.app.box.com/s/jrmg33n2dcljxadenvx47ji8c6azg02s/file/2575321878803 "data/upstream/final_tables/Supp_Table_2_AC_master.csv") |
| `Supp_Table_3_ME3_master.csv` | Final annotated H3K27me3 intervals and modules | [Box](https://ucdavis.app.box.com/s/jrmg33n2dcljxadenvx47ji8c6azg02s/file/2575315736751 "data/upstream/final_tables/Supp_Table_3_ME3_master.csv") |
| <code>Supp_Table_4_WGCNA_</code><br><code>interactions_and_</code><br><code>Programs.csv</code> | Integrated WGCNA interval interactions and developmental programs | [Box](https://ucdavis.app.box.com/s/jrmg33n2dcljxadenvx47ji8c6azg02s/file/2575312699106 "data/upstream/final_tables/Supp_Table_4_WGCNA_interactions_and_Programs.csv") |
| <code>all_overlaps_total_aug_with_</code><br><code>Programs.csv</code> | Linked ATAC/AC/ME3 interval catalog with program assignments | [Box](https://ucdavis.app.box.com/s/jrmg33n2dcljxadenvx47ji8c6azg02s/file/2575313800586 "data/upstream/integration/all_overlaps_total_aug_with_Programs.csv") |
| `AC.linked_peaks.csv` | H3K27ac peak-to-gene links and statistics | [File](data/upstream/multiome/AC.linked_peaks.csv "data/upstream/multiome/AC.linked_peaks.csv") |
| `ME3.linked_peaks.csv` | H3K27me3 peak-to-gene links and statistics | [File](data/upstream/multiome/ME3.linked_peaks.csv "data/upstream/multiome/ME3.linked_peaks.csv") |
| `network_nodes.csv` | Integrated module-network node annotations | [File](data/upstream/network_nodes.csv "data/upstream/network_nodes.csv") |
| `network_edges.csv` | Integrated module-network connections | [File](data/upstream/network_edges.csv "data/upstream/network_edges.csv") |
| `4f_Program_counts.csv` | Counts assigned to each integrated developmental program | [File](data/source_tables/Fig4_KC/4f_Program_counts.csv "data/source_tables/Fig4_KC/4f_Program_counts.csv") |
| <code>4f_Quantized_</code><br><code>trajectories.csv</code> | Quantized ATAC/AC/ME3 program trajectories | [File](data/source_tables/Fig4_KC/4f_Quantized_trajectories.csv "data/source_tables/Fig4_KC/4f_Quantized_trajectories.csv") |

### Kdm6b cKO DiffBind and RNA-seq results

These are derived statistical results; the chromatin tables also retain per-sample signal columns. They are not raw read-count matrices.

| File name | Contents | Repository |
|---|---|---|
| <code>Ac_FTneg_Modality_H3K27Ac_</code><br><code>PeakSet_ac_expand_bed_5000_</code><br><code>Condition_FlashtagNeg__</code><br><code>genotype_DBA_DESEQ2.csv</code> | H3K27ac Kdm6b cKO versus WT DiffBind/DESeq2 results, FlashTag-negative cells | [Box](https://ucdavis.app.box.com/s/jrmg33n2dcljxadenvx47ji8c6azg02s/file/2575319214240 "data/upstream/final_tables/Ac_FTneg_Modality_H3K27Ac_PeakSet_ac_expand_bed_5000_Condition_FlashtagNeg__genotype_DBA_DESEQ2.csv") |
| <code>Ac_FTpos_Modality_H3K27Ac_</code><br><code>PeakSet_ac_expand_bed_5000_</code><br><code>Condition_FlashtagPos__</code><br><code>genotype_DBA_DESEQ2.csv</code> | H3K27ac Kdm6b cKO versus WT DiffBind/DESeq2 results, FlashTag-positive cells | [Box](https://ucdavis.app.box.com/s/jrmg33n2dcljxadenvx47ji8c6azg02s/file/2575320037788 "data/upstream/final_tables/Ac_FTpos_Modality_H3K27Ac_PeakSet_ac_expand_bed_5000_Condition_FlashtagPos__genotype_DBA_DESEQ2.csv") |
| <code>ATAC_FTneg_Modality_ATAC_</code><br><code>PeakSet_atac_expand_bed_</code><br><code>5000_Condition_FlashtagNeg__</code><br><code>genotype_DBA_DESEQ2.csv</code> | ATAC Kdm6b cKO versus WT DiffBind/DESeq2 results, FlashTag-negative cells | [Box](https://ucdavis.app.box.com/s/jrmg33n2dcljxadenvx47ji8c6azg02s/file/2575311270143 "data/upstream/final_tables/ATAC_FTneg_Modality_ATAC_PeakSet_atac_expand_bed_5000_Condition_FlashtagNeg__genotype_DBA_DESEQ2.csv") |
| <code>ATAC_FTpos_Modality_ATAC_</code><br><code>PeakSet_atac_expand_bed_</code><br><code>5000_Condition_FlashtagPos__</code><br><code>genotype_DBA_DESEQ2.csv</code> | ATAC Kdm6b cKO versus WT DiffBind/DESeq2 results, FlashTag-positive cells | [Box](https://ucdavis.app.box.com/s/jrmg33n2dcljxadenvx47ji8c6azg02s/file/2575321875790 "data/upstream/final_tables/ATAC_FTpos_Modality_ATAC_PeakSet_atac_expand_bed_5000_Condition_FlashtagPos__genotype_DBA_DESEQ2.csv") |
| <code>H3K27Me3_FTneg_Modality_</code><br><code>H3K27Me3_PeakSet_me3_expand_</code><br><code>bed_5000_Condition_</code><br><code>FlashtagNeg__genotype_DBA_</code><br><code>DESEQ2.csv</code> | H3K27me3 Kdm6b cKO versus WT DiffBind/DESeq2 results, FlashTag-negative cells | [Box](https://ucdavis.app.box.com/s/jrmg33n2dcljxadenvx47ji8c6azg02s/file/2575311414133 "data/upstream/final_tables/H3K27Me3_FTneg_Modality_H3K27Me3_PeakSet_me3_expand_bed_5000_Condition_FlashtagNeg__genotype_DBA_DESEQ2.csv") |
| <code>H3K27Me3_FTpos_Modality_</code><br><code>H3K27Me3_PeakSet_me3_expand_</code><br><code>bed_5000_Condition_</code><br><code>FlashtagPos__genotype_DBA_</code><br><code>DESEQ2.csv</code> | H3K27me3 Kdm6b cKO versus WT DiffBind/DESeq2 results, FlashTag-positive cells | [Box](https://ucdavis.app.box.com/s/jrmg33n2dcljxadenvx47ji8c6azg02s/file/2575311219711 "data/upstream/final_tables/H3K27Me3_FTpos_Modality_H3K27Me3_PeakSet_me3_expand_bed_5000_Condition_FlashtagPos__genotype_DBA_DESEQ2.csv") |
| <code>VZ_wt_vs_mut_DE_Litter_Sex_</code><br><code>cor.csv</code> | VZ Kdm6b cKO versus WT RNA-seq differential expression, adjusted for litter and sex | [File](data/upstream/RNA/VZ_wt_vs_mut_DE_Litter_Sex_cor.csv "data/upstream/RNA/VZ_wt_vs_mut_DE_Litter_Sex_cor.csv") |
| <code>nonVZ_wt_vs_mut_DE_Litter_</code><br><code>Sex_cor.csv</code> | nonVZ Kdm6b cKO versus WT RNA-seq differential expression, adjusted for litter and sex | [File](data/upstream/RNA/nonVZ_wt_vs_mut_DE_Litter_Sex_cor.csv "data/upstream/RNA/nonVZ_wt_vs_mut_DE_Litter_Sex_cor.csv") |

Other figure-level tables, IHC/thickness measurements and saved-object references: [Source_tables](SOURCE_TABLES.md).

Shared R functions: [R/](R/).
