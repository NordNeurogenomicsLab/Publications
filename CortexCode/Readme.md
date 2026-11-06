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

### Primary RNA-seq and chromatin count tables

| File name | Contents | GitHub repository path |
|---|---|---|
| `featureCounts_final_15samples.txt` | Bulk RNA-seq gene counts and feature lengths for 15 libraries | **Not deposited**; expected path: `data/upstream/RNA/featureCounts_final_15samples.txt` |
| `hac_overlap_peaks_merged-mm10.tsv` | H3K27ac ChIP-seq peak-by-sample counts | [data/upstream/WGCNA/WGCNA_input_data/hac_overlap_peaks_merged-mm10.tsv](data/upstream/WGCNA/WGCNA_input_data/hac_overlap_peaks_merged-mm10.tsv) |
| `htk_overlap_peaks_merged-mm10.tsv` | H3K27me3 ChIP-seq peak-by-sample counts | [data/upstream/WGCNA/WGCNA_input_data/htk_overlap_peaks_merged-mm10.tsv](data/upstream/WGCNA/WGCNA_input_data/htk_overlap_peaks_merged-mm10.tsv) |
| `atac_idr_peaks_merged_uniform-mm10.tsv` | ATAC-seq peak-by-sample counts | [data/upstream/WGCNA/WGCNA_input_data/atac_idr_peaks_merged_uniform-mm10.tsv](data/upstream/WGCNA/WGCNA_input_data/atac_idr_peaks_merged_uniform-mm10.tsv) |

### Metadata

| File name | Contents | GitHub repository path |
|---|---|---|
| `sample_metadata.csv` | Bulk RNA library IDs, genotype, litter, sex and FlashTag population | [data/upstream/RNA/sample_metadata.csv](data/upstream/RNA/sample_metadata.csv) |
| `sample_metadata.csv` | ATAC/ChIP-seq sample metadata and BAM basenames for DiffBind | **Not deposited**; expected path: `data/upstream/diffbind/sample_metadata.csv` |
| `multiome_cell_metadata.csv` | Multiome cell identities and annotations | **Not deposited**; expected path: `data/upstream/multiome/multiome_cell_metadata.csv` |

### WGCNA annotations and integrated data

| File name | Contents | GitHub repository path |
|---|---|---|
| `Module_renaming_dictionary.csv` | Original module labels mapped to manuscript names | [data/upstream/Module_renaming_dictionary.csv](data/upstream/Module_renaming_dictionary.csv) |
| `atac_modules.csv` | ATAC interval-to-module assignments | [data/upstream/multiome/atac_modules.csv](data/upstream/multiome/atac_modules.csv) |
| `Supp_Table_1_ATAC_master.csv` | Final annotated ATAC intervals and modules | **Not deposited**; expected path: `data/upstream/final_tables/Supp_Table_1_ATAC_master.csv` |
| `Supp_Table_2_AC_master.csv` | Final annotated H3K27ac intervals and modules | **Not deposited**; expected path: `data/upstream/final_tables/Supp_Table_2_AC_master.csv` |
| `Supp_Table_3_ME3_master.csv` | Final annotated H3K27me3 intervals and modules | **Not deposited**; expected path: `data/upstream/final_tables/Supp_Table_3_ME3_master.csv` |
| `Supp_Table_4_WGCNA_interactions_and_Programs.csv` | Integrated WGCNA interval interactions and developmental programs | **Not deposited**; expected path: `data/upstream/final_tables/Supp_Table_4_WGCNA_interactions_and_Programs.csv` |
| `all_overlaps_total_aug_with_Programs.csv` | Linked ATAC/AC/ME3 interval catalog with program assignments | **Not deposited**; expected path: `data/upstream/integration/all_overlaps_total_aug_with_Programs.csv` |
| `AC.linked_peaks.csv` | H3K27ac peak-to-gene links and statistics | [data/upstream/multiome/AC.linked_peaks.csv](data/upstream/multiome/AC.linked_peaks.csv) |
| `ME3.linked_peaks.csv` | H3K27me3 peak-to-gene links and statistics | [data/upstream/multiome/ME3.linked_peaks.csv](data/upstream/multiome/ME3.linked_peaks.csv) |
| `network_nodes.csv` | Integrated module-network node annotations | [data/upstream/network_nodes.csv](data/upstream/network_nodes.csv) |
| `network_edges.csv` | Integrated module-network connections | [data/upstream/network_edges.csv](data/upstream/network_edges.csv) |
| `4f_Program_counts.csv` | Counts assigned to each integrated developmental program | [data/source_tables/Fig4_KC/4f_Program_counts.csv](data/source_tables/Fig4_KC/4f_Program_counts.csv) |
| `4f_Quantized_trajectories.csv` | Quantized ATAC/AC/ME3 program trajectories | [data/source_tables/Fig4_KC/4f_Quantized_trajectories.csv](data/source_tables/Fig4_KC/4f_Quantized_trajectories.csv) |

### DiffBind and RNA differential-expression results

| File name | Contents | GitHub repository path |
|---|---|---|
| `Ac_FTneg_Modality_H3K27Ac_PeakSet_ac_expand_bed_5000_Condition_FlashtagNeg__genotype_DBA_DESEQ2.csv` | H3K27ac DiffBind/DESeq2 results, FlashTag-negative cells | **Not deposited**; expected path: `data/upstream/final_tables/Ac_FTneg_Modality_H3K27Ac_PeakSet_ac_expand_bed_5000_Condition_FlashtagNeg__genotype_DBA_DESEQ2.csv` |
| `Ac_FTpos_Modality_H3K27Ac_PeakSet_ac_expand_bed_5000_Condition_FlashtagPos__genotype_DBA_DESEQ2.csv` | H3K27ac DiffBind/DESeq2 results, FlashTag-positive cells | **Not deposited**; expected path: `data/upstream/final_tables/Ac_FTpos_Modality_H3K27Ac_PeakSet_ac_expand_bed_5000_Condition_FlashtagPos__genotype_DBA_DESEQ2.csv` |
| `ATAC_FTneg_Modality_ATAC_PeakSet_atac_expand_bed_5000_Condition_FlashtagNeg__genotype_DBA_DESEQ2.csv` | ATAC DiffBind/DESeq2 results, FlashTag-negative cells | **Not deposited**; expected path: `data/upstream/final_tables/ATAC_FTneg_Modality_ATAC_PeakSet_atac_expand_bed_5000_Condition_FlashtagNeg__genotype_DBA_DESEQ2.csv` |
| `ATAC_FTpos_Modality_ATAC_PeakSet_atac_expand_bed_5000_Condition_FlashtagPos__genotype_DBA_DESEQ2.csv` | ATAC DiffBind/DESeq2 results, FlashTag-positive cells | **Not deposited**; expected path: `data/upstream/final_tables/ATAC_FTpos_Modality_ATAC_PeakSet_atac_expand_bed_5000_Condition_FlashtagPos__genotype_DBA_DESEQ2.csv` |
| `H3K27Me3_FTneg_Modality_H3K27Me3_PeakSet_me3_expand_bed_5000_Condition_FlashtagNeg__genotype_DBA_DESEQ2.csv` | H3K27me3 DiffBind/DESeq2 results, FlashTag-negative cells | **Not deposited**; expected path: `data/upstream/final_tables/H3K27Me3_FTneg_Modality_H3K27Me3_PeakSet_me3_expand_bed_5000_Condition_FlashtagNeg__genotype_DBA_DESEQ2.csv` |
| `H3K27Me3_FTpos_Modality_H3K27Me3_PeakSet_me3_expand_bed_5000_Condition_FlashtagPos__genotype_DBA_DESEQ2.csv` | H3K27me3 DiffBind/DESeq2 results, FlashTag-positive cells | **Not deposited**; expected path: `data/upstream/final_tables/H3K27Me3_FTpos_Modality_H3K27Me3_PeakSet_me3_expand_bed_5000_Condition_FlashtagPos__genotype_DBA_DESEQ2.csv` |
| `VZ_wt_vs_mut_DE_Litter_Sex_cor.csv` | VZ RNA-seq differential expression, adjusted for litter and sex | [data/upstream/RNA/VZ_wt_vs_mut_DE_Litter_Sex_cor.csv](data/upstream/RNA/VZ_wt_vs_mut_DE_Litter_Sex_cor.csv) |
| `nonVZ_wt_vs_mut_DE_Litter_Sex_cor.csv` | nonVZ RNA-seq differential expression, adjusted for litter and sex | [data/upstream/RNA/nonVZ_wt_vs_mut_DE_Litter_Sex_cor.csv](data/upstream/RNA/nonVZ_wt_vs_mut_DE_Litter_Sex_cor.csv) |

Other figure-level tables, IHC/thickness measurements and saved-object references: [Source_tables](SOURCE_TABLES.md).

Shared R functions: [R/](R/).
