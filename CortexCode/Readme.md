## Manuscript

**Profiling dynamic epigenetic programs in the developing mouse cortex reveals the role of ASD-associated Kdm6b in deep-layer neuron differentiation.**

Athéna R. Ypsilanti<sup>†,*</sup>, Eirene Markenscoff-Papadimitriou<sup>†,*</sup>, Karol Cichewicz<sup>†</sup>, Sean Whalen, Nicolas Seban, Ethan Fenton, Kenneth Lim, Shan Dong, Jada Baldoza, Yurong Wang, Matt Barraza, Fadya Binyameen, Diane Dickel, Len A. Pennacchio, Axel Visel, Stephan J. Sanders, Katherine S. Pollard, Alex S. Nord<sup>*</sup>, John Rubenstein<sup>*</sup>

† Contributed equally: Athéna R. Ypsilanti, Eirene Markenscoff-Papadimitriou and Karol Cichewicz.

\* Corresponding authors: Athéna R. Ypsilanti, Eirene Markenscoff-Papadimitriou, Alex S. Nord and John Rubenstein.

## Inputs and environment

[Tables on GitHub](data/TABLES_ON_GITHUB.csv) · [Input files and run order](INPUTS.md) · [Saved RDS/RData objects](SAVED_OBJECTS.md) · [R session and package versions](environment/README.md)

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

### Primary count tables, measurements and sample metadata

Count tables are assay-level inputs; sequencing reads and BAMs are not included.

| File name | Contents | GitHub repository path |
|---|---|---|
| `Ext_7j_Tle4_thickness.csv` | Original TLE4 layer-thickness measurements | [data/raw_sheets/Extended_Fig7_KC/Ext_7j_Tle4_thickness.csv](data/raw_sheets/Extended_Fig7_KC/Ext_7j_Tle4_thickness.csv) |
| `Ext_7kTbr1_thickness.csv` | Original TBR1 layer-thickness measurements | [data/raw_sheets/Extended_Fig7_KC/Ext_7kTbr1_thickness.csv](data/raw_sheets/Extended_Fig7_KC/Ext_7kTbr1_thickness.csv) |
| `9f_Tbr1_Nkx2-1.csv` | Original TBR1/Nkx2-1 intensity measurements and supplied normalization | [data/raw_sheets/Extended_Fig9_KC/9f_Tbr1_Nkx2-1.csv](data/raw_sheets/Extended_Fig9_KC/9f_Tbr1_Nkx2-1.csv) |
| `5e_TLE4_E14.csv` | Original E14 TLE4 histology measurement grid | [data/raw_sheets/Fig5_KC/5e_TLE4_E14.csv](data/raw_sheets/Fig5_KC/5e_TLE4_E14.csv) |
| `5f_Tbr1_E14.csv` | Original E14 TBR1 histology measurement grid | [data/raw_sheets/Fig5_KC/5f_Tbr1_E14.csv](data/raw_sheets/Fig5_KC/5f_Tbr1_E14.csv) |
| `5g_Tbr1_P10.csv` | Original P10 TBR1 histology measurement grid | [data/raw_sheets/Fig5_KC/5g_Tbr1_P10.csv](data/raw_sheets/Fig5_KC/5g_Tbr1_P10.csv) |
| `5h_Bcl11b_P10.csv` | Original P10 BCL11B histology measurement grid | [data/raw_sheets/Fig5_KC/5h_Bcl11b_P10.csv](data/raw_sheets/Fig5_KC/5h_Bcl11b_P10.csv) |
| `5i_Fezf2_P10.csv` | Original P10 FEZF2 histology measurement grid | [data/raw_sheets/Fig5_KC/5i_Fezf2_P10.csv](data/raw_sheets/Fig5_KC/5i_Fezf2_P10.csv) |
| `7k_l_Tbr1_416Del.csv` | Original enhancer-deletion TBR1 intensity measurements and supplied normalization | [data/raw_sheets/Fig7_KC/7k_l_Tbr1_416Del.csv](data/raw_sheets/Fig7_KC/7k_l_Tbr1_416Del.csv) |
| `sample_metadata.csv` | Bulk RNA library IDs, genotype, litter, sex and FlashTag population | [data/upstream/RNA/sample_metadata.csv](data/upstream/RNA/sample_metadata.csv) |
| `atac_idr_peaks_merged_uniform-mm10.tsv` | Developmental ATAC peak-by-sample counts used as WGCNA input | [data/upstream/WGCNA/WGCNA_input_data/atac_idr_peaks_merged_uniform-mm10.tsv](data/upstream/WGCNA/WGCNA_input_data/atac_idr_peaks_merged_uniform-mm10.tsv) |
| `hac_overlap_peaks_merged-mm10.tsv` | Developmental H3K27ac peak-by-sample counts used as WGCNA input | [data/upstream/WGCNA/WGCNA_input_data/hac_overlap_peaks_merged-mm10.tsv](data/upstream/WGCNA/WGCNA_input_data/hac_overlap_peaks_merged-mm10.tsv) |
| `htk_overlap_peaks_merged-mm10.tsv` | Developmental H3K27me3 peak-by-sample counts used as WGCNA input | [data/upstream/WGCNA/WGCNA_input_data/htk_overlap_peaks_merged-mm10.tsv](data/upstream/WGCNA/WGCNA_input_data/htk_overlap_peaks_merged-mm10.tsv) |
| `featureCounts_final_15samples.txt` | Original bulk RNA gene counts and feature lengths for 15 libraries | **Not deposited**; expected path: `data/upstream/RNA/featureCounts_final_15samples.txt` |

### Derived analysis and plotting tables

| File name | Contents | GitHub repository path |
|---|---|---|
| `VZ_E12_vs_ImN_E13.csv` | Original footprinting contrast worksheet: VZ E12 vs ImN E13 | [data/raw_sheets/Extended_Fig2_KC/VZ_E12_vs_ImN_E13.csv](data/raw_sheets/Extended_Fig2_KC/VZ_E12_vs_ImN_E13.csv) |
| `VZ_E12_vs_L5N_P2.csv` | Original footprinting contrast worksheet: VZ E12 vs L5N P2 | [data/raw_sheets/Extended_Fig2_KC/VZ_E12_vs_L5N_P2.csv](data/raw_sheets/Extended_Fig2_KC/VZ_E12_vs_L5N_P2.csv) |
| `VZ_E12_vs_L6N_P2.csv` | Original footprinting contrast worksheet: VZ E12 vs L6N P2 | [data/raw_sheets/Extended_Fig2_KC/VZ_E12_vs_L6N_P2.csv](data/raw_sheets/Extended_Fig2_KC/VZ_E12_vs_L6N_P2.csv) |
| `VZ_E12_vs_SVZ_E12.csv` | Original footprinting contrast worksheet: VZ E12 vs SVZ E12 | [data/raw_sheets/Extended_Fig2_KC/VZ_E12_vs_SVZ_E12.csv](data/raw_sheets/Extended_Fig2_KC/VZ_E12_vs_SVZ_E12.csv) |
| `1e_DiffBind_counts.csv` | Extended Data 1: DiffBind counts | [data/source_tables/Extended_Fig1_KC/1e_DiffBind_counts.csv](data/source_tables/Extended_Fig1_KC/1e_DiffBind_counts.csv) |
| `1f_DiffBind_counts.csv` | Extended Data 1: DiffBind counts | [data/source_tables/Extended_Fig1_KC/1f_DiffBind_counts.csv](data/source_tables/Extended_Fig1_KC/1f_DiffBind_counts.csv) |
| `1g_DiffBind_counts.csv` | Extended Data 1: DiffBind counts | [data/source_tables/Extended_Fig1_KC/1g_DiffBind_counts.csv](data/source_tables/Extended_Fig1_KC/1g_DiffBind_counts.csv) |
| `1g_Significant_intervals.csv` | Extended Data 1: Significant intervals | [data/source_tables/Extended_Fig1_KC/1g_Significant_intervals.csv](data/source_tables/Extended_Fig1_KC/1g_Significant_intervals.csv) |
| `VZ_E12_vs_ImN_E13.csv` | Footprinting statistics and derived display values: VZ E12 vs ImN E13 | [data/source_tables/Extended_Fig2_KC/VZ_E12_vs_ImN_E13.csv](data/source_tables/Extended_Fig2_KC/VZ_E12_vs_ImN_E13.csv) |
| `VZ_E12_vs_L5N_P2.csv` | Footprinting statistics and derived display values: VZ E12 vs L5N P2 | [data/source_tables/Extended_Fig2_KC/VZ_E12_vs_L5N_P2.csv](data/source_tables/Extended_Fig2_KC/VZ_E12_vs_L5N_P2.csv) |
| `VZ_E12_vs_L6N_P2.csv` | Footprinting statistics and derived display values: VZ E12 vs L6N P2 | [data/source_tables/Extended_Fig2_KC/VZ_E12_vs_L6N_P2.csv](data/source_tables/Extended_Fig2_KC/VZ_E12_vs_L6N_P2.csv) |
| `VZ_E12_vs_SVZ_E12.csv` | Footprinting statistics and derived display values: VZ E12 vs SVZ E12 | [data/source_tables/Extended_Fig2_KC/VZ_E12_vs_SVZ_E12.csv](data/source_tables/Extended_Fig2_KC/VZ_E12_vs_SVZ_E12.csv) |
| `4c_GENOA_distal_counts.csv` | Extended Data 4: GENOA distal counts | [data/source_tables/Extended_Fig4_KC/4c_GENOA_distal_counts.csv](data/source_tables/Extended_Fig4_KC/4c_GENOA_distal_counts.csv) |
| `5a_Per_cell_QC.csv` | Extended Data 5: Per cell QC | [data/source_tables/Extended_Fig5_KC/5a_Per_cell_QC.csv](data/source_tables/Extended_Fig5_KC/5a_Per_cell_QC.csv) |
| `5b_RNA_marker_expression.csv` | Extended Data 5: RNA marker expression | [data/source_tables/Extended_Fig5_KC/5b_RNA_marker_expression.csv](data/source_tables/Extended_Fig5_KC/5b_RNA_marker_expression.csv) |
| `5c_Bulk_median_CPM.csv` | Extended Data 5: Bulk median CPM | [data/source_tables/Extended_Fig5_KC/5c_Bulk_median_CPM.csv](data/source_tables/Extended_Fig5_KC/5c_Bulk_median_CPM.csv) |
| `5g_Marker_results.csv` | Extended Data 5: Marker results | [data/source_tables/Extended_Fig5_KC/5g_Marker_results.csv](data/source_tables/Extended_Fig5_KC/5g_Marker_results.csv) |
| `5h_PCA_coordinates.csv` | Extended Data 5: PCA coordinates | [data/source_tables/Extended_Fig5_KC/5h_PCA_coordinates.csv](data/source_tables/Extended_Fig5_KC/5h_PCA_coordinates.csv) |
| `5h_selected_counts.csv` | Extended Data 5: selected counts | [data/source_tables/Extended_Fig5_KC/5h_selected_counts.csv](data/source_tables/Extended_Fig5_KC/5h_selected_counts.csv) |
| `5h_variance.csv` | Extended Data 5: variance | [data/source_tables/Extended_Fig5_KC/5h_variance.csv](data/source_tables/Extended_Fig5_KC/5h_variance.csv) |
| `5h_VST_inputs.csv` | Extended Data 5: VST inputs | [data/source_tables/Extended_Fig5_KC/5h_VST_inputs.csv](data/source_tables/Extended_Fig5_KC/5h_VST_inputs.csv) |
| `5j_PCA_coordinates.csv` | Extended Data 5: PCA coordinates | [data/source_tables/Extended_Fig5_KC/5j_PCA_coordinates.csv](data/source_tables/Extended_Fig5_KC/5j_PCA_coordinates.csv) |
| `5j_selected_counts.csv` | Extended Data 5: selected counts | [data/source_tables/Extended_Fig5_KC/5j_selected_counts.csv](data/source_tables/Extended_Fig5_KC/5j_selected_counts.csv) |
| `5j_variance.csv` | Extended Data 5: variance | [data/source_tables/Extended_Fig5_KC/5j_variance.csv](data/source_tables/Extended_Fig5_KC/5j_variance.csv) |
| `5j_VST_inputs.csv` | Extended Data 5: VST inputs | [data/source_tables/Extended_Fig5_KC/5j_VST_inputs.csv](data/source_tables/Extended_Fig5_KC/5j_VST_inputs.csv) |
| `5l_DA_module_proportions.csv` | Positive-DA interval counts and pooled module proportions within each cluster | [data/source_tables/Extended_Fig5_KC/5l_DA_module_proportions.csv](data/source_tables/Extended_Fig5_KC/5l_DA_module_proportions.csv) |
| `5l_Positive_DA_module_proportions.csv` | Explicit positive-DA export; same validated summary as 5l_DA_module_proportions.csv | [data/source_tables/Extended_Fig5_KC/5l_Positive_DA_module_proportions.csv](data/source_tables/Extended_Fig5_KC/5l_Positive_DA_module_proportions.csv) |
| `6a_TMK3_ATAC_link_membership.csv` | Extended Data 6: TMK3 ATAC link membership | [data/source_tables/Extended_Fig6_KC/6a_TMK3_ATAC_link_membership.csv](data/source_tables/Extended_Fig6_KC/6a_TMK3_ATAC_link_membership.csv) |
| `6a_TMK3_linked_ATAC_means.csv` | Extended Data 6: TMK3 linked ATAC means | [data/source_tables/Extended_Fig6_KC/6a_TMK3_linked_ATAC_means.csv](data/source_tables/Extended_Fig6_KC/6a_TMK3_linked_ATAC_means.csv) |
| `6c_1a_Quantized_rows.csv` | Extended Data 6: 1a Quantized rows | [data/source_tables/Extended_Fig6_KC/6c_1a_Quantized_rows.csv](data/source_tables/Extended_Fig6_KC/6c_1a_Quantized_rows.csv) |
| `6c_2a_Quantized_rows.csv` | Extended Data 6: 2a Quantized rows | [data/source_tables/Extended_Fig6_KC/6c_2a_Quantized_rows.csv](data/source_tables/Extended_Fig6_KC/6c_2a_Quantized_rows.csv) |
| `6c_2b_Quantized_rows.csv` | Extended Data 6: 2b Quantized rows | [data/source_tables/Extended_Fig6_KC/6c_2b_Quantized_rows.csv](data/source_tables/Extended_Fig6_KC/6c_2b_Quantized_rows.csv) |
| `7ab_Cluster_RNA_summary.csv` | Extended Data 7: Cluster RNA summary | [data/source_tables/Extended_Fig7_KC/7ab_Cluster_RNA_summary.csv](data/source_tables/Extended_Fig7_KC/7ab_Cluster_RNA_summary.csv) |
| `7ab_Per_cell_RNA_expression.csv` | Extended Data 7: Per cell RNA expression | [data/source_tables/Extended_Fig7_KC/7ab_Per_cell_RNA_expression.csv](data/source_tables/Extended_Fig7_KC/7ab_Per_cell_RNA_expression.csv) |
| `7cde_jk_RPKM.csv` | Sample-level RPKM for the Extended Data 7 RNA panels | [data/source_tables/Extended_Fig7_KC/7cde_jk_RPKM.csv](data/source_tables/Extended_Fig7_KC/7cde_jk_RPKM.csv) |
| `8f_Overlap_counts.csv` | Extended Data 8: Overlap counts | [data/source_tables/Extended_Fig8_KC/8f_Overlap_counts.csv](data/source_tables/Extended_Fig8_KC/8f_Overlap_counts.csv) |
| `8f_Overlap_pairs.csv` | Extended Data 8: Overlap pairs | [data/source_tables/Extended_Fig8_KC/8f_Overlap_pairs.csv](data/source_tables/Extended_Fig8_KC/8f_Overlap_pairs.csv) |
| `8g_GO_plotted.csv` | Extended Data 8: GO plotted | [data/source_tables/Extended_Fig8_KC/8g_GO_plotted.csv](data/source_tables/Extended_Fig8_KC/8g_GO_plotted.csv) |
| `8h_Module_enrichment.csv` | Extended Data 8: Module enrichment | [data/source_tables/Extended_Fig8_KC/8h_Module_enrichment.csv](data/source_tables/Extended_Fig8_KC/8h_Module_enrichment.csv) |
| `8h_Selected_membership.csv` | Extended Data 8: Selected membership | [data/source_tables/Extended_Fig8_KC/8h_Selected_membership.csv](data/source_tables/Extended_Fig8_KC/8h_Selected_membership.csv) |
| `8i_Variant_intersection.csv` | Extended Data 8: Variant intersection | [data/source_tables/Extended_Fig8_KC/8i_Variant_intersection.csv](data/source_tables/Extended_Fig8_KC/8i_Variant_intersection.csv) |
| `c1_Module_counts.csv` | Figure 1: ATAC Module counts | [data/source_tables/Fig1_KC/c1_Module_counts.csv](data/source_tables/Fig1_KC/c1_Module_counts.csv) |
| `c2_constitutive_heatmap.csv` | Figure 1: ATAC constitutive heatmap | [data/source_tables/Fig1_KC/c2_constitutive_heatmap.csv](data/source_tables/Fig1_KC/c2_constitutive_heatmap.csv) |
| `c2_dynamic_heatmap.csv` | Figure 1: ATAC dynamic heatmap | [data/source_tables/Fig1_KC/c2_dynamic_heatmap.csv](data/source_tables/Fig1_KC/c2_dynamic_heatmap.csv) |
| `c2_Median_signal.csv` | Figure 1: ATAC Median signal | [data/source_tables/Fig1_KC/c2_Median_signal.csv](data/source_tables/Fig1_KC/c2_Median_signal.csv) |
| `c2_Saved_zscore_matrix.csv` | Figure 1: ATAC Saved zscore matrix | [data/source_tables/Fig1_KC/c2_Saved_zscore_matrix.csv](data/source_tables/Fig1_KC/c2_Saved_zscore_matrix.csv) |
| `c3_Module_summaries.csv` | Figure 1: ATAC Module summaries | [data/source_tables/Fig1_KC/c3_Module_summaries.csv](data/source_tables/Fig1_KC/c3_Module_summaries.csv) |
| `d1_Module_counts.csv` | Figure 1: H3K27ac Module counts | [data/source_tables/Fig1_KC/d1_Module_counts.csv](data/source_tables/Fig1_KC/d1_Module_counts.csv) |
| `d2_Median_signal.csv` | Figure 1: H3K27ac Median signal | [data/source_tables/Fig1_KC/d2_Median_signal.csv](data/source_tables/Fig1_KC/d2_Median_signal.csv) |
| `d3_Module_summaries.csv` | Figure 1: H3K27ac Module summaries | [data/source_tables/Fig1_KC/d3_Module_summaries.csv](data/source_tables/Fig1_KC/d3_Module_summaries.csv) |
| `e1_Module_counts.csv` | Figure 1: H3K27me3 Module counts | [data/source_tables/Fig1_KC/e1_Module_counts.csv](data/source_tables/Fig1_KC/e1_Module_counts.csv) |
| `e2_Median_signal.csv` | Figure 1: H3K27me3 Median signal | [data/source_tables/Fig1_KC/e2_Median_signal.csv](data/source_tables/Fig1_KC/e2_Median_signal.csv) |
| `e3_Module_summaries.csv` | Figure 1: H3K27me3 Module summaries | [data/source_tables/Fig1_KC/e3_Module_summaries.csv](data/source_tables/Fig1_KC/e3_Module_summaries.csv) |
| `2a_GO_plot_matrix.csv` | Figure 2: GO plot matrix | [data/source_tables/Fig2_KC/2a_GO_plot_matrix.csv](data/source_tables/Fig2_KC/2a_GO_plot_matrix.csv) |
| `2a_GO_statistics.csv` | Figure 2: GO statistics | [data/source_tables/Fig2_KC/2a_GO_statistics.csv](data/source_tables/Fig2_KC/2a_GO_statistics.csv) |
| `2b_GO_plot_matrix.csv` | Figure 2: GO plot matrix | [data/source_tables/Fig2_KC/2b_GO_plot_matrix.csv](data/source_tables/Fig2_KC/2b_GO_plot_matrix.csv) |
| `2b_GO_statistics.csv` | Figure 2: GO statistics | [data/source_tables/Fig2_KC/2b_GO_statistics.csv](data/source_tables/Fig2_KC/2b_GO_statistics.csv) |
| `2c_GO_plot_matrix.csv` | Figure 2: GO plot matrix | [data/source_tables/Fig2_KC/2c_GO_plot_matrix.csv](data/source_tables/Fig2_KC/2c_GO_plot_matrix.csv) |
| `2c_GO_statistics.csv` | Figure 2: GO statistics | [data/source_tables/Fig2_KC/2c_GO_statistics.csv](data/source_tables/Fig2_KC/2c_GO_statistics.csv) |
| `2d_VISTA_counts.csv` | Figure 2: VISTA counts | [data/source_tables/Fig2_KC/2d_VISTA_counts.csv](data/source_tables/Fig2_KC/2d_VISTA_counts.csv) |
| `2d_VISTA_tests.csv` | Figure 2: VISTA tests | [data/source_tables/Fig2_KC/2d_VISTA_tests.csv](data/source_tables/Fig2_KC/2d_VISTA_tests.csv) |
| `2dprime_VISTA_counts.csv` | Figure 2: VISTA counts | [data/source_tables/Fig2_KC/2dprime_VISTA_counts.csv](data/source_tables/Fig2_KC/2dprime_VISTA_counts.csv) |
| `2dprime_VISTA_tests.csv` | Figure 2: VISTA tests | [data/source_tables/Fig2_KC/2dprime_VISTA_tests.csv](data/source_tables/Fig2_KC/2dprime_VISTA_tests.csv) |
| `2g_Consensus_membership.csv` | Figure 2: Consensus membership | [data/source_tables/Fig2_KC/2g_Consensus_membership.csv](data/source_tables/Fig2_KC/2g_Consensus_membership.csv) |
| `2g_Curated_plot_scores.csv` | Figure 2: Curated plot scores | [data/source_tables/Fig2_KC/2g_Curated_plot_scores.csv](data/source_tables/Fig2_KC/2g_Curated_plot_scores.csv) |
| `2g_HOMER_enrichment.csv` | Figure 2: HOMER enrichment | [data/source_tables/Fig2_KC/2g_HOMER_enrichment.csv](data/source_tables/Fig2_KC/2g_HOMER_enrichment.csv) |
| `2g_Module_key.csv` | Figure 2: Module key | [data/source_tables/Fig2_KC/2g_Module_key.csv](data/source_tables/Fig2_KC/2g_Module_key.csv) |
| `2g_Selected_motif_annotations.csv` | Figure 2: Selected motif annotations | [data/source_tables/Fig2_KC/2g_Selected_motif_annotations.csv](data/source_tables/Fig2_KC/2g_Selected_motif_annotations.csv) |
| `2g_Uncapped_group_statistics.csv` | Figure 2: Uncapped group statistics | [data/source_tables/Fig2_KC/2g_Uncapped_group_statistics.csv](data/source_tables/Fig2_KC/2g_Uncapped_group_statistics.csv) |
| `3_Analysis_parameters.csv` | Figure 3: Analysis parameters | [data/source_tables/Fig3_KC/3_Analysis_parameters.csv](data/source_tables/Fig3_KC/3_Analysis_parameters.csv) |
| `3_Module_key.csv` | Figure 3: Module key | [data/source_tables/Fig3_KC/3_Module_key.csv](data/source_tables/Fig3_KC/3_Module_key.csv) |
| `3a_Cell_data.csv` | Multiome cell IDs, cluster labels and UMAP coordinates | [data/source_tables/Fig3_KC/3a_Cell_data.csv](data/source_tables/Fig3_KC/3a_Cell_data.csv) |
| `3b_Dotplot_values.csv` | Module-by-cluster dot-plot scores and positive-score percentages | [data/source_tables/Fig3_KC/3b_Dotplot_values.csv](data/source_tables/Fig3_KC/3b_Dotplot_values.csv) |
| `3c_Plot_values.csv` | Per-cell module scores and UMAP coordinates | [data/source_tables/Fig3_KC/3c_Plot_values.csv](data/source_tables/Fig3_KC/3c_Plot_values.csv) |
| `3d_Gene_cluster_means.csv` | Linked-gene expression summarized by multiome cluster | [data/source_tables/Fig3_KC/3d_Gene_cluster_means.csv](data/source_tables/Fig3_KC/3d_Gene_cluster_means.csv) |
| `3d_Module_gene_sets.csv` | Linked genes assigned to each module | [data/source_tables/Fig3_KC/3d_Module_gene_sets.csv](data/source_tables/Fig3_KC/3d_Module_gene_sets.csv) |
| `3d_Scaling_parameters.csv` | Recorded linked-gene score centering/scaling parameters | [data/source_tables/Fig3_KC/3d_Scaling_parameters.csv](data/source_tables/Fig3_KC/3d_Scaling_parameters.csv) |
| `3d_Selected_links.csv` | Selected peak-to-gene links used for module gene sets | [data/source_tables/Fig3_KC/3d_Selected_links.csv](data/source_tables/Fig3_KC/3d_Selected_links.csv) |
| `3ef_Cell_accessibility.csv` | Per-cell accessibility scores for selected module groups | [data/source_tables/Fig3_KC/3ef_Cell_accessibility.csv](data/source_tables/Fig3_KC/3ef_Cell_accessibility.csv) |
| `4a_Interaction_counts.csv` | Figure 4: Interaction counts | [data/source_tables/Fig4_KC/4a_Interaction_counts.csv](data/source_tables/Fig4_KC/4a_Interaction_counts.csv) |
| `4b_A1.csv` | Saved developmental trajectory values for ATAC module A1 | [data/source_tables/Fig4_KC/4b_A1.csv](data/source_tables/Fig4_KC/4b_A1.csv) |
| `4c_Ac7.csv` | Saved developmental trajectory values for H3K27ac module Ac7 | [data/source_tables/Fig4_KC/4c_Ac7.csv](data/source_tables/Fig4_KC/4c_Ac7.csv) |
| `4d_TM3.csv` | Saved developmental trajectory values for H3K27me3 module TM3 | [data/source_tables/Fig4_KC/4d_TM3.csv](data/source_tables/Fig4_KC/4d_TM3.csv) |
| `4f_Program_counts.csv` | Figure 4: Program counts | [data/source_tables/Fig4_KC/4f_Program_counts.csv](data/source_tables/Fig4_KC/4f_Program_counts.csv) |
| `4f_Quantized_trajectories.csv` | Figure 4: Quantized trajectories | [data/source_tables/Fig4_KC/4f_Quantized_trajectories.csv](data/source_tables/Fig4_KC/4f_Quantized_trajectories.csv) |
| `6bo_Biological_sample_sizes.csv` | Figure 6: Biological sample sizes | [data/source_tables/Fig6_KC/6bo_Biological_sample_sizes.csv](data/source_tables/Fig6_KC/6bo_Biological_sample_sizes.csv) |
| `6bprime_RPKM.csv` | Sample-level Neurod2/Ppp1r1b RPKM and biological metadata | [data/source_tables/Fig6_KC/6bprime_RPKM.csv](data/source_tables/Fig6_KC/6bprime_RPKM.csv) |
| `6d_GO_plotted.csv` | Figure 6: GO plotted | [data/source_tables/Fig6_KC/6d_GO_plotted.csv](data/source_tables/Fig6_KC/6d_GO_plotted.csv) |
| `6e_Module_enrichment.csv` | Figure 6: Module enrichment | [data/source_tables/Fig6_KC/6e_Module_enrichment.csv](data/source_tables/Fig6_KC/6e_Module_enrichment.csv) |
| `6e_Selected_membership.csv` | Figure 6: Selected membership | [data/source_tables/Fig6_KC/6e_Selected_membership.csv](data/source_tables/Fig6_KC/6e_Selected_membership.csv) |
| `6f_Counting_units.csv` | Definitions of counting units for the recovered gene–peak analysis | [data/source_tables/Fig6_KC/6f_Counting_units.csv](data/source_tables/Fig6_KC/6f_Counting_units.csv) |
| `6f_Gene_peak_counts.csv` | Recovered historical gene–peak overlap counts; final artwork totals differ | [data/source_tables/Fig6_KC/6f_Gene_peak_counts.csv](data/source_tables/Fig6_KC/6f_Gene_peak_counts.csv) |
| `6f_Gene_peak_membership.csv` | Membership underlying the recovered historical gene–peak counts | [data/source_tables/Fig6_KC/6f_Gene_peak_membership.csv](data/source_tables/Fig6_KC/6f_Gene_peak_membership.csv) |
| `6g_nonVZ_Null_overlaps.csv` | Figure 6: nonVZ Null overlaps | [data/source_tables/Fig6_KC/6g_nonVZ_Null_overlaps.csv](data/source_tables/Fig6_KC/6g_nonVZ_Null_overlaps.csv) |
| `6g_nonVZ_Statistics.csv` | Figure 6: nonVZ Statistics | [data/source_tables/Fig6_KC/6g_nonVZ_Statistics.csv](data/source_tables/Fig6_KC/6g_nonVZ_Statistics.csv) |
| `6g_VZ_Null_overlaps.csv` | Figure 6: VZ Null overlaps | [data/source_tables/Fig6_KC/6g_VZ_Null_overlaps.csv](data/source_tables/Fig6_KC/6g_VZ_Null_overlaps.csv) |
| `6g_VZ_Statistics.csv` | Figure 6: VZ Statistics | [data/source_tables/Fig6_KC/6g_VZ_Statistics.csv](data/source_tables/Fig6_KC/6g_VZ_Statistics.csv) |
| `6h_case_Variant_overlaps.csv` | Case variant–region overlap records from the saved enrichment analysis | [data/source_tables/Fig6_KC/6h_case_Variant_overlaps.csv](data/source_tables/Fig6_KC/6h_case_Variant_overlaps.csv) |
| `6h_control_Variant_overlaps.csv` | Control variant–region overlap records from the saved enrichment analysis | [data/source_tables/Fig6_KC/6h_control_Variant_overlaps.csv](data/source_tables/Fig6_KC/6h_control_Variant_overlaps.csv) |
| `6h_ME3_regions.csv` | Figure 6: ME3 regions | [data/source_tables/Fig6_KC/6h_ME3_regions.csv](data/source_tables/Fig6_KC/6h_ME3_regions.csv) |
| `6h_Per_region_counts.csv` | Figure 6: Per region counts | [data/source_tables/Fig6_KC/6h_Per_region_counts.csv](data/source_tables/Fig6_KC/6h_Per_region_counts.csv) |
| `6h_Resampling_parameters.csv` | Figure 6: Resampling parameters | [data/source_tables/Fig6_KC/6h_Resampling_parameters.csv](data/source_tables/Fig6_KC/6h_Resampling_parameters.csv) |
| `6h_Variant_null_replicates.csv` | Figure 6: Variant null replicates | [data/source_tables/Fig6_KC/6h_Variant_null_replicates.csv](data/source_tables/Fig6_KC/6h_Variant_null_replicates.csv) |
| `6h_Variant_plot_summary.csv` | Figure 6: Variant plot summary | [data/source_tables/Fig6_KC/6h_Variant_plot_summary.csv](data/source_tables/Fig6_KC/6h_Variant_plot_summary.csv) |
| `6h_Variant_statistics.csv` | Figure 6: Variant statistics | [data/source_tables/Fig6_KC/6h_Variant_statistics.csv](data/source_tables/Fig6_KC/6h_Variant_statistics.csv) |
| `6o_RPKM.csv` | Sample-level selected neuronal-gene RPKM and biological metadata | [data/source_tables/Fig6_KC/6o_RPKM.csv](data/source_tables/Fig6_KC/6o_RPKM.csv) |
| `Module_renaming_dictionary.csv` | Mapping from original module labels to manuscript names | [data/upstream/Module_renaming_dictionary.csv](data/upstream/Module_renaming_dictionary.csv) |
| `AC.linked_peaks.csv` | Saved H3K27ac peak-to-gene links and link statistics | [data/upstream/multiome/AC.linked_peaks.csv](data/upstream/multiome/AC.linked_peaks.csv) |
| `atac_modules.csv` | Derived interval-to-ATAC-module assignments | [data/upstream/multiome/atac_modules.csv](data/upstream/multiome/atac_modules.csv) |
| `ME3.linked_peaks.csv` | Saved H3K27me3 peak-to-gene links and link statistics | [data/upstream/multiome/ME3.linked_peaks.csv](data/upstream/multiome/ME3.linked_peaks.csv) |
| `RNA_pseudobulk_counts.csv` | RNA counts aggregated across cells by multiome cluster (derived pseudobulk) | [data/upstream/multiome/RNA_pseudobulk_counts.csv](data/upstream/multiome/RNA_pseudobulk_counts.csv) |
| `network_edges.csv` | Derived module-network connections | [data/upstream/network_edges.csv](data/upstream/network_edges.csv) |
| `network_nodes.csv` | Derived module-network node annotations | [data/upstream/network_nodes.csv](data/upstream/network_nodes.csv) |
| `nonVZ_wt_vs_mut_DE_Litter_Sex_cor.csv` | Original nonVZ bulk RNA differential-expression results, adjusted for litter and sex | [data/upstream/RNA/nonVZ_wt_vs_mut_DE_Litter_Sex_cor.csv](data/upstream/RNA/nonVZ_wt_vs_mut_DE_Litter_Sex_cor.csv) |
| `VZ_wt_vs_mut_DE_Litter_Sex_cor.csv` | Original VZ bulk RNA differential-expression results, adjusted for litter and sex | [data/upstream/RNA/VZ_wt_vs_mut_DE_Litter_Sex_cor.csv](data/upstream/RNA/VZ_wt_vs_mut_DE_Litter_Sex_cor.csv) |

### Saved R objects

These checkpoints are **not yet deposited on GitHub**. Paths below show where the notebooks expect them.

| File name | Contents | Repository path / availability |
|---|---|---|
| `net_WGCNA_ATAC_sp6.Rds` | Fitted ATAC WGCNA network, module colors and eigengenes | Not deposited: `data/upstream/WGCNA/WGCNA_objects/net_WGCNA_ATAC_sp6.Rds` |
| `net_WGCNA_H3K27ac_merged_by_1kb_sp5.Rds` | Fitted H3K27ac WGCNA network | Not deposited: `data/upstream/WGCNA/WGCNA_objects/net_WGCNA_H3K27ac_merged_by_1kb_sp5.Rds` |
| `net_WGCNA_H3K27me3_merged_by_1kb_sp9.Rds` | Fitted H3K27me3 WGCNA network | Not deposited: `data/upstream/WGCNA/WGCNA_objects/net_WGCNA_H3K27me3_merged_by_1kb_sp9.Rds` |
| `exact_V2_cache.rds` | Stored normalized multiome assays, cell metadata and links | Not deposited: `objects/Fig3/exact_V2_cache.rds` |
| `ATAC_scoring.rds` | Saved ATAC module features, control features and scores | Not deposited: `objects/Fig3/ATAC_scoring.rds` |
| `RNA_scoring.rds` | Saved linked-gene features, control features and scores | Not deposited: `objects/Fig3/RNA_scoring.rds` |
| `umap_reclustering_bundles.RData` | Saved integration embeddings and cluster assignments | Not deposited: `objects/integration/umap_reclustering_bundles.RData` |
| `extended5_linked_scores.rds` | AC/ME3 linked-gene scores and controls for Extended Data 5e/f | Not deposited: `objects/multiome/extended5_linked_scores.rds` |
| `marker_heatmap_inputs.rds` | Saved expression/accessibility layers for Extended Data 5g/i/k | Not deposited: `objects/multiome/marker_heatmap_inputs.rds` |
| `5l_positive_DA_composition.rds` | Derived counts, selected intervals, filtering and pooling for panel 5l | Not deposited: `objects/multiome/5l_positive_DA_composition.rds` |

[Complete saved-object inventory](SAVED_OBJECTS.md).

Shared R functions: [R/](R/).
