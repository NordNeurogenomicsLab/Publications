## Manuscript

**Profiling dynamic epigenetic programs in the developing mouse cortex reveals the role of ASD-associated Kdm6b in deep-layer neuron differentiation.**

Athéna R. Ypsilanti<sup>†,*</sup>, Eirene Markenscoff-Papadimitriou<sup>†,*</sup>, Karol Cichewicz<sup>†</sup>, Sean Whalen, Nicolas Seban, Ethan Fenton, Kenneth Lim, Shan Dong, Jada Baldoza, Yurong Wang, Matt Barraza, Fadya Binyameen, Diane Dickel, Len A. Pennacchio, Axel Visel, Stephan J. Sanders, Katherine S. Pollard, Alex S. Nord<sup>*</sup>, John Rubenstein<sup>*</sup>

† Contributed equally: Athéna R. Ypsilanti, Eirene Markenscoff-Papadimitriou and Karol Cichewicz.

\* Corresponding authors: Athéna R. Ypsilanti, Eirene Markenscoff-Papadimitriou, Alex S. Nord and John Rubenstein.

## Inputs and environment

[Source_tables](SOURCE_TABLES.md) · [All data files](config/BOX_RELEASE_FILES.csv) · [Input files and run order](INPUTS.md) · [Saved RDS/RData objects](SAVED_OBJECTS.md) · [R session and package versions](environment/README.md)

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

[Box data and saved objects](https://ucdavis.box.com/s/jrmg33n2dcljxadenvx47ji8c6azg02s) (436 files, 5.3 GB). Download the release and merge its `data/` and `objects/` folders into this code directory. [Paths and checksums](config/BOX_RELEASE_FILES.csv).

### Primary RNA-seq and chromatin count tables

<table width="100%">
<thead><tr><th width="40%">File name</th><th width="52%">Contents</th><th width="8%">File</th></tr></thead>
<tbody>
<tr><td><code>featureCounts_final_15samples.txt</code></td><td>Kdm6b cKO and WT bulk RNA-seq gene counts and feature lengths (15 libraries)</td><td><a href="https://ucdavis.app.box.com/s/jrmg33n2dcljxadenvx47ji8c6azg02s/file/2575320512459" title="data/upstream/RNA/featureCounts_final_15samples.txt">File</a></td></tr>
<tr><td><code>hac_overlap_peaks_merged-mm10.tsv</code></td><td>Developmental H3K27ac ChIP-seq peak-by-sample counts</td><td><a href="https://ucdavis.app.box.com/s/jrmg33n2dcljxadenvx47ji8c6azg02s/file/2575362263134" title="data/upstream/WGCNA/WGCNA_input_data/hac_overlap_peaks_merged-mm10.tsv">File</a></td></tr>
<tr><td><code>htk_overlap_peaks_merged-mm10.tsv</code></td><td>Developmental H3K27me3 ChIP-seq peak-by-sample counts</td><td><a href="https://ucdavis.app.box.com/s/jrmg33n2dcljxadenvx47ji8c6azg02s/file/2575366744761" title="data/upstream/WGCNA/WGCNA_input_data/htk_overlap_peaks_merged-mm10.tsv">File</a></td></tr>
<tr><td><code>atac_idr_peaks_merged_uniform-mm10.tsv</code></td><td>Developmental ATAC-seq peak-by-sample counts</td><td><a href="https://ucdavis.app.box.com/s/jrmg33n2dcljxadenvx47ji8c6azg02s/file/2575367425504" title="data/upstream/WGCNA/WGCNA_input_data/atac_idr_peaks_merged_uniform-mm10.tsv">File</a></td></tr>
</tbody>
</table>

Kdm6b cKO chromatin and RNA differential results: [ATAC/ChIP-seq and RNA-seq results](#kdm6b-cko-diffbind-and-rna-seq-results).

### Metadata

<table width="100%">
<thead><tr><th width="40%">File name</th><th width="52%">Contents</th><th width="8%">File</th></tr></thead>
<tbody>
<tr><td><code>sample_metadata.csv</code></td><td>Bulk RNA library IDs, genotype, litter, sex and FlashTag population</td><td><a href="https://ucdavis.app.box.com/s/jrmg33n2dcljxadenvx47ji8c6azg02s/file/2575361884550" title="data/upstream/RNA/sample_metadata.csv">File</a></td></tr>
<tr><td><code>sample_metadata.csv</code></td><td>ATAC/ChIP-seq sample metadata and BAM basenames for DiffBind</td><td><a href="https://ucdavis.app.box.com/s/jrmg33n2dcljxadenvx47ji8c6azg02s/file/2575320687685" title="data/upstream/diffbind/sample_metadata.csv">File</a></td></tr>
<tr><td><code>multiome_cell_metadata.csv</code></td><td>Multiome cell identities and annotations</td><td><a href="https://ucdavis.app.box.com/s/jrmg33n2dcljxadenvx47ji8c6azg02s/file/2575321679668" title="data/upstream/multiome/multiome_cell_metadata.csv">File</a></td></tr>
</tbody>
</table>

### WGCNA annotations and integrated data

<table width="100%">
<thead><tr><th width="40%">File name</th><th width="52%">Contents</th><th width="8%">File</th></tr></thead>
<tbody>
<tr><td><code>Module_renaming_dictionary.csv</code></td><td>Original module labels mapped to manuscript names</td><td><a href="https://ucdavis.app.box.com/s/jrmg33n2dcljxadenvx47ji8c6azg02s/file/2575359049988" title="data/upstream/Module_renaming_dictionary.csv">File</a></td></tr>
<tr><td><code>atac_modules.csv</code></td><td>ATAC interval-to-module assignments</td><td><a href="https://ucdavis.app.box.com/s/jrmg33n2dcljxadenvx47ji8c6azg02s/file/2575355675103" title="data/upstream/multiome/atac_modules.csv">File</a></td></tr>
<tr><td><code>Supp_Table_1_ATAC_master.csv</code></td><td>Final annotated ATAC intervals and modules</td><td><a href="https://ucdavis.app.box.com/s/jrmg33n2dcljxadenvx47ji8c6azg02s/file/2575313834073" title="data/upstream/final_tables/Supp_Table_1_ATAC_master.csv">File</a></td></tr>
<tr><td><code>Supp_Table_2_AC_master.csv</code></td><td>Final annotated H3K27ac intervals and modules</td><td><a href="https://ucdavis.app.box.com/s/jrmg33n2dcljxadenvx47ji8c6azg02s/file/2575321878803" title="data/upstream/final_tables/Supp_Table_2_AC_master.csv">File</a></td></tr>
<tr><td><code>Supp_Table_3_ME3_master.csv</code></td><td>Final annotated H3K27me3 intervals and modules</td><td><a href="https://ucdavis.app.box.com/s/jrmg33n2dcljxadenvx47ji8c6azg02s/file/2575315736751" title="data/upstream/final_tables/Supp_Table_3_ME3_master.csv">File</a></td></tr>
<tr><td><code>Supp_Table_4_WGCNA_interactions_and_Programs.csv</code></td><td>Integrated WGCNA interval interactions and developmental programs</td><td><a href="https://ucdavis.app.box.com/s/jrmg33n2dcljxadenvx47ji8c6azg02s/file/2575312699106" title="data/upstream/final_tables/Supp_Table_4_WGCNA_interactions_and_Programs.csv">File</a></td></tr>
<tr><td><code>all_overlaps_total_aug_with_Programs.csv</code></td><td>Linked ATAC/AC/ME3 interval catalog with program assignments</td><td><a href="https://ucdavis.app.box.com/s/jrmg33n2dcljxadenvx47ji8c6azg02s/file/2575313800586" title="data/upstream/integration/all_overlaps_total_aug_with_Programs.csv">File</a></td></tr>
<tr><td><code>AC.linked_peaks.csv</code></td><td>H3K27ac peak-to-gene links and statistics</td><td><a href="https://ucdavis.app.box.com/s/jrmg33n2dcljxadenvx47ji8c6azg02s/file/2575360382339" title="data/upstream/multiome/AC.linked_peaks.csv">File</a></td></tr>
<tr><td><code>ME3.linked_peaks.csv</code></td><td>H3K27me3 peak-to-gene links and statistics</td><td><a href="https://ucdavis.app.box.com/s/jrmg33n2dcljxadenvx47ji8c6azg02s/file/2575360900035" title="data/upstream/multiome/ME3.linked_peaks.csv">File</a></td></tr>
<tr><td><code>network_nodes.csv</code></td><td>Integrated module-network node annotations</td><td><a href="https://ucdavis.app.box.com/s/jrmg33n2dcljxadenvx47ji8c6azg02s/file/2575360585769" title="data/upstream/network_nodes.csv">File</a></td></tr>
<tr><td><code>network_edges.csv</code></td><td>Integrated module-network connections</td><td><a href="https://ucdavis.app.box.com/s/jrmg33n2dcljxadenvx47ji8c6azg02s/file/2575366835464" title="data/upstream/network_edges.csv">File</a></td></tr>
<tr><td><code>4f_Program_counts.csv</code></td><td>Counts assigned to each integrated developmental program</td><td><a href="https://ucdavis.app.box.com/s/jrmg33n2dcljxadenvx47ji8c6azg02s/file/2575361211869" title="data/source_tables/Fig4_KC/4f_Program_counts.csv">File</a></td></tr>
<tr><td><code>4f_Quantized_trajectories.csv</code></td><td>Quantized ATAC/AC/ME3 program trajectories</td><td><a href="https://ucdavis.app.box.com/s/jrmg33n2dcljxadenvx47ji8c6azg02s/file/2575360362045" title="data/source_tables/Fig4_KC/4f_Quantized_trajectories.csv">File</a></td></tr>
</tbody>
</table>

### Kdm6b cKO DiffBind and RNA-seq results

These are derived statistical results; the chromatin tables also retain per-sample signal columns. They are not raw read-count matrices.

<table width="100%">
<thead><tr><th width="40%">File name</th><th width="52%">Contents</th><th width="8%">File</th></tr></thead>
<tbody>
<tr><td><code>Ac_FTneg_Modality_H3K27Ac_<br>PeakSet_ac_expand_bed_5000_<br>Condition_FlashtagNeg__genotype_DBA_DESEQ2.csv</code></td><td>H3K27ac Kdm6b cKO versus WT DiffBind/DESeq2 results, FlashTag-negative cells</td><td><a href="https://ucdavis.app.box.com/s/jrmg33n2dcljxadenvx47ji8c6azg02s/file/2575319214240" title="data/upstream/final_tables/Ac_FTneg_Modality_H3K27Ac_PeakSet_ac_expand_bed_5000_Condition_FlashtagNeg__genotype_DBA_DESEQ2.csv">File</a></td></tr>
<tr><td><code>Ac_FTpos_Modality_H3K27Ac_<br>PeakSet_ac_expand_bed_5000_<br>Condition_FlashtagPos__genotype_DBA_DESEQ2.csv</code></td><td>H3K27ac Kdm6b cKO versus WT DiffBind/DESeq2 results, FlashTag-positive cells</td><td><a href="https://ucdavis.app.box.com/s/jrmg33n2dcljxadenvx47ji8c6azg02s/file/2575320037788" title="data/upstream/final_tables/Ac_FTpos_Modality_H3K27Ac_PeakSet_ac_expand_bed_5000_Condition_FlashtagPos__genotype_DBA_DESEQ2.csv">File</a></td></tr>
<tr><td><code>ATAC_FTneg_Modality_ATAC_<br>PeakSet_atac_expand_bed_5000_<br>Condition_FlashtagNeg__genotype_DBA_DESEQ2.csv</code></td><td>ATAC Kdm6b cKO versus WT DiffBind/DESeq2 results, FlashTag-negative cells</td><td><a href="https://ucdavis.app.box.com/s/jrmg33n2dcljxadenvx47ji8c6azg02s/file/2575311270143" title="data/upstream/final_tables/ATAC_FTneg_Modality_ATAC_PeakSet_atac_expand_bed_5000_Condition_FlashtagNeg__genotype_DBA_DESEQ2.csv">File</a></td></tr>
<tr><td><code>ATAC_FTpos_Modality_ATAC_<br>PeakSet_atac_expand_bed_5000_<br>Condition_FlashtagPos__genotype_DBA_DESEQ2.csv</code></td><td>ATAC Kdm6b cKO versus WT DiffBind/DESeq2 results, FlashTag-positive cells</td><td><a href="https://ucdavis.app.box.com/s/jrmg33n2dcljxadenvx47ji8c6azg02s/file/2575321875790" title="data/upstream/final_tables/ATAC_FTpos_Modality_ATAC_PeakSet_atac_expand_bed_5000_Condition_FlashtagPos__genotype_DBA_DESEQ2.csv">File</a></td></tr>
<tr><td><code>H3K27Me3_FTneg_Modality_H3K27Me3_<br>PeakSet_me3_expand_bed_5000_<br>Condition_FlashtagNeg__genotype_DBA_DESEQ2.csv</code></td><td>H3K27me3 Kdm6b cKO versus WT DiffBind/DESeq2 results, FlashTag-negative cells</td><td><a href="https://ucdavis.app.box.com/s/jrmg33n2dcljxadenvx47ji8c6azg02s/file/2575311414133" title="data/upstream/final_tables/H3K27Me3_FTneg_Modality_H3K27Me3_PeakSet_me3_expand_bed_5000_Condition_FlashtagNeg__genotype_DBA_DESEQ2.csv">File</a></td></tr>
<tr><td><code>H3K27Me3_FTpos_Modality_H3K27Me3_<br>PeakSet_me3_expand_bed_5000_<br>Condition_FlashtagPos__genotype_DBA_DESEQ2.csv</code></td><td>H3K27me3 Kdm6b cKO versus WT DiffBind/DESeq2 results, FlashTag-positive cells</td><td><a href="https://ucdavis.app.box.com/s/jrmg33n2dcljxadenvx47ji8c6azg02s/file/2575311219711" title="data/upstream/final_tables/H3K27Me3_FTpos_Modality_H3K27Me3_PeakSet_me3_expand_bed_5000_Condition_FlashtagPos__genotype_DBA_DESEQ2.csv">File</a></td></tr>
<tr><td><code>VZ_wt_vs_mut_DE_Litter_Sex_cor.csv</code></td><td>VZ Kdm6b cKO versus WT RNA-seq differential expression, adjusted for litter and sex</td><td><a href="https://ucdavis.app.box.com/s/jrmg33n2dcljxadenvx47ji8c6azg02s/file/2575366734159" title="data/upstream/RNA/VZ_wt_vs_mut_DE_Litter_Sex_cor.csv">File</a></td></tr>
<tr><td><code>nonVZ_wt_vs_mut_DE_Litter_Sex_cor.csv</code></td><td>nonVZ Kdm6b cKO versus WT RNA-seq differential expression, adjusted for litter and sex</td><td><a href="https://ucdavis.app.box.com/s/jrmg33n2dcljxadenvx47ji8c6azg02s/file/2575361031850" title="data/upstream/RNA/nonVZ_wt_vs_mut_DE_Litter_Sex_cor.csv">File</a></td></tr>
</tbody>
</table>

Other figure-level tables, IHC/thickness measurements and saved-object references: [Source_tables](SOURCE_TABLES.md).

Shared R functions: [R/](R/).
