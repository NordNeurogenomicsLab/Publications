# Saved R objects

This lists **all 57 local RDS/RData/RDA files**, not only those required by routine rendering. These files are available in the [public Box release](https://ucdavis.box.com/s/jrmg33n2dcljxadenvx47ji8c6azg02s); retain the relative paths below when copying them into the code directory. See [INPUTS.md](INPUTS.md). `readRDS()` loads an RDS object; `load(..., envir=...)` loads named objects from an RData/RDA workspace. Internal object schemas are inspected in the notebooks that consume them.

Size is decimal MB. “No declared reader” means retained/generated material, not a required input inferred from its filename.

[Release file paths and checksums](config/BOX_RELEASE_FILES.csv).

## Main checkpoints

| Object | Contents / role |
|---|---|
| [`data/upstream/WGCNA/WGCNA_objects/net_WGCNA_ATAC_sp6.Rds`](https://ucdavis.app.box.com/s/jrmg33n2dcljxadenvx47ji8c6azg02s/file/2575319031742) | Fitted ATAC network; original module colors and eigengenes |
| [`data/upstream/WGCNA/WGCNA_objects/net_WGCNA_H3K27ac_merged_by_1kb_sp5.Rds`](https://ucdavis.app.box.com/s/jrmg33n2dcljxadenvx47ji8c6azg02s/file/2575320489170) | Fitted AC network |
| [`data/upstream/WGCNA/WGCNA_objects/net_WGCNA_H3K27me3_merged_by_1kb_sp9.Rds`](https://ucdavis.app.box.com/s/jrmg33n2dcljxadenvx47ji8c6azg02s/file/2575312089542) | Fitted ME3 network |
| [`objects/Fig3/exact_V2_cache.rds`](https://ucdavis.app.box.com/s/jrmg33n2dcljxadenvx47ji8c6azg02s/file/2575315640767) | Compact stored multiome assays, cell metadata and links; used for optional module-score reconstruction |
| [`objects/Fig3/ATAC_scoring.rds`](https://ucdavis.app.box.com/s/jrmg33n2dcljxadenvx47ji8c6azg02s/file/2575315470293), [`RNA_scoring.rds`](https://ucdavis.app.box.com/s/jrmg33n2dcljxadenvx47ji8c6azg02s/file/2575324450997) | Gene/interval feature sets, control sets and saved module scores; inspected in pipeline notebook 03 |
| [`objects/integration/umap_reclustering_bundles.RData`](https://ucdavis.app.box.com/s/jrmg33n2dcljxadenvx47ji8c6azg02s/file/2575312064873) | Saved linked-interval embeddings and cluster assignments; loaded by `R/program_summary_saved.R` |
| [`objects/multiome/extended5_linked_scores.rds`](https://ucdavis.app.box.com/s/jrmg33n2dcljxadenvx47ji8c6azg02s/file/2575323333107) | AC/ME3 linked-gene scores, features, controls and cell metadata for Extended Data 5e/f |
| [`objects/multiome/marker_heatmap_inputs.rds`](https://ucdavis.app.box.com/s/jrmg33n2dcljxadenvx47ji8c6azg02s/file/2575321480370) | Stored SCT scaled expression and interval/peak layers for Extended Data 5g/i/k |
| [`objects/multiome/5l_positive_DA_composition.rds`](https://ucdavis.app.box.com/s/jrmg33n2dcljxadenvx47ji8c6azg02s/file/2575314445837) | Derived panel 5l counts, selected intervals, filter, pooling and palette; created by notebook 03 |
| `objects/source_tables/*.rds` | Figure-specific table caches created by notebook 00; CSVs remain the direct plotting inputs |
| `objects/author_checkpoints/*.rds` | Retained historical session objects; not required merely because they are archived |

## `data/upstream/WGCNA/WGCNA_input_data`

| File | MB | Reader / role |
|---|---:|---|
| [`data_h3k27ac_merged_by_1kb.RData`](https://ucdavis.app.box.com/s/jrmg33n2dcljxadenvx47ji8c6azg02s/file/2575316869726) | 2.05 | 00_inputs_and_provenance.Rmd (saved-results) |
| [`data_h3k27me3_merged_by_1kb.RData`](https://ucdavis.app.box.com/s/jrmg33n2dcljxadenvx47ji8c6azg02s/file/2575320438770) | 1.96 | 00_inputs_and_provenance.Rmd (saved-results) |

## `data/upstream/WGCNA/WGCNA_objects`

| File | MB | Reader / role |
|---|---:|---|
| [`net_WGCNA_ATAC_sp6.Rds`](https://ucdavis.app.box.com/s/jrmg33n2dcljxadenvx47ji8c6azg02s/file/2575319031742) | 2.50 | 00_inputs_and_provenance.Rmd; 01_WGCNA_checkpoints.Rmd; pipeline/01_WGCNA_network_generation.Rmd (saved-results) |
| [`net_WGCNA_H3K27ac_merged_by_1kb_sp5.Rds`](https://ucdavis.app.box.com/s/jrmg33n2dcljxadenvx47ji8c6azg02s/file/2575320489170) | 2.43 | 00_inputs_and_provenance.Rmd; 01_WGCNA_checkpoints.Rmd; pipeline/01_WGCNA_network_generation.Rmd (saved-results) |
| [`net_WGCNA_H3K27me3_merged_by_1kb_sp9.Rds`](https://ucdavis.app.box.com/s/jrmg33n2dcljxadenvx47ji8c6azg02s/file/2575312089542) | 2.63 | 00_inputs_and_provenance.Rmd; 01_WGCNA_checkpoints.Rmd; pipeline/01_WGCNA_network_generation.Rmd (saved-results) |

## `objects/Fig3`

| File | MB | Reader / role |
|---|---:|---|
| [`ATAC_scoring.rds`](https://ucdavis.app.box.com/s/jrmg33n2dcljxadenvx47ji8c6azg02s/file/2575315470293) | 29.41 | 00_inputs_and_provenance.Rmd; pipeline/03_multiome_module_scores.Rmd (saved-results) |
| [`exact_V2_cache.rds`](https://ucdavis.app.box.com/s/jrmg33n2dcljxadenvx47ji8c6azg02s/file/2575315640767) | 2313.96 | 00_inputs_and_provenance.Rmd; 03_Figure_3_multiome.Rmd; 11_Extended_multiome_heatmaps.Rmd; pipeline/03_multiome_module_scores.Rmd (optional-rerun; saved-results) |
| [`heatmap_plot_object.rds`](https://ucdavis.app.box.com/s/jrmg33n2dcljxadenvx47ji8c6azg02s/file/2575322560020) | 2.01 | No declared reader; retained/generated checkpoint (retained or generated artifact; not declared as a direct input) |
| [`reconstruction_results.rds`](https://ucdavis.app.box.com/s/jrmg33n2dcljxadenvx47ji8c6azg02s/file/2575314870148) | 48.96 | 00_inputs_and_provenance.Rmd (saved-results) |
| [`RNA_scoring.rds`](https://ucdavis.app.box.com/s/jrmg33n2dcljxadenvx47ji8c6azg02s/file/2575324450997) | 9.18 | 00_inputs_and_provenance.Rmd; pipeline/03_multiome_module_scores.Rmd (saved-results) |

## `objects/RNA`

| File | MB | Reader / role |
|---|---:|---|
| [`verified_plotting_inputs.rds`](https://ucdavis.app.box.com/s/jrmg33n2dcljxadenvx47ji8c6azg02s/file/2575314258978) | 0.01 | 00_inputs_and_provenance.Rmd (saved-results) |

## `objects/WGCNA`

| File | MB | Reader / role |
|---|---:|---|
| [`final_networks.rds`](https://ucdavis.app.box.com/s/jrmg33n2dcljxadenvx47ji8c6azg02s/file/2575313975192) | 7.56 | No declared reader; retained/generated checkpoint (retained or generated artifact; not declared as a direct input) |

## `objects/author_checkpoints`

| File | MB | Reader / role |
|---|---:|---|
| [`Integration2_12_15_2025.rds`](https://ucdavis.app.box.com/s/jrmg33n2dcljxadenvx47ji8c6azg02s/file/2575316615949) | 14.25 | No declared reader; retained/generated checkpoint (retained or generated artifact; not declared as a direct input) |
| [`Session_dump_annotations_Fig1_and_2..rds`](https://ucdavis.app.box.com/s/jrmg33n2dcljxadenvx47ji8c6azg02s/file/2575318960543) | 0.09 | No declared reader; retained/generated checkpoint (retained or generated artifact; not declared as a direct input) |
| [`Session_dump_FIg6_12_06_2025.rds`](https://ucdavis.app.box.com/s/jrmg33n2dcljxadenvx47ji8c6azg02s/file/2575322658329) | 64.12 | No declared reader; retained/generated checkpoint (retained or generated artifact; not declared as a direct input) |
| [`Session_dump_Integration_10_30_2025v2.rds`](https://ucdavis.app.box.com/s/jrmg33n2dcljxadenvx47ji8c6azg02s/file/2575317428954) | 50.99 | No declared reader; retained/generated checkpoint (retained or generated artifact; not declared as a direct input) |

## `objects/histology`

| File | MB | Reader / role |
|---|---:|---|
| [`5e_TLE4_E14.rds`](https://ucdavis.app.box.com/s/jrmg33n2dcljxadenvx47ji8c6azg02s/file/2575321295310) | 0.00 | No declared reader; retained/generated checkpoint (retained or generated artifact; not declared as a direct input) |
| [`5f_Tbr1_E14.rds`](https://ucdavis.app.box.com/s/jrmg33n2dcljxadenvx47ji8c6azg02s/file/2575322797723) | 0.00 | No declared reader; retained/generated checkpoint (retained or generated artifact; not declared as a direct input) |
| [`5g_Tbr1_P10.rds`](https://ucdavis.app.box.com/s/jrmg33n2dcljxadenvx47ji8c6azg02s/file/2575320294264) | 0.00 | No declared reader; retained/generated checkpoint (retained or generated artifact; not declared as a direct input) |
| [`5h_Bcl11b_P10.rds`](https://ucdavis.app.box.com/s/jrmg33n2dcljxadenvx47ji8c6azg02s/file/2575320160185) | 0.00 | No declared reader; retained/generated checkpoint (retained or generated artifact; not declared as a direct input) |
| [`5i_Fezf2_P10.rds`](https://ucdavis.app.box.com/s/jrmg33n2dcljxadenvx47ji8c6azg02s/file/2575310386842) | 0.00 | No declared reader; retained/generated checkpoint (retained or generated artifact; not declared as a direct input) |
| [`7j_CP_micrometres.rds`](https://ucdavis.app.box.com/s/jrmg33n2dcljxadenvx47ji8c6azg02s/file/2575318568857) | 0.00 | No declared reader; retained/generated checkpoint (retained or generated artifact; not declared as a direct input) |
| [`7j_TLE4_micrometres.rds`](https://ucdavis.app.box.com/s/jrmg33n2dcljxadenvx47ji8c6azg02s/file/2575323652675) | 0.00 | No declared reader; retained/generated checkpoint (retained or generated artifact; not declared as a direct input) |
| [`7j_TLE4_thickness.rds`](https://ucdavis.app.box.com/s/jrmg33n2dcljxadenvx47ji8c6azg02s/file/2575310919920) | 0.00 | No declared reader; retained/generated checkpoint (retained or generated artifact; not declared as a direct input) |
| [`7k_TBR1_layer_thickness.rds`](https://ucdavis.app.box.com/s/jrmg33n2dcljxadenvx47ji8c6azg02s/file/2575316197927) | 0.00 | No declared reader; retained/generated checkpoint (retained or generated artifact; not declared as a direct input) |
| [`7kl_Tbr1_416Del.rds`](https://ucdavis.app.box.com/s/jrmg33n2dcljxadenvx47ji8c6azg02s/file/2575320412035) | 0.00 | No declared reader; retained/generated checkpoint (retained or generated artifact; not declared as a direct input) |
| [`9f_Tbr1_Nkx2-1.rds`](https://ucdavis.app.box.com/s/jrmg33n2dcljxadenvx47ji8c6azg02s/file/2575307828527) | 0.00 | No declared reader; retained/generated checkpoint (retained or generated artifact; not declared as a direct input) |

## `objects/integration`

| File | MB | Reader / role |
|---|---:|---|
| [`umap_reclustering_bundles.RData`](https://ucdavis.app.box.com/s/jrmg33n2dcljxadenvx47ji8c6azg02s/file/2575312064873) | 571.00 | 04_Figure_4_integration.Rmd (saved-results) |

## `objects/multiome`

| File | MB | Reader / role |
|---|---:|---|
| [`5d_ATAC_cluster_means_matched.rds`](https://ucdavis.app.box.com/s/jrmg33n2dcljxadenvx47ji8c6azg02s/file/2575318377546) | 3.56 | No declared reader; retained/generated checkpoint (retained or generated artifact; not declared as a direct input) |
| [`5e_AC_linked_gene_scores.rds`](https://ucdavis.app.box.com/s/jrmg33n2dcljxadenvx47ji8c6azg02s/file/2575321119790) | 1.73 | No declared reader; retained/generated checkpoint (retained or generated artifact; not declared as a direct input) |
| [`5f_ME3_linked_gene_scores.rds`](https://ucdavis.app.box.com/s/jrmg33n2dcljxadenvx47ji8c6azg02s/file/2575320038086) | 1.73 | No declared reader; retained/generated checkpoint (retained or generated artifact; not declared as a direct input) |
| [`5g_marker_heatmap.rds`](https://ucdavis.app.box.com/s/jrmg33n2dcljxadenvx47ji8c6azg02s/file/2575312026041) | 14.68 | No declared reader; retained/generated checkpoint (retained or generated artifact; not declared as a direct input) |
| [`5h_PCA.rds`](https://ucdavis.app.box.com/s/jrmg33n2dcljxadenvx47ji8c6azg02s/file/2575317985817) | 0.07 | No declared reader; retained/generated checkpoint (retained or generated artifact; not declared as a direct input) |
| [`5i_marker_heatmap.rds`](https://ucdavis.app.box.com/s/jrmg33n2dcljxadenvx47ji8c6azg02s/file/2575322557907) | 2.50 | No declared reader; retained/generated checkpoint (retained or generated artifact; not declared as a direct input) |
| [`5j_PCA.rds`](https://ucdavis.app.box.com/s/jrmg33n2dcljxadenvx47ji8c6azg02s/file/2575315557339) | 0.07 | No declared reader; retained/generated checkpoint (retained or generated artifact; not declared as a direct input) |
| [`5k_marker_heatmap.rds`](https://ucdavis.app.box.com/s/jrmg33n2dcljxadenvx47ji8c6azg02s/file/2575321456727) | 2.52 | No declared reader; retained/generated checkpoint (retained or generated artifact; not declared as a direct input) |
| [`5l_positive_DA_composition.rds`](https://ucdavis.app.box.com/s/jrmg33n2dcljxadenvx47ji8c6azg02s/file/2575314445837) | 1.11 | No declared reader; retained/generated checkpoint (retained or generated artifact; not declared as a direct input) |
| [`6a_TMK3_linked_ATAC_matched.rds`](https://ucdavis.app.box.com/s/jrmg33n2dcljxadenvx47ji8c6azg02s/file/2575321487925) | 0.49 | No declared reader; retained/generated checkpoint (retained or generated artifact; not declared as a direct input) |
| [`extended5_linked_scores.rds`](https://ucdavis.app.box.com/s/jrmg33n2dcljxadenvx47ji8c6azg02s/file/2575323333107) | 10.11 | 11_Extended_multiome_heatmaps.Rmd (saved-results) |
| [`marker_heatmap_inputs.rds`](https://ucdavis.app.box.com/s/jrmg33n2dcljxadenvx47ji8c6azg02s/file/2575321480370) | 1173.28 | 11_Extended_multiome_heatmaps.Rmd (saved-results) |

## `objects/source_tables`

| File | MB | Reader / role |
|---|---:|---|
| [`Extended_Fig1_KC.rds`](https://ucdavis.app.box.com/s/jrmg33n2dcljxadenvx47ji8c6azg02s/file/2575319934250) | 8.46 | No declared reader; retained/generated checkpoint (retained or generated artifact; not declared as a direct input) |
| [`Extended_Fig2_KC.rds`](https://ucdavis.app.box.com/s/jrmg33n2dcljxadenvx47ji8c6azg02s/file/2575319790940) | 0.13 | No declared reader; retained/generated checkpoint (retained or generated artifact; not declared as a direct input) |
| [`Extended_Fig4_KC.rds`](https://ucdavis.app.box.com/s/jrmg33n2dcljxadenvx47ji8c6azg02s/file/2575321532800) | 0.00 | No declared reader; retained/generated checkpoint (retained or generated artifact; not declared as a direct input) |
| [`Extended_Fig5_KC.rds`](https://ucdavis.app.box.com/s/jrmg33n2dcljxadenvx47ji8c6azg02s/file/2575320815740) | 12.47 | No declared reader; retained/generated checkpoint (retained or generated artifact; not declared as a direct input) |
| [`Extended_Fig6_KC.rds`](https://ucdavis.app.box.com/s/jrmg33n2dcljxadenvx47ji8c6azg02s/file/2575315542330) | 0.86 | No declared reader; retained/generated checkpoint (retained or generated artifact; not declared as a direct input) |
| [`Extended_Fig7_KC.rds`](https://ucdavis.app.box.com/s/jrmg33n2dcljxadenvx47ji8c6azg02s/file/2575314851244) | 0.38 | No declared reader; retained/generated checkpoint (retained or generated artifact; not declared as a direct input) |
| [`Extended_Fig8_KC.rds`](https://ucdavis.app.box.com/s/jrmg33n2dcljxadenvx47ji8c6azg02s/file/2575314195712) | 11.48 | No declared reader; retained/generated checkpoint (retained or generated artifact; not declared as a direct input) |
| [`Extended_Fig9_KC.rds`](https://ucdavis.app.box.com/s/jrmg33n2dcljxadenvx47ji8c6azg02s/file/2575313989986) | 0.00 | No declared reader; retained/generated checkpoint (retained or generated artifact; not declared as a direct input) |
| [`Fig1_KC.rds`](https://ucdavis.app.box.com/s/jrmg33n2dcljxadenvx47ji8c6azg02s/file/2575310718143) | 5.44 | No declared reader; retained/generated checkpoint (retained or generated artifact; not declared as a direct input) |
| [`Fig2_KC.rds`](https://ucdavis.app.box.com/s/jrmg33n2dcljxadenvx47ji8c6azg02s/file/2575321377525) | 0.08 | No declared reader; retained/generated checkpoint (retained or generated artifact; not declared as a direct input) |
| [`Fig3_KC.rds`](https://ucdavis.app.box.com/s/jrmg33n2dcljxadenvx47ji8c6azg02s/file/2575316224257) | 8.83 | No declared reader; retained/generated checkpoint (retained or generated artifact; not declared as a direct input) |
| [`Fig4_KC.rds`](https://ucdavis.app.box.com/s/jrmg33n2dcljxadenvx47ji8c6azg02s/file/2575311167038) | 0.02 | No declared reader; retained/generated checkpoint (retained or generated artifact; not declared as a direct input) |
| [`Fig5_KC.rds`](https://ucdavis.app.box.com/s/jrmg33n2dcljxadenvx47ji8c6azg02s/file/2575321605043) | 1.60 | No declared reader; retained/generated checkpoint (retained or generated artifact; not declared as a direct input) |
| [`Fig6_KC.rds`](https://ucdavis.app.box.com/s/jrmg33n2dcljxadenvx47ji8c6azg02s/file/2575310694143) | 4.17 | No declared reader; retained/generated checkpoint (retained or generated artifact; not declared as a direct input) |
| [`Fig7_KC.rds`](https://ucdavis.app.box.com/s/jrmg33n2dcljxadenvx47ji8c6azg02s/file/2575319913539) | 0.00 | No declared reader; retained/generated checkpoint (retained or generated artifact; not declared as a direct input) |

## `objects/variants`

| File | MB | Reader / role |
|---|---:|---|
| [`Fu_and_Kdm6b_DE_enr.RDS`](https://ucdavis.app.box.com/s/jrmg33n2dcljxadenvx47ji8c6azg02s/file/2575314116512) | 8.19 | 00_inputs_and_provenance.Rmd (saved-results) |
| [`saved_resampling.rds`](https://ucdavis.app.box.com/s/jrmg33n2dcljxadenvx47ji8c6azg02s/file/2575318153440) | 0.01 | No declared reader; retained/generated checkpoint (retained or generated artifact; not declared as a direct input) |

The [machine-readable object inventory](config/saved_objects.csv) includes full paths, original source identities where recorded, and copy hashes where available.
