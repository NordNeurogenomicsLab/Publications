# Data availability and historical execution details

The local delivery includes final plotting CSVs, upstream count tables and catalogs, original differential statistics, full GREAT results, original HOMER known-motif tables/BEDs, saved WGCNA networks, selected integration bundles, compact multiome assays and score/control checkpoints, and original resampling distributions. These large files are excluded from ordinary Git commits and need an accompanying archived data release.

## External inputs for optional reruns

| Step | Required external input |
|---|---|
| DiffBind | Original BAMs and indices; set `CORTEX_BAM_DIR`. Portable metadata retains BAM basenames. All selected samples must be available. |
| HOMER | Installed `findMotifsGenome.pl` and mm10 genome; optionally set `CORTEX_HOMER_BIN`. Target and background BEDs are bundled. |
| GREAT | Working rGREAT service and its annotation release. Default figures use archived results. |
| RNA topGO | Explicit analyst-supplied configuration CSV via `CORTEX_TOPGO_CONFIG`: Population, Direction, P_field, P_cutoff, abs_logFC_min, node_size. No historical Figure 5d configuration is inferred. |
| Variant analysis | Original variant TSV and hg38-to-mm10 chain; set `CORTEX_VARIANTS_TSV` and `CORTEX_LIFTOVER_CHAIN`. |
| Orthology | Original orthogene database release was not recorded. New mappings must be saved and checked against the original gene-overlap counts. |

## Details that remain unresolved

1. **Figure 5d:** Original topGO numerical results, plotting selections, gene universe, thresholds, algorithm settings and annotation version were not recovered. The explicit mouse topGO implementation is adapted from Neurona and is not certified as the historical Figure 5d run.
2. **Bulk RNA differential fitting:** Original litter/sex-adjusted result tables are retained, but the original fitting script, design matrix and filter/normalization settings were not recovered. A model is not reverse-engineered from the result filenames and presented as original code.
3. **Grey-module splitting:** Recovered annotation code uses hierarchical dendrogram cuts; manuscript wording refers to k-means. The notebooks show the recovered operation and retain final interval assignments.
4. **HOMER:** Selected `motifFindingParameters` records use other modules as background. A GC-matched global-background branch also exists historically, but has not been established as the source of the selected known-motif matrices. Native HOMER `Log P-value` is ln(P); the stored display is −ln(P), not −log10(P).
5. **Figure 6f and selected graphics:** The saved original donut input counts 484 down / 250 up gene–peak overlaps, whereas final artwork uses 362 / 281. Final Figure 4e network version, GENOA fitting details and some final heatmap selections remain uncertified. `provenance/panel_coverage.csv` describes panel-specific limitations.
6. **Historical environments:** A complete historical package lockfile and all database versions were not supplied. Newly captured session information is not represented as proof of those historical versions.

Default rendering preserves archived numerical results. It does not query live annotation services, refit differential models, construct new networks or draw new null simulations.


### Figure revision audit (2026-11-05)

The manuscript-style draft reserves Figure 6f: neither the recovered nonVZ gene-body overlap analysis (484 down / 250 up), nor the VZ alternative (316 / 102), reproduces the final 362 / 281 totals. Unique-peak and unique-gene counting also do not establish the final input. The recovered donut is a separate diagnostic export, not the final manuscript result.

Figure 6b-prime, Figure 6o and the bulk-RNA panels of Extended Data 7 use the explicit 15-library metadata table recovered from `GEO_submission.Rmd` (lines 3540–3556). The filename suffixes reverse between litters: minus=VZ/FTpos in litter 1, plus=VZ/FTpos in litter 2. Population, RNA-derived sex and population-specific DE annotations were corrected in the plotting/source tables. All counts, lengths, library totals and RPKM measurements are unchanged. Every one of the 285 affected RPKM rows independently agrees with the original featureCounts file. No plotted outliers or analyzed libraries were removed.

Figure 6o displays FTneg for Tbr1, Bcl11b, Wnt7b, Nrxn2 and Shank1, and FTpos for Slc17a6. Correct group sizes are FTpos WT=4/cKO=4 and FTneg WT=3/cKO=4. The high Shank1 cKO observation (22.2913 RPKM, Kdm6b_08plus_mut_S10) is retained and is also present in the manuscript PDF. The original author analysis uses 15 libraries; this correction does not change that inclusion set. `provenance/RNA_sample_metadata_validation.json` records the source and numerical checks.

Figure 2d-prime retains verified Ac9 and Ac14 source values, whose displayed heights appear exchanged in the manuscript artwork. Figure 4e uses the original dynamic-only Integration2 rules on the bundled overlap catalog; canonical labels and Illustrator layout can differ from the draft.


### Module-selection and palette follow-up

Figure 2b now uses the exact nine-module manuscript selection: Ac1, Ac3, Ac5, Ac6, Ac7, Ac9, Ac10, Ac11 and Ac14. The full source matrix is retained; excluded columns are not deleted from source data. Ac16 remains present in VISTA and network panels where it is shown in the manuscript. Figure 4c-prime was corrected to the existing AC-blue/Ac7 chord export (n=4794); the previously selected AC-yellow/Ac9 export was the wrong focal module.

Extended Data 5l: the original all-significant summary has 22 raw module groups. The final plotted subset now uses the recovered positive-DA filter and native pooling, producing the 16 displayed groups without arbitrary module removal. Figure 4e canonical source labels differ from Illustrator relabeling. The detailed module audit is in Figure_drafts/Module_panel_audit.csv. Figure 3 category bars and cluster strips now use the same categories/palette as its UMAP points, without changing scores, coordinates or cell inclusion.

Extended Data 2b uses log2 ratios of footprint activity (Protection Score + TC), retaining the original P values and all 528 raw motif records per contrast. OTX1 has non-positive activity in two comparisons and therefore no defined log2 ratio; no pseudocount was added.



Extended Data Figure 5e/f now have compact plotting checkpoints derived from saved AC/ME3 positive peak-to-gene links and the exact V2 SCT layer. The computation applies the explicitly supplied Figure 3 scoring method to these two assays. The supplied Nick notebook does not establish the historical supplemental marker count or all supplemental display choices. Notebook 11 records those choices and contains the complete regeneration path. Panels 5g/i/k now use the original saved SCT scaled residuals, interval TF-IDF layer and peak TF-IDF layer extracted from the full V2 object. Every one of the 11,614 interval marker features and 12,237 peak marker features matches its original assay exactly. The compact extraction checkpoint is included under objects/multiome. No marker test or normalization model was refitted.

Extended Data 5d displays all 25,218 stored interval rows as five equal-height module strips. Display row ordering is explicit in notebook 11 and is not asserted to reproduce the original Illustrator ordering. Extended Data 6a uses all saved TM3-linked interval means, fixed developmental columns, row z scores and a complete-linkage display dendrogram. No network or model was refitted.


Extended Data 5l is resolved as a plotting-definition discrepancy. The previous summary counted both increased and decreased accessibility with the unpooled bulk module catalog. The corrected plot selects positive DA intervals (adjusted P < 0.001, average log2FC > 0.5, pct.1 >= 0.25), applies Nic's module pooling, and divides by retained intervals within each cluster. All 26,020 retained cluster/interval records contribute to the 16 displayed modules. The original all-significant summary is preserved in the private local archive. The Excel summary, canonical CSV and RDS now match the corrected derived table, 5l_Positive_DA_module_proportions.csv. The historical plotting statement was not recovered, and nearby P-value/boundary choices cannot be uniquely distinguished from a raster, but the reconstructed definition matches the native bar image within approximately two pixels.
