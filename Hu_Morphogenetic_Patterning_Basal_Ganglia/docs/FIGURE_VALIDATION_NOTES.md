# Figure validation notes

Saved original integrated object and saved MZ object both pass validity checks. Counts and mean detected genes match the manuscript captions exactly to rounding.

Recreated UMAPs use stored coordinates; feature plots use saved SCT values; volcano plots use existing DE tables with the original script's sign reversal. Manual labels, point sizes/order and expression scaling can differ. No integration, clustering or DE tests were repeated.

The older Dev_classes column is missing for the entire later JH4327 dataset. Selecting solely on that column therefore excludes those later cells. A second VZ recreation uses cluster IDs 2,7,5,14 explicitly documented as VZ in the original notebook (1950 cells: 1512 E11 and 438 E13). This is a documented candidate subset, not proof of the exact manually selected manuscript subset. VZ plot differences remain visible in the deck.

Figure 4D has 15 saved clusters but 8 named anatomical groups. MZ_candidate_annotations.tsv reports expression-based candidate phenotypes, not a recovered eight-group mapping. Markers and age proportions do not uniquely identify every anatomical type. Several cluster combinations numerically resemble the printed proportions, which is insufficient evidence to assign identities. The displayed eight E11 bar labels sum to about 90.9% and E13 labels to about 92%; they do not form a complete 100% partition as printed. This suggests omitted/unassigned cells or a different categorization; it must be clarified rather than renormalized silently.

The comparison deck intentionally shows the known mismatch and the verified 15-cluster proportions. Final anatomical cell labels remain unresolved. Histology, diagrams and non-single-cell experiments are excluded from recreation. Supplementary PDF panels were not provided.
