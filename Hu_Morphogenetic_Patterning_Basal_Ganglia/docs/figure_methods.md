# Single-cell figure methods

| Manuscript panels | Analysis input | Output |
| --- | --- | --- |
| Figure 3A | Global SCT-integrated MGE cells | Developmental-age UMAP and Hes1/Ascl1/Dcx expression |
| Figure 3B | Broad VZ E11.5 versus E13.5 SCT differential expression | Volcano plot with E13.5-positive horizontal axis |
| Figure 3C–D, computational panels | Global VZ cell subset | Id4, Mest, Dlx2, and Nfib expression split by age |
| Figure 4A | Broad MZ E11.5 versus E13.5 SCT differential expression | Volcano plot |
| Figure 4B–C, computational panels | Refined neuronal MZ object | Lhx8, Mest, Npy, and Mafb expression split by age |
| Figure 4D and Figure S7A | MZ 52-gene embedding and cluster identities | Age/cluster UMAPs; Figure S7 numerical-cluster proportions |
| Figure S7B′–C′, computational panels | Full-gene SCT expression in MZ | Nr2f1, Id4, Ntrk1, Isl1, Gbx1, Cdh8, and Th expression maps |

Scripts 04–05 specify the expression assay, contrast direction, color scales, label selection, dimensions, and output resolution. Script 03 supplies the corresponding statistical analyses; scripts 01–02 supply the count-to-object methods. The saved-object entry point preserves the manuscript embeddings while the count entry point calculates new embeddings.

The original cluster annotation is a reference input, not a stable identifier produced by any future clustering run. The six-cluster VZ developmental compartment is used consistently for expression maps and broad VZ differential expression. The source’s four-cluster request concerned individual cluster contrasts and is not substituted for the broad compartment definition. Figure 4D's eight anatomical labels still require an author-confirmed key from the 15 saved clusters; the current scripts report the stored numerical identities and complete age-specific denominator tables. Histology panels are outside this code's computational scope.
