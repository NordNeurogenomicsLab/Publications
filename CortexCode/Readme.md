## Manuscript

**Profiling dynamic epigenetic programs in the developing mouse cortex reveals the role of ASD-associated Kdm6b in deep-layer neuron differentiation.**

Athéna R. Ypsilanti<sup>†,*</sup>, Eirene Markenscoff-Papadimitriou<sup>†,*</sup>, Karol Cichewicz<sup>†</sup>, Sean Whalen, Nicolas Seban, Ethan Fenton, Kenneth Lim, Shan Dong, Jada Baldoza, Yurong Wang, Matt Barraza, Fadya Binyameen, Diane Dickel, Len A. Pennacchio, Axel Visel, Stephan J. Sanders, Katherine S. Pollard, Alex S. Nord<sup>*</sup>, John Rubenstein<sup>*</sup>

† Contributed equally: Athéna R. Ypsilanti, Eirene Markenscoff-Papadimitriou and Karol Cichewicz.

\* Corresponding authors: Athéna R. Ypsilanti, Eirene Markenscoff-Papadimitriou, Alex S. Nord and John Rubenstein.

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
| [10_manuscript_graphics.Rmd](10_manuscript_graphics.Rmd) | Panel styling and exports |
| [11_additional_manuscript_panels.Rmd](11_additional_manuscript_panels.Rmd) | Additional computational panels |
| [12_Extended_multiome_heatmaps.Rmd](12_Extended_multiome_heatmaps.Rmd) | Extended Data multiome heatmaps |

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

Shared R functions: [R/](R/).
