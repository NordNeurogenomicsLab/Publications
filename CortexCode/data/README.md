# Analysis data

Selected numerical inputs are included in GitHub. [TABLES_ON_GITHUB.csv](TABLES_ON_GITHUB.csv) lists deposited files, sizes and SHA-256 checksums. Larger data files remain excluded from Git.

| Subdirectory | Contents |
|---|---|
| `upstream/WGCNA/` | Original count TSVs, merged histone count objects and fitted networks |
| `upstream/RNA/` | Bulk RNA featureCounts, sample metadata and original DE tables |
| `upstream/multiome/` | Module definitions, peak-to-gene links and pseudobulk counts |
| `upstream/final_tables/` | Final ATAC/AC/ME3 catalogs and differential peak results |
| `upstream/diffbind/` | DiffBind sample sheet |
| `upstream/enrichment/`, `upstream/HOMER/` | Saved GREAT and HOMER results and motif regions |
| `upstream/integration/`, `upstream/GENOA/`, `upstream/variants/` | Integration and variant-analysis inputs |
| `source_tables/` | Figure-specific numerical CSVs |
| `raw_sheets/` | Original histology worksheet grids |

[Exact filenames and input requirements](../INPUTS.md) · [Complete inventory](../config/file_inventory.csv)

Files not listed in TABLES_ON_GITHUB.csv still require a separate input archive.
