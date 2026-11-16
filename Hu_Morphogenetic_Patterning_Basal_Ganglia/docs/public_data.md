# Public single-cell analysis data

The [public Box release](https://ucdavis.box.com/s/2ij8dibr0u8djy93kis82okf42xfnjf4) contains **46 scientific files**, **4,405,155,826 bytes** (4.41 GB; 4.10 GiB). Downloads are available without authentication. The original release is dated November 13, 2026; the same scientific files were transferred to Box on November 16, 2026.

Sequencing data: **[GSE354638](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE354638)**, scheduled for public release on **November 18, 2026**.

[File manifest](../config/PUBLIC_DATA_FILES.tsv) · [Hosted manifest](https://ucdavis.box.com/s/fkebyyhtq9pol7dlqijt9a6aomkyevla) · [SHA-256 checksums](https://ucdavis.box.com/s/6h94sngxe2j1b5ryv42wznw4iecodxgh)

## Download

From the manuscript code directory:

```bash
python3 scripts/download_public_data.py --output /path/to/jia_release
export JIA_DATA_DIR=/path/to/jia_release/data
export JIA_OUTPUT_DIR=/path/to/analysis_outputs
```

The downloader checks file sizes and SHA-256 checksums and skips files that already match. To download only analysis inputs, add `--include data`. For result tables and checkpoints, use `--include results`; for the original Tables S1–S4, use `--include supplementary_tables`. Multiple `--include` options may be combined. Individual files can also be downloaded from the links below. Count matrices need all three files in their library directory. See [input configuration and execution](../config/inputs.md) for the analysis commands.

## File contents

`data/objects/` contains the three supplied manuscript Seurat objects, retaining their stored cell assignments, expression assays, and embeddings. `data/counts/` contains RNA UMI counts for the accepted MGE cells; `data/original_counts/` contains the pre-exclusion Cell Ranger-called cells for 1538TE10H11 and 1538TE10H13. The annotation and count manifests use paths relative to `JIA_DATA_DIR`.

`results/count_outputs/` contains newly calculated results and checkpoints from the count-based scripts. `results/saved_outputs/` contains figure data calculated from the supplied manuscript objects. These folders are identified separately because fitted models and numerical cluster labels may differ between runs. The [scientific methods](scientific_methods.md) describe the analysis and annotation requirements.

`TableS1.xlsx`–`TableS4.xlsx` are unchanged source workbooks. The [gene-symbol correction record](supplement_gene_symbol_dates.tsv) documents symbols that Excel converted to dates; analysis CSVs retain gene symbols.

## Analysis inputs

| File | MB | Contents |
| --- | ---: | --- |
| [`data/cell_annotations.tsv`](https://ucdavis.box.com/s/sbwgo4aaimghmxr2betxg5ge75m5u9ol) | 0.53 | Manuscript analysis input |
| [`data/count_inputs.tsv`](https://ucdavis.box.com/s/iugu3cji5iwunleief5ox08ihleg9mde) | 0.00 | Manuscript analysis input |
| [`data/counts/1538TE10H11/barcodes.tsv`](https://ucdavis.box.com/s/sr14460v442rhrkflsexdz9z31qwd9ew) | 0.06 | Manuscript analysis input |
| [`data/counts/1538TE10H11/features.tsv`](https://ucdavis.box.com/s/9rqo2f2dcb2m0p828b7z5oa152em13dx) | 1.01 | Manuscript analysis input |
| [`data/counts/1538TE10H11/matrix.mtx`](https://ucdavis.box.com/s/k985yqhof4aqzy0hmz6cl8hhxydi68ou) | 162.75 | Manuscript analysis input |
| [`data/counts/1538TE10H13/barcodes.tsv`](https://ucdavis.box.com/s/raec2byfhigq2uh7i8mr0z22waii2ihf) | 0.03 | Manuscript analysis input |
| [`data/counts/1538TE10H13/features.tsv`](https://ucdavis.box.com/s/2scq3jsblye5e61dbudufohn3c3rkz41) | 1.01 | Manuscript analysis input |
| [`data/counts/1538TE10H13/matrix.mtx`](https://ucdavis.box.com/s/enn3ivh9dtcpjrhxrov2bzw46v83wbxv) | 58.81 | Manuscript analysis input |
| [`data/counts/JH4327_01_TE10H13/barcodes.tsv`](https://ucdavis.box.com/s/mw2x7160in4dmjpo3k5wr6v3rj12k7so) | 0.06 | Manuscript analysis input |
| [`data/counts/JH4327_01_TE10H13/features.tsv`](https://ucdavis.box.com/s/j4l0guexejmlcqn58b86gm1211olixbf) | 1.01 | Manuscript analysis input |
| [`data/counts/JH4327_01_TE10H13/matrix.mtx`](https://ucdavis.box.com/s/ot7plh0hi910pahzoag35knu4cc2s7e5) | 220.50 | Manuscript analysis input |
| [`data/objects/MGE_1538_E11_and_E13_with_new_E13.Rds`](https://ucdavis.box.com/s/lv3ph1bgx940fckm17ho223awepgt4bl) | 880.36 | Manuscript analysis input |
| [`data/objects/MZ_reintegrated_on_52_PN_genes_08_15_2025.rds`](https://ucdavis.box.com/s/m1yc3nrivfuce626uoxq0b342p7mnde7) | 452.55 | Manuscript analysis input |
| [`data/objects/VZ_SVZ_reintegrated_on_221_progenitor_genes_08_15_2025.rds`](https://ucdavis.box.com/s/5qr5qc7r5fr6onec8eg44trc97q3b4yq) | 378.76 | Manuscript analysis input |
| [`data/original_counts/1538TE10H11/barcodes.tsv.gz`](https://ucdavis.box.com/s/37hlm0a86zeta3g71tg4wzrxeyi2rs1s) | 0.03 | Pre-exclusion Cell Ranger counts |
| [`data/original_counts/1538TE10H11/features.tsv.gz`](https://ucdavis.box.com/s/ti1qqi6yy8mhgo5d01797e4ovf87wm3u) | 0.29 | Pre-exclusion Cell Ranger counts |
| [`data/original_counts/1538TE10H11/matrix.mtx.gz`](https://ucdavis.box.com/s/a5j1lartbktoxp066mawf8t3wvlvdavv) | 81.26 | Pre-exclusion Cell Ranger counts |
| [`data/original_counts/1538TE10H13/barcodes.tsv.gz`](https://ucdavis.box.com/s/gq09m1ho6rtitabdnbgp8ijb06ikx7gq) | 0.02 | Pre-exclusion Cell Ranger counts |
| [`data/original_counts/1538TE10H13/features.tsv.gz`](https://ucdavis.box.com/s/l5xhsz2rjq3kicohxzo6pf0bmdeudycj) | 0.29 | Pre-exclusion Cell Ranger counts |
| [`data/original_counts/1538TE10H13/matrix.mtx.gz`](https://ucdavis.box.com/s/syk73j158lr4shlpmf3ftytgupxnvwhh) | 32.93 | Pre-exclusion Cell Ranger counts |
| [`data/count_inputs_preselection.tsv`](https://ucdavis.box.com/s/cw93q02s3fr5tvseyl9xz00vjkbrr13e) | 0.00 | Pre-exclusion count manifest |

## Count-based analysis results

| File | MB | Contents |
| --- | ---: | --- |
| [`results/count_outputs/MZ_cell_metadata.csv`](https://ucdavis.box.com/s/i23lapw62ioazebew827egzps2yk2jm1) | 0.49 | Count-based analysis result |
| [`results/count_outputs/MZ_cluster_age_proportions.csv`](https://ucdavis.box.com/s/r6ivj77wkqw8yl2r3p3ugltfjklw4hfx) | 0.00 | Count-based analysis result |
| [`results/count_outputs/MZ_cluster_markers.csv`](https://ucdavis.box.com/s/7igddj6xwaln3mimlcfef3p13822hkid) | 1.30 | Count-based analysis result |
| [`results/count_outputs/MZ_markers.csv`](https://ucdavis.box.com/s/u37eld88aj1vid0yfyy33b78bdmk1qgt) | 0.41 | Count-based analysis result |
| [`results/count_outputs/MZ_volcano_data.csv`](https://ucdavis.box.com/s/ccwgtdv9uyr9hcitqed3tr2is6us6d0z) | 0.59 | Count-based analysis result |
| [`results/count_outputs/MZ_within_cluster_age_markers.csv`](https://ucdavis.box.com/s/x4bqyg8yb8bxic0gyoizctkvdl7ddsf8) | 5.39 | Count-based analysis result |
| [`results/count_outputs/SVZ_markers.csv`](https://ucdavis.box.com/s/r436st1jomxf3fdwdzwhuzospdk7ppt6) | 0.38 | Count-based analysis result |
| [`results/count_outputs/VZ_SVZ_cell_metadata.csv`](https://ucdavis.box.com/s/58b75su3jxukj8ctf9onb34jypfndkxw) | 0.65 | Count-based analysis result |
| [`results/count_outputs/VZ_SVZ_cluster_markers.csv`](https://ucdavis.box.com/s/r416l1vbdvw3dx8q2crh5w635nv31at5) | 1.01 | Count-based analysis result |
| [`results/count_outputs/VZ_SVZ_within_cluster_age_markers.csv`](https://ucdavis.box.com/s/zxv7q4az0zbnz3v1naw6illntm0ahk93) | 3.81 | Count-based analysis result |
| [`results/count_outputs/VZ_markers.csv`](https://ucdavis.box.com/s/xkcm7c35yuwploaf7eyuyv78acqxvpk2) | 0.43 | Count-based analysis result |
| [`results/count_outputs/VZ_volcano_data.csv`](https://ucdavis.box.com/s/r0o7y6dt8qihfiwhzitt4haugpkzhydg) | 0.61 | Count-based analysis result |
| [`results/count_outputs/count_import_qc.csv`](https://ucdavis.box.com/s/jte579omg7nlj7n8358f1w0njaat6zwv) | 0.00 | Count-based analysis result |
| [`results/count_outputs/global_cluster_correspondence.csv`](https://ucdavis.box.com/s/82xf733w3llvxhny9t3m2evnj5ojvawh) | 0.00 | Count-based analysis result |
| [`results/count_outputs/objects/MGE_1538_E11_and_E13_with_new_E13.Rds`](https://ucdavis.box.com/s/ip96ifam8xfftqw4rus80otnhsi1mpmj) | 880.27 | Count-based analysis checkpoint |
| [`results/count_outputs/objects/MZ_reintegrated_on_52_PN_genes_08_15_2025.rds`](https://ucdavis.box.com/s/dzxxnff40arokowww42lnuyr2yufct4u) | 455.35 | Count-based analysis checkpoint |
| [`results/count_outputs/objects/VZ_SVZ_cluster_marker_model.rds`](https://ucdavis.box.com/s/vcr1xr80xb26ryf8ddfhicyo2y2ag38l) | 391.49 | Count-based analysis checkpoint |
| [`results/count_outputs/objects/VZ_SVZ_reintegrated_on_221_progenitor_genes_08_15_2025.rds`](https://ucdavis.box.com/s/huz1gmxafzbv8fpp9tcmgvcgsou2acb4) | 389.30 | Count-based analysis checkpoint |

## Manuscript-object figure data

| File | MB | Contents |
| --- | ---: | --- |
| [`results/saved_outputs/MZ_cluster_age_proportions.csv`](https://ucdavis.box.com/s/8atauix1gjwt01mhse3xzrv46554owtj) | 0.00 | Manuscript-object figure data |
| [`results/saved_outputs/MZ_volcano_data.csv`](https://ucdavis.box.com/s/koqxgxo5w567lbcp0lch0opzcaqalvkn) | 0.40 | Manuscript-object figure data |
| [`results/saved_outputs/VZ_volcano_data.csv`](https://ucdavis.box.com/s/0lpoub9mno0gw2gvgaqhy9luvn2nh058) | 0.37 | Manuscript-object figure data |

## Original supplementary tables

| File | MB | Contents |
| --- | ---: | --- |
| [`supplementary_tables/TableS1.xlsx`](https://ucdavis.box.com/s/gpmsg6chz1ydw3l2cwye36j7h2656a86) | 0.30 | Original manuscript supplementary table |
| [`supplementary_tables/TableS2.xlsx`](https://ucdavis.box.com/s/uguea6m91ya4hhwkchfsp0c3374kenqp) | 0.32 | Original manuscript supplementary table |
| [`supplementary_tables/TableS3.xlsx`](https://ucdavis.box.com/s/zcxbuggmd9cdeht9az6kqbbmqip1gp08) | 0.01 | Original manuscript supplementary table |
| [`supplementary_tables/TableS4.xlsx`](https://ucdavis.box.com/s/2pxj4yurqxfv7flmmt32spt22ynnavb7) | 0.01 | Original manuscript supplementary table |

## Release verification

All 46 scientific files retain the original SHA-256 checksums. All 49 files in the Box release (46 scientific files plus the manifest, checksum file, and release description) were checked against the server copies by file size and SHA-1. Anonymous downloads were tested for every file, and complete supplementary-table downloads passed SHA-256 verification on November 16, 2026.
