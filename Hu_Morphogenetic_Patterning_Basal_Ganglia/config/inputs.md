# Required inputs under JIA_WORK_DIR

- inputs/seurat_objects/MGE_1538_E11_and_E13_with_new_E13.Rds
- inputs/seurat_objects/MZ_reintegrated_on_52_PN_genes_08_15_2025.rds
- inputs/DE_tables/From_previous_integrations/VZ_markers.csv
- inputs/DE_tables/From_previous_integrations/MZ_markers.csv

Create figures/ and audit/. Use rclone copy, preserve sources, verify hashes. Raw FASTQs are not required to render saved plots. The staging script is specific to the existing Jia server and Box layout.
