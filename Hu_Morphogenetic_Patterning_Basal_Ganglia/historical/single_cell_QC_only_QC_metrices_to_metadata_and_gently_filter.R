
 single_cell_QC <- function(sample){
    
    # Load RNA counts
    counts <- Read10X(paste0("cellranger_output/", sample,  "/outs/filtered_feature_bc_matrix"))
    
    # Create a Seurat object containing the RNA data
    pbmc <- CreateSeuratObject(
        counts = counts,
        assay = "RNA"
    )
    
    # Check names of the mitochondrial genes
    # rownames(pbmc)[grepl("^mt-", rownames(pbmc), ignore.case = TRUE)]
    
    # Mitochondrial genes: 
    DefaultAssay(pbmc) <- "RNA"
    pbmc[["percent.mt"]] <- PercentageFeatureSet(pbmc, pattern = "^mt-")
    
    # Violin plot - pre-filtering        
    p1 <- VlnPlot(
        object = pbmc,
        features = c("nCount_RNA", "nFeature_RNA", "percent.mt"),
        ncol = 3,
        pt.size = 0.1)& 
        theme(axis.text.x = element_blank())&
        plot_annotation(title = paste0(sample, ", Cellranger filtered cells.", " Features: ", nrow(pbmc),
                                       " Cells: ", ncol(pbmc)), 
                        theme = theme(plot.title = element_text(hjust = 0.5, size = 16, face = "bold")))
    
    
    p3 <- ggplot(pbmc@meta.data, aes(x = log2(nFeature_RNA),  y = log(nCount_RNA)))+
        geom_point(alpha = 0.1)+
        theme_bw()+
        labs(title = "Cellranger filtered cells")
    
    
    
    
    # Define QC filtering thresholds using 3 median absolute deviations. 
    # This idea borrowed from the scater package
    
    
    
    #count_threshold_high <- median(pbmc@meta.data$nCount_RNA) + mad_distance*mad(pbmc@meta.data$nCount_RNA)
    #features_threshold_high <- median(pbmc@meta.data$nFeature_RNA) + mad_distance*mad(pbmc@meta.data$nFeature_RNA)
    #mito_threshold_high <- median(pbmc@meta.data$percent.mt) + mad_distance*mad(pbmc@meta.data$percent.mt)
    
    #
    mad_distance = 5
    
    pbmc[["nCount_RNA_abv_1000"]] <- pbmc@meta.data$nCount_RNA > 1000
    pbmc[["nCount_RNA_less_than_5_MADs"]] <- pbmc@meta.data$nCount_RNA < 
                                             median(pbmc@meta.data$nCount_RNA) + mad_distance*mad(pbmc@meta.data$nCount_RNA)
    
    pbmc[["nCount_RNA_abv_1000_and_less_than_5_MADs"]] <- pbmc[["nCount_RNA_abv_1000"]] & 
                                                          pbmc[["nCount_RNA_less_than_5_MADs"]]
    
    
    pbmc[["nFeature_RNA_abv_200"]] <- pbmc@meta.data$nFeature_RNA > 200
    pbmc[["nFeature_RNA_less_than_5_MADs"]] <- pbmc@meta.data$nFeature_RNA < 
                                               median(pbmc@meta.data$nFeature_RNA) + mad_distance*mad(pbmc@meta.data$nFeature_RNA) 
    
    pbmc[["nFeature_RNA_abv_200_and_less_than_5_MADs"]] <- pbmc[["nFeature_RNA_abv_200"]] &
                                                           pbmc[["nFeature_RNA_less_than_5_MADs"]]
    
    pbmc[["percent.mt_less_than_15"]] <- pbmc@meta.data$percent.mt < 15
    
    
    ######
    pbmc[["Cell_passing_QC"]] <- pbmc[["nCount_RNA_abv_1000_and_less_than_5_MADs"]] & 
                                 pbmc[["nFeature_RNA_abv_200_and_less_than_5_MADs"]] &
                                 pbmc[["percent.mt_less_than_15"]]
        
    
    
    ## Filter out low quality cells
    pbmc <- subset(
        x = pbmc,
        subset =
            #nCount_RNA < count_threshold_high &
            nCount_RNA > 1000 &
            #nFeature_RNA < features_threshold_high &  
            #nFeature_RNA > 200 &
            percent.mt < 20 
    )
    
    
    # Violin plot - post-filtering    
    p2 <- VlnPlot(
        object = pbmc,
        features = c("nCount_RNA", "nFeature_RNA", "percent.mt"),
        ncol = 3,
        pt.size = 0.1)& 
        theme(axis.text.x = element_blank())&
        plot_annotation(title = paste0(sample, ", Quality filtered cells.", " Features: ", nrow(pbmc),
                                       " Cells: ", ncol(pbmc)), 
                        theme = theme(plot.title = element_text(hjust = 0.5, size = 16, face = "bold")))
    
    
    p4 <- ggplot(pbmc@meta.data, aes(x = log2(nFeature_RNA),  y = log(nCount_RNA)))+
        geom_point(alpha = 0.1)+
        theme_bw()+
        labs(title = "Quality filtered cells")
    
    
    
    # Save post filtering object    
    saveRDS(pbmc, file = paste0("./QCed Seurat objects/", sample, "_QCed_Seurat_obj.rds"))
    
    
   top <-  plot_grid(p1, p2, ncol = 1)
    
   bottom <- plot_grid(p3, p4, ncol = 2)
   
   plot_grid(top, bottom, ncol = 1)
    
    
    }