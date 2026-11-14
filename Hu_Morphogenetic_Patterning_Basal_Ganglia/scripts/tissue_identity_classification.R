# Marker-guided tissue identity: original VZ H11 vs H13 notebook (July 2023).
# POA markers Igfbp5/Dbx1/Hmx3 span adjacent progenitor and neuronal identities;
# Spry1/Fgf8/Zic4 identify septal contaminants. A union across markers avoids
# requiring simultaneous detection of all three genes in sparse single-cell data.
# Initial labels use SCT expression; SVM fitting uses original RNA UMI counts.
# The source deliberately used RNA counts to reduce low-expression classifications.
# The random 50/50 split prevents a sequential split from separating whole libraries.
# The exact manuscript retained-cell table remains available for deterministic
# count-based analysis; refitting classifiers must start with pre-exclusion cells.
classify_tissue_identity <- function(o, output_dir) {
  stopifnot(requireNamespace("scAnnotatR", quietly = TRUE))
  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
  positive <- function(gene, q) {
    expr <- GetAssayData(o, assay = "SCT", slot = "data")[gene, , drop = TRUE]
    expr <- expr[expr > 0]
    if (!length(expr)) stop("No positive expression for tissue marker ", gene)
    names(expr)[expr > quantile(expr, q)]
  }
  poa_genes <- c("Igfbp5", "Dbx1", "Hmx3")
  septum_genes <- c("Spry1", "Fgf8", "Zic4")
  poa_seed <- unique(unlist(lapply(poa_genes, positive, q = 0.8)))
  septum_seed <- unique(unlist(lapply(septum_genes, positive, q = 0.5)))
  o$Most_likely_cell_type <- ifelse(colnames(o) %in% poa_seed, "POA", "MGE")
  fit <- function(o, genes, target) {
    set.seed(1234)
    train_idx <- sample(seq_len(ncol(o)), ncol(o) %/% 2, replace = FALSE)
    test_idx <- setdiff(seq_len(ncol(o)), train_idx)
    # scAnnotatR renamed slot to layer in recent releases. Specify counts with
    # either API; changing package defaults must not switch the expression input.
    expression_arg <- if ("layer" %in% names(formals(scAnnotatR::train_classifier))) "layer" else "slot"
    args <- list(train_obj = o[, train_idx], assay = "RNA", cell_type = target,
      marker_genes = genes, tag_slot = "Most_likely_cell_type", zscore = TRUE)
    args[[expression_arg]] <- "counts"
    model <- do.call(scAnnotatR::train_classifier, args)
    predict <- function(z) {
      args <- list(classifier = model, test_obj = z, assay = "RNA", tag_slot = "Most_likely_cell_type")
      expression_arg <- if ("layer" %in% names(formals(scAnnotatR::test_classifier))) "layer" else "slot"
      args[[expression_arg]] <- "counts"
      do.call(scAnnotatR::test_classifier, args)
    }
    heldout <- predict(o[, test_idx]); total <- predict(o)
    # scAnnotatR's binary SVM balances the training classes and estimates class
    # probabilities. Save the actual fitted caret model/tuning and held-out AUC,
    # rather than hard-coding the notebook's reported validation values.
    results <- list(model = model, heldout = heldout, total = total,
      train_cells = colnames(o)[train_idx], test_cells = colnames(o)[test_idx])
    saveRDS(results, file.path(output_dir, paste0(target, "_classifier_and_predictions.rds")))
    writeLines(capture.output(scAnnotatR::caret_model(model)), file.path(output_dir, paste0(target, "_model_settings.txt")))
    rownames(total$pred)[total$pred$class == 1]
  }
  poa <- fit(o, poa_genes, "POA")
  o$Most_likely_cell_type <- ifelse(colnames(o) %in% septum_seed, "Septum", o$Most_likely_cell_type)
  septum <- fit(o, septum_genes, "Septum")
  labels <- setNames(rep("MGE", ncol(o)), colnames(o))
  labels[septum] <- "Septum"; labels[poa] <- "POA"
  labels[intersect(poa, septum)] <- "POA_or_Septum"
  o$SVM_cell_type <- labels[colnames(o)]
  write.table(data.frame(cell_id = colnames(o), tissue = o$SVM_cell_type),
    file.path(output_dir, "tissue_classifications.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
  writeLines(capture.output(sessionInfo()), file.path(output_dir, "classifier_sessionInfo.txt"))
  list(classified = o, MGE = subset(o, subset = SVM_cell_type == "MGE"))
}
