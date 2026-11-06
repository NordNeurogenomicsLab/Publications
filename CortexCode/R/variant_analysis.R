# Original final balanced-downsampling implementation.
# Source: Var_enr_short_obj_oriented.Rmd, ME3_Exact_Enrichment_Run.
# No analysis is executed by sourcing this class definition.

ME3_Exact_Enrichment_Run <- R6::R6Class(
  classname = "ME3_Exact_Enrichment_Run",
  public = list(

    paths = NULL,
    cfg = NULL,
    log_df = NULL,

    me3_df_raw = NULL,
    me3_regions = NULL,
    atac_gr = NULL,
    genoa_gr = NULL,

    var_human_df = NULL,
    var_human_gr = NULL,
    var_lifted_gr = NULL,
    var_filtered_gr = NULL,
    cases_gr = NULL,
    ctrls_gr = NULL,

    case_overlap_table = NULL,
    control_overlap_table = NULL,
    per_region_counts = NULL,

    null_reps = NULL,
    melted_stats = NULL,
    p_values = NULL,
    violin_plot = NULL,
    violin_plot_mut1 = NULL,

    initialize = function(
      me3_diffbind_path,
      atac_path,
      genoa_path,
      variants_tsv_path,
      liftover_chain_path,
      n_null_reps = 1000L,
      downsample_to = 5000L,
      downsample_frac = 0.8,
      rand_seed = 1234L,
      verbose = TRUE
    ) {
      self$paths <- list(
        me3_diffbind_path   = me3_diffbind_path,
        atac_path           = atac_path,
        genoa_path          = genoa_path,
        variants_tsv_path   = variants_tsv_path,
        liftover_chain_path = liftover_chain_path
      )

      self$cfg <- list(
        variant_filter = "atac_and_genoa",
        me3_p_thr = 0.05,
        me3_log2fc_thr = 0,
        me3_expand_lr_bp = 0L,
        atac_expand_lr_bp = 0L,
        merge_overlapping_ME3_intervals = TRUE,
        max_me3_region_size = 50000L,
        thresholds_k = c(0L, 1L, 2L),   # >0, >1, >2 mutations per region
        n_null_reps = as.integer(n_null_reps),
        downsample_to = as.integer(downsample_to),
        downsample_frac = as.numeric(downsample_frac),
        rand_seed = as.integer(rand_seed),
        verbose = isTRUE(verbose)
      )

      self$log_df <- data.frame(
        step = character(),
        message = character(),
        stringsAsFactors = FALSE
      )
    },

    logit = function(step, message) {
      self$log_df <- bind_rows(
        self$log_df,
        data.frame(step = step,
          message = message,
          stringsAsFactors = FALSE)
      )
    },

    compute_log2fc = function(df) {
      df <- as.data.frame(df)

      if ("unshrunk_l2fc" %in% names(df)) {
        df$log2FC <- suppressWarnings(as.numeric(df$unshrunk_l2fc))
        return(df)
      }

      if (all(c("Conc_expr", "Conc_ctrl") %in% names(df))) {
        df$log2FC <- suppressWarnings(as.numeric(df$Conc_expr) - as.numeric(df$Conc_ctrl))
        return(df)
      }

      if (all(c("Conc_Kdm6b_expr", "Conc_Kdm6b_ctrl") %in% names(df))) {
        df$log2FC <- suppressWarnings(as.numeric(df$Conc_Kdm6b_expr) - as.numeric(df$Conc_Kdm6b_ctrl))
        return(df)
      }

      stop("Cannot compute log2FC from supplied ME3 table.")
    },

    expand_lr = function(gr, bp) {
      bp <- as.integer(bp)
      if (is.na(bp) || bp <= 0L) return(gr)

      gr2 <- GenomicRanges::GRanges(
        seqnames = seqnames(gr),
        ranges = IRanges::IRanges(
          start = pmax(1L, start(gr) - bp),
          end   = end(gr) + bp
        ),
        strand = strand(gr)
      )
      mcols(gr2) <- mcols(gr)
      genome(gr2) <- genome(gr)
      gr2
    },

    harmonize_pair = function(x, y) {
      common <- intersect(seqlevels(x), seqlevels(y))
      list(
        x = keepSeqlevels(x, common, pruning.mode = "coarse"),

        y = keepSeqlevels(y,
        common,
        pruning.mode = "coarse")
      )
    },

    read_variants_robust = function(path) {
      needed <- c("chr","start","ID","REF","ALT","TYPE","LEN","FAM","SAMPLE","STATUS")
      df <- readr::read_tsv(path, show_col_types = FALSE)

      if (all(needed %in% names(df))) {
        return(
          df %>%
            select(all_of(needed)) %>%
            mutate(
              start = as.integer(start),
              LEN   = as.integer(LEN),
              end   = start + LEN
            ) %>%
            select(chr, start, end, ID, REF, ALT, TYPE, LEN, FAM, SAMPLE, STATUS)
        )
      }

      if (ncol(df) == 10) {
        colnames(df) <- needed
        df$start <- as.integer(df$start)
        df$LEN   <- as.integer(df$LEN)
        df$end   <- df$start + df$LEN
        return(df[, c("chr","start","end","ID","REF","ALT","TYPE","LEN","FAM","SAMPLE","STATUS")])
      }

      stop("Variant TSV has unexpected format.")
    },

    load_me3_regions = function() {
      self$logit("ME3", "Loading ME3 DiffBind table and selecting regions")

      df <- read.csv(self$paths$me3_diffbind_path, stringsAsFactors = FALSE)
      df <- self$compute_log2fc(df)

      df$p.value <- suppressWarnings(as.numeric(df$p.value))
      df$width   <- suppressWarnings(as.integer(df$width))

      self$me3_df_raw <- df

      df2 <- df %>%
        filter(
          is.finite(p.value), p.value < self$cfg$me3_p_thr,
          is.finite(log2FC),  log2FC  > self$cfg$me3_log2fc_thr,
          is.finite(width),   width   < self$cfg$max_me3_region_size
        )

      if (!nrow(df2)) stop("No ME3 intervals pass filters.")

      gr <- GenomicRanges::GRanges(
        seqnames = df2$seqnames,
        ranges   = IRanges::IRanges(as.integer(df2$start), as.integer(df2$end)),
        strand   = if ("strand" %in% names(df2)) df2$strand else "*"
      )
      genome(gr) <- "mm10"
      mcols(gr)$p.value <- df2$p.value
      mcols(gr)$log2FC  <- df2$log2FC

      gr <- self$expand_lr(gr, self$cfg$me3_expand_lr_bp)
      gr <- gr[width(gr) < self$cfg$max_me3_region_size]

      if (isTRUE(self$cfg$merge_overlapping_ME3_intervals)) {
        gr <- GenomicRanges::reduce(gr, ignore.strand = TRUE)
        genome(gr) <- "mm10"
      }

      mcols(gr)$region_id <- paste0("ME3_", seq_along(gr))
      self$me3_regions <- gr

      self$logit("ME3", paste0("Selected ", length(gr), " ME3 loci"))
      invisible(gr)
    },

    load_atac_and_genoa = function() {
      self$logit("SETS", "Loading ATAC and GENOA")

      atac <- read.csv(self$paths$atac_path, stringsAsFactors = FALSE)
      gr_atac <- GenomicRanges::GRanges(
        seqnames = atac$seqnames,
        ranges   = IRanges::IRanges(as.integer(atac$start), as.integer(atac$end)),
        strand   = "*"
      )
      genome(gr_atac) <- "mm10"
      gr_atac <- self$expand_lr(gr_atac, self$cfg$atac_expand_lr_bp)
      self$atac_gr <- gr_atac

      genoa <- read.table(self$paths$genoa_path, stringsAsFactors = FALSE)
      colnames(genoa)[1:3] <- c("chr", "start", "end")
      gr_genoa <- GenomicRanges::GRanges(
        seqnames = genoa$chr,
        ranges   = IRanges::IRanges(as.integer(genoa$start), as.integer(genoa$end)),
        strand   = "*"
      )
      genome(gr_genoa) <- "mm10"
      self$genoa_gr <- gr_genoa

      self$logit("SETS",
        paste0("ATAC=", length(self$atac_gr), " GENOA=", length(self$genoa_gr)))
      invisible(TRUE)
    },

    load_and_filter_variants = function() {
      self$logit("VAR", "Loading, lifting over, and filtering variants")

      var_df <- self$read_variants_robust(self$paths$variants_tsv_path)
      self$var_human_df <- var_df
      self$var_human_gr <- GenomicRanges::GRanges(var_df)

      chain <- rtracklayer::import.chain(self$paths$liftover_chain_path)
      lifted <- unlist(rtracklayer::liftOver(self$var_human_gr, chain))
      genome(lifted) <- "mm10"
      self$var_lifted_gr <- lifted

      v0 <- lifted[lifted$STATUS != "missing"]
      genome(v0) <- "mm10"
      v0$..idx.. <- seq_along(v0)

      # exact requested filter: atac_and_genoa
      h1 <- self$harmonize_pair(v0, self$atac_gr)
      hits1 <- suppressWarnings(findOverlaps(h1$x, h1$y, ignore.strand = TRUE))
      hit_atac <- rep(FALSE, length(v0))
      hit_atac[h1$x$..idx..[unique(queryHits(hits1))]] <- TRUE

      h2 <- self$harmonize_pair(v0, self$genoa_gr)
      hits2 <- suppressWarnings(findOverlaps(h2$x, h2$y, ignore.strand = TRUE))
      hit_genoa <- rep(FALSE, length(v0))
      hit_genoa[h2$x$..idx..[unique(queryHits(hits2))]] <- TRUE

      keep <- hit_atac & hit_genoa
      v0$..idx.. <- NULL

      v_filt <- v0[keep]
      self$var_filtered_gr <- v_filt
      self$cases_gr <- v_filt[v_filt$STATUS == "case"]
      self$ctrls_gr <- v_filt[v_filt$STATUS == "control"]

      self$logit(
        "VAR",
        paste0(
          "Filtered variants: case=", length(self$cases_gr),
          " control=", length(self$ctrls_gr)
        )
      )
      invisible(v_filt)
    },

    build_overlap_tables = function() {
      self$logit("OVERLAP", "Building overlap tables against final ME3 loci")

      # harmonize for case
      h_case <- self$harmonize_pair(self$me3_regions, self$cases_gr)
      regs_case <- h_case$x
      cases2 <- h_case$y

      ov_case <- suppressWarnings(findOverlaps(regs_case, cases2, ignore.strand = TRUE))
      case_tbl <- tibble(
        region_id      = mcols(regs_case)$region_id[queryHits(ov_case)],
        seqnames       = as.character(seqnames(regs_case))[queryHits(ov_case)],
        start          = start(regs_case)[queryHits(ov_case)],
        end            = end(regs_case)[queryHits(ov_case)],
        width          = width(regs_case)[queryHits(ov_case)],
        variant_chr    = as.character(seqnames(cases2))[subjectHits(ov_case)],
        variant_start  = start(cases2)[subjectHits(ov_case)],
        variant_end    = end(cases2)[subjectHits(ov_case)],
        ID             = mcols(cases2)$ID[subjectHits(ov_case)],
        REF            = mcols(cases2)$REF[subjectHits(ov_case)],
        ALT            = mcols(cases2)$ALT[subjectHits(ov_case)],
        TYPE           = mcols(cases2)$TYPE[subjectHits(ov_case)],
        LEN            = mcols(cases2)$LEN[subjectHits(ov_case)],
        FAM            = mcols(cases2)$FAM[subjectHits(ov_case)],
        SAMPLE         = mcols(cases2)$SAMPLE[subjectHits(ov_case)],
        STATUS         = mcols(cases2)$STATUS[subjectHits(ov_case)]
      )

      # harmonize for control
      h_ctrl <- self$harmonize_pair(self$me3_regions, self$ctrls_gr)
      regs_ctrl <- h_ctrl$x
      ctrls2 <- h_ctrl$y

      ov_ctrl <- suppressWarnings(findOverlaps(regs_ctrl, ctrls2, ignore.strand = TRUE))
      ctrl_tbl <- tibble(
        region_id      = mcols(regs_ctrl)$region_id[queryHits(ov_ctrl)],
        seqnames       = as.character(seqnames(regs_ctrl))[queryHits(ov_ctrl)],
        start          = start(regs_ctrl)[queryHits(ov_ctrl)],
        end            = end(regs_ctrl)[queryHits(ov_ctrl)],
        width          = width(regs_ctrl)[queryHits(ov_ctrl)],
        variant_chr    = as.character(seqnames(ctrls2))[subjectHits(ov_ctrl)],
        variant_start  = start(ctrls2)[subjectHits(ov_ctrl)],
        variant_end    = end(ctrls2)[subjectHits(ov_ctrl)],
        ID             = mcols(ctrls2)$ID[subjectHits(ov_ctrl)],
        REF            = mcols(ctrls2)$REF[subjectHits(ov_ctrl)],
        ALT            = mcols(ctrls2)$ALT[subjectHits(ov_ctrl)],
        TYPE           = mcols(ctrls2)$TYPE[subjectHits(ov_ctrl)],
        LEN            = mcols(ctrls2)$LEN[subjectHits(ov_ctrl)],
        FAM            = mcols(ctrls2)$FAM[subjectHits(ov_ctrl)],
        SAMPLE         = mcols(ctrls2)$SAMPLE[subjectHits(ov_ctrl)],
        STATUS         = mcols(ctrls2)$STATUS[subjectHits(ov_ctrl)]
      )

      self$case_overlap_table <- case_tbl
      self$control_overlap_table <- ctrl_tbl

      case_counts <- case_tbl %>%
        dplyr::count(region_id, name = "n_case_mut") %>%
        right_join(tibble(region_id = mcols(self$me3_regions)$region_id),
          by = "region_id") %>%
        mutate(n_case_mut = replace_na(n_case_mut, 0L))

      ctrl_counts <- ctrl_tbl %>%
        dplyr::count(region_id, name = "n_ctrl_mut") %>%
        right_join(tibble(region_id = mcols(self$me3_regions)$region_id),
          by = "region_id") %>%
        mutate(n_ctrl_mut = replace_na(n_ctrl_mut, 0L))

      self$per_region_counts <- case_counts %>%
        left_join(ctrl_counts, by = "region_id") %>%
        mutate(n_ctrl_mut = replace_na(n_ctrl_mut, 0L))

      invisible(self$per_region_counts)
    },

    resolve_downsample_to = function() {
      n_min <- min(length(self$cases_gr), length(self$ctrls_gr))
      used_frac <- floor(n_min * self$cfg$downsample_frac)
      ds_used <- min(self$cfg$downsample_to, used_frac)
      if (ds_used < 1L) stop("downsample_to_used < 1")
      as.integer(ds_used)
    },

    calc_pvalues_from_reps = function(reps_df) {
      reps_df %>%
        group_by(k, MinMutPerRegion) %>%
        summarise(
          mean_n_case      = mean(n_case),
          mean_n_ctrl      = mean(n_ctrl),
          median_n_case    = median(n_case),
          median_n_ctrl    = median(n_ctrl),
          sd_n_case        = sd(n_case),
          sd_n_ctrl        = sd(n_ctrl),
          mean_diff        = mean(diff),
          median_diff      = median(diff),
          sd_diff          = sd(diff),
          enrichment_ratio = ifelse(mean(n_ctrl) > 0, mean(n_case) / mean(n_ctrl), NA_real_),
          P_null1_1sided   = (sum(diff <= 0) + 1) / (n() + 1),
          P_null1_2sided   = pmin(
            1,
            2 * pmin(
              (sum(diff <= 0) + 1) / (n() + 1),
              1 - ((sum(diff <= 0) + 1) / (n() + 1))
            )
          ),
          .groups = "drop"
        ) %>%
        arrange(k)
    },

    build_violin_plot_all = function(stat_long) {
      ggplot(
        stat_long,
        aes(x = factor(MinMutPerRegion),
          y = n_regions,
          fill = group)
      ) +
        geom_violin(trim = FALSE,
          alpha = 0.5,
          position = position_dodge(width = 0.8)) +
        labs(
          title = "Null resampling distributions of ME3-overlapped regions",
          subtitle = "Exact condition: atac_and_genoa, p<0.05, log2FC>0, ME3 expand=0, ATAC expand=0, merge=TRUE",

          x = "Minimum mutations per region",
          y = "# regions passing threshold",
          fill = NULL
        ) +
        theme_bw(base_size = 12)
    },

    build_violin_plot_mut1 = function(stat_long, p_values_tbl) {
      plot_df <- stat_long %>%
        filter(MinMutPerRegion == 1) %>%
        mutate(group = factor(group,
          levels = c("control", "case")))

      p_row <- p_values_tbl %>% filter(MinMutPerRegion == 1)
      p_lab <- if (nrow(p_row) == 1) {
          paste0("P = ", sprintf("%.5f", p_row$P_null1_1sided), " (Permutation test)")
            } else {
              "P = NA"
    }

      ggplot(plot_df,
        aes(x = group, y = n_regions, fill = group, color = group)) +
        geom_violin(trim = FALSE,
          alpha = 0.22,
          linewidth = 0.7,
          width = 0.7) +
        geom_jitter(width = 0.16,
          height = 0,
          alpha = 0.18,
          size = 1.8,
          stroke = 0) +
        scale_fill_manual(values = c("control" = "#00BFC4", "case" = "#F8766D")) +
        scale_color_manual(values = c("control" = "#00BFC4", "case" = "#F8766D")) +
        labs(
          title = paste0(">0 mut. per region\n[", p_lab, "]"),
          x = NULL,
          y = "Regions with intersecting mutations [counts]"
        ) +
        theme_bw(base_size = 18) +
        theme(
          legend.position = "none",
          plot.title = element_text(hjust = 0.5, face = "plain"),

          axis.text.x = element_text(size = 16),

          axis.title.y = element_text(size = 18),

          panel.grid.minor = element_blank()
        )
    },

    run_null = function() {
      self$logit("NULL", "Running case/control resampling null")

      set.seed(self$cfg$rand_seed)

      h_case <- self$harmonize_pair(self$cases_gr, self$me3_regions)
      cases2 <- h_case$x
      regs2  <- h_case$y

      h_ctrl <- self$harmonize_pair(self$ctrls_gr, regs2)
      ctrls2 <- h_ctrl$x
      regs2  <- h_ctrl$y

      ds_used <- self$resolve_downsample_to()
      ks <- self$cfg$thresholds_k

      one_rep <- function(rep_id) {
        case_hits <- countOverlaps(
          regs2,
          sample(cases2, ds_used, replace = FALSE),
          ignore.strand = TRUE
        )
        ctrl_hits <- countOverlaps(
          regs2,
          sample(ctrls2, ds_used, replace = FALSE),
          ignore.strand = TRUE
        )

        bind_rows(lapply(ks, function(k) {
          tibble(
            rep = rep_id,
            k = k,
            MinMutPerRegion = k + 1L,
            n_case = sum(case_hits > k),
            n_ctrl = sum(ctrl_hits > k),
            diff = sum(case_hits > k) - sum(ctrl_hits > k)
          )
        }))
      }

      reps_df <- bind_rows(lapply(seq_len(self$cfg$n_null_reps), one_rep))
      self$null_reps <- reps_df

      stat_long <- reps_df %>%
        pivot_longer(
          cols = c(n_case, n_ctrl),
          names_to = "group",
          values_to = "n_regions"
        ) %>%
        mutate(
          group = recode(group,
          n_case = "case",
          n_ctrl = "control")
        )

      self$melted_stats <- stat_long %>%
        group_by(group, k, MinMutPerRegion) %>%
        summarise(
          mean_n_regions   = mean(n_regions),
          median_n_regions = median(n_regions),
          sd_n_regions     = sd(n_regions),
          min_n_regions    = min(n_regions),
          max_n_regions    = max(n_regions),
          .groups = "drop"
        ) %>%
        arrange(k, group)

      self$p_values <- self$calc_pvalues_from_reps(reps_df)

      self$violin_plot <- self$build_violin_plot_all(stat_long)
      self$violin_plot_mut1 <- self$build_violin_plot_mut1(stat_long, self$p_values)

      self$logit("NULL", paste0("Completed ", self$cfg$n_null_reps, " null reps"))
      invisible(TRUE)
    },

    run_all = function() {
      self$load_me3_regions()
      self$load_atac_and_genoa()
      self$load_and_filter_variants()
      self$build_overlap_tables()
      self$run_null()
      invisible(TRUE)
    },

    get_results = function() {
      list(
        config = self$cfg,
        me3_regions = self$me3_regions,
        case_overlap_table = self$case_overlap_table,
        control_overlap_table = self$control_overlap_table,
        per_region_counts = self$per_region_counts,
        null_reps = self$null_reps,
        melted_stats = self$melted_stats,
        p_values = self$p_values,
        violin_plot = self$violin_plot,
        violin_plot_mut1 = self$violin_plot_mut1,
        log = self$log_df
      )
    }
  )
)
