############################################################
## SCRIPT 5
## State-specific DEG analysis + PM vs miBrain comparison
##
## Purpose:
##   - run PM state-specific DEG for selected scVI states
##   - support two PM model modes:
##       (1) additive / GA-included model
##       (2) homozygote-only AA vs GG model
##   - compare PM state-specific DEG to matched miBrain state DEG
##   - make state-specific concordance scatter plots
##
## This script is separate from abundance analysis.
## It does NOT:
##   - rebuild pseudobulk
##   - fit abundance models
##   - compute state proportions
############################################################

############################################################
## 0. Configuration
############################################################

config <- list(
  out_dir = file.path(getwd(), "results", "script5_state_specific_deg"),
  
  # PM input object can be either:
  #   (1) already loaded in the environment under pm_object_name, or
  #   (2) read from pm_rds_path if provided.
  pm_object_name = "pb_MG",
  pm_rds_path    = NULL,
  
  # State assays to analyze (edit once here)
  states_to_run = c("M3: SPP1", "M5: PCDH9"),
  
  # Optional matched miBrain state-DEG tables for comparison plots.
  # These should be the saved DEG outputs for the same states.
  # Example names might be:
  #   miBrain_state_spp1_path  = ".../MB_DREAM_M3_SPP1.tsv"
  #   miBrain_state_pcdh9_path = ".../MB_DREAM_M5_PCDH9.tsv"
  miBrain_state_paths = list(
    "M3: SPP1"  = NULL,
    "M5: PCDH9" = NULL
  ),
  
  # Output subfolders
  out_dir_ga    = file.path(getwd(), "results", "script5_state_specific_deg", "GA_included"),
  out_dir_homo  = file.path(getwd(), "results", "script5_state_specific_deg", "homo_only"),
  out_dir_cmpga = file.path(getwd(), "results", "script5_state_specific_deg", "GA_included", "comparison"),
  out_dir_cmpho = file.path(getwd(), "results", "script5_state_specific_deg", "homo_only", "comparison"),
  
  # Column names in colData(pb_MG)
  genotype_num_col      = "rs3732765_num",
  genotype_homo_col     = "rs3732765",
  collapsed_col         = "rs3732765_collapsed",
  pathology_col         = "pathology_group",
  braak_col             = "braak_bins",
  dataset_col           = "dataset",
  age_col               = "age",
  sex_col               = "sex",
  brain_region_col      = "brain_region",
  pmi_col               = "PMI",
  participant_col       = "participant",
  n_genes_col           = "n_genes",
  n_counts_col          = "n_counts",
  percent_mito_col      = "percent_mito",
  mito_ribo_col         = "mito_ribo",
  
  # Model settings
  min_samples_after_filter = 15,
  min_lib_size             = 100,
  min_genes_after_filter   = 100,
  
  # Labels/colors
  state_colors = c(
    "M3: SPP1"  = "#53B400",
    "M5: PCDH9" = "#00A6FF"
  ),
  
  # Key genes to inspect in PM results
  key_genes = c(
    "MIF", "RPL17", "RPS17", "RPS10", "MT2A",
    "NDUFC2", "NDUFB8", "NDUFA13", "MGP", "IFITM1",
    "BLOC1S1", "KRTCAP2", "CRABP1", "SGCZ", "TIMP1", "TIMP3"
  ),
  
  # Scatter-plot tuning
  highlight_fdr_cutoff = 0.05,
  label_all_if_small = TRUE
)

dir.create(config$out_dir, showWarnings = FALSE, recursive = TRUE)
dir.create(config$out_dir_ga, showWarnings = FALSE, recursive = TRUE)
dir.create(config$out_dir_homo, showWarnings = FALSE, recursive = TRUE)
dir.create(config$out_dir_cmpga, showWarnings = FALSE, recursive = TRUE)
dir.create(config$out_dir_cmpho, showWarnings = FALSE, recursive = TRUE)

############################################################
## 1. Packages
############################################################

required_pkgs <- c(
  "SummarizedExperiment",
  "variancePartition",
  "limma",
  "edgeR",
  "dplyr",
  "data.table",
  "ggplot2",
  "ggrepel",
  "patchwork",
  "readr",
  "tidyr",
  "forcats"
)

missing_pkgs <- required_pkgs[!vapply(required_pkgs, requireNamespace, logical(1), quietly = TRUE)]
if (length(missing_pkgs) > 0) {
  stop(
    "Missing packages: ", paste(missing_pkgs, collapse = ", "),
    "\nInstall them first (ideally with renv) before running this script."
  )
}

suppressPackageStartupMessages({
  library(SummarizedExperiment)
  library(variancePartition)
  library(limma)
  library(edgeR)
  library(dplyr)
  library(data.table)
  library(ggplot2)
  library(ggrepel)
  library(patchwork)
  library(readr)
  library(tidyr)
  library(forcats)
})

############################################################
## 2. Helper functions
############################################################

`%||%` <- function(x, y) if (!is.null(x)) x else y

load_object <- function(object_name = NULL, rds_path = NULL, label = "object") {
  if (!is.null(object_name) && exists(object_name, inherits = TRUE)) {
    return(get(object_name, inherits = TRUE))
  }
  if (!is.null(rds_path) && file.exists(rds_path)) {
    return(readRDS(rds_path))
  }
  stop("Could not load ", label, ": provide an object in the environment or an RDS path.")
}

save_both <- function(plot_obj, file_base, out_dir, width = 7, height = 6, dpi = 300) {
  ggsave(file.path(out_dir, paste0(file_base, ".png")), plot_obj, width = width, height = height, dpi = dpi)
  ggsave(file.path(out_dir, paste0(file_base, ".pdf")), plot_obj, width = width, height = height)
}

clean_file_tag <- function(x) {
  x <- gsub("[: ]", "_", x)
  x <- gsub("__+", "_", x)
  x <- gsub("^_+|_+$", "", x)
  x
}

make_sample_info_pm <- function(pb, config, homo_only = FALSE) {
  meta <- as.data.frame(colData(pb))
  meta$sample <- rownames(meta)
  
  needed <- c(
    config$genotype_num_col,
    config$pathology_col,
    config$braak_col,
    config$dataset_col,
    config$age_col,
    config$sex_col,
    config$brain_region_col,
    config$pmi_col,
    config$participant_col,
    config$n_genes_col,
    config$n_counts_col,
    config$percent_mito_col,
    config$mito_ribo_col
  )
  
  missing_cols <- setdiff(needed, colnames(meta))
  if (length(missing_cols) > 0) {
    stop("PM object is missing required metadata columns: ", paste(missing_cols, collapse = ", "))
  }
  
  meta <- meta %>%
    filter(
      !is.na(.data[[config$genotype_num_col]]),
      !is.na(.data[[config$pathology_col]]),
      !is.na(.data[[config$braak_col]]),
      !is.na(.data[[config$dataset_col]]),
      !is.na(.data[[config$age_col]]),
      !is.na(.data[[config$sex_col]]),
      !is.na(.data[[config$brain_region_col]]),
      !is.na(.data[[config$pmi_col]]),
      !is.na(.data[[config$participant_col]]),
      !is.na(.data[[config$n_genes_col]]),
      !is.na(.data[[config$n_counts_col]]),
      !is.na(.data[[config$percent_mito_col]]),
      !is.na(.data[[config$mito_ribo_col]])
    )
  
  if (homo_only) {
    meta <- meta %>% filter(.data[[config$genotype_homo_col]] %in% c("G/G", "A/A"))
  }
  
  meta <- meta %>%
    mutate(
      rs3732765_num = as.numeric(.data[[config$genotype_num_col]]),
      genotype_binary = factor(.data[[config$genotype_homo_col]], levels = c("G/G", "A/A")),
      pathology_group = factor(.data[[config$pathology_col]]),
      braak_bins      = factor(.data[[config$braak_col]]),
      dataset         = factor(.data[[config$dataset_col]]),
      sex             = factor(.data[[config$sex_col]]),
      brain_region    = factor(.data[[config$brain_region_col]]),
      participant     = factor(.data[[config$participant_col]])
    )
  
  if (homo_only) {
    meta <- meta %>% filter(genotype_binary %in% c("G/G", "A/G", "A/A"))
    meta$genotype_binary <- droplevels(meta$genotype_binary)
  }
  
  meta
}

pm_formula_full <- function(config) {
  as.formula(paste0(
    "~ rs3732765_num + ",
    "(1 | pathology_group) + ",
    "(1 | braak_bins) + ",
    "scale(age) + ",
    "(1 | sex) + ",
    "(1 | brain_region) + ",
    "scale(PMI) + ",
    "(1 | participant) + ",
    "log(n_genes) + ",
    "log(n_counts) + ",
    "scale(percent_mito) + ",
    "scale(mito_ribo)"
  ))
}

pm_formula_homo <- function() {
  ~ genotype_binary + brain_region + (1 | participant)
}

prepare_state_counts <- function(pb, state_name, sample_info) {
  counts_state <- assay(pb, state_name)
  common_samples <- intersect(colnames(counts_state), sample_info$sample)
  if (length(common_samples) == 0) {
    stop("No shared samples between state assay and sample info for ", state_name)
  }
  
  counts_state <- counts_state[, common_samples, drop = FALSE]
  si <- sample_info[match(common_samples, sample_info$sample), , drop = FALSE]
  
  # Remove low-library-size samples for this specific state assay.
  lib_sizes <- colSums(counts_state)
  remove_samples <- names(lib_sizes)[lib_sizes < config$min_lib_size]
  if (length(remove_samples) > 0) {
    counts_state <- counts_state[, !colnames(counts_state) %in% remove_samples, drop = FALSE]
    si <- si %>% filter(sample %in% colnames(counts_state))
  }
  
  list(counts_state = counts_state, sample_info = si)
}

fit_state_dream <- function(counts_state, si, formula, coef_pattern, min_genes_after_filter,
                            out_dir, state_name, file_tag, genotype_label = "rs3732765_num") {
  if (ncol(counts_state) < 2) return(NULL)
  
  si <- si %>% droplevels()
  
  dge <- DGEList(counts = counts_state)
  keep <- filterByExpr(dge, group = si[[1]])
  dge <- dge[keep, , keep.lib.sizes = FALSE]
  if (nrow(dge) < min_genes_after_filter) return(NULL)
  dge <- calcNormFactors(dge)
  
  v <- tryCatch(
    voomWithDreamWeights(dge, formula, si, BPPARAM = SerialParam()),
    error = function(e) {
      message("voomWithDreamWeights failed for ", state_name, ": ", conditionMessage(e))
      NULL
    }
  )
  if (is.null(v)) return(NULL)
  
  fit <- tryCatch(
    dream(v, formula, si, BPPARAM = SerialParam()),
    error = function(e) {
      message("dream failed for ", state_name, ": ", conditionMessage(e))
      NULL
    }
  )
  if (is.null(fit)) return(NULL)
  
  fit <- eBayes(fit)
  coef_name <- grep(coef_pattern, colnames(coef(fit)), value = TRUE)
  coef_name <- setdiff(coef_name, "(Intercept)")
  if (length(coef_name) != 1) {
    message("Could not uniquely identify coefficient for ", state_name, ": ", paste(coef_name, collapse = ", "))
    return(NULL)
  }
  
  res <- topTable(fit, coef = coef_name, number = Inf, sort.by = "P")
  res$symbol <- rownames(res)
  
  out_file <- file.path(out_dir, paste0(file_tag, "_", clean_file_tag(state_name), ".tsv"))
  fwrite(res, out_file, sep = "\t")
  
  list(
    state_name = state_name,
    counts_state = counts_state,
    sample_info = si,
    dge = dge,
    v = v,
    fit = fit,
    res = res,
    coef_name = coef_name,
    out_file = out_file
  )
}

state_scatter_plot <- function(pm_res, mb_res, state_label, state_color,
                               out_dir, fname_suffix = "",
                               pm_model_title = "PM model") {
  shared <- inner_join(
    mb_res %>% select(symbol, logFC, P.Value, adj.P.Val) %>%
      rename(logFC_miBrain = logFC, pval_miBrain = P.Value, fdr_miBrain = adj.P.Val),
    pm_res %>% select(symbol, logFC, P.Value, adj.P.Val) %>%
      rename(logFC_PM = logFC, pval_PM = P.Value, fdr_PM = adj.P.Val),
    by = "symbol"
  )
  
  r_pearson  <- cor(shared$logFC_PM, shared$logFC_miBrain, method = "pearson")
  r_spearman <- cor(shared$logFC_PM, shared$logFC_miBrain, method = "spearman")
  
  shared <- shared %>%
    mutate(
      concordant = sign(logFC_PM) == sign(logFC_miBrain),
      pm_sig_tier = case_when(
        fdr_PM < 0.05 ~ "FDR",
        pval_PM < 0.05 ~ "Nominal",
        TRUE          ~ "NS"
      ),
      point_cat = case_when(
        pm_sig_tier == "FDR"     & concordant  ~ "Concordant, sig in PM",
        pm_sig_tier == "FDR"     & !concordant ~ "Discordant, sig in PM",
        pm_sig_tier == "Nominal"  & concordant  ~ "Concordant, sig in PM",
        pm_sig_tier == "Nominal"  & !concordant ~ "Discordant, sig in PM",
        pm_sig_tier == "NS"       & concordant  ~ "Concordant, NS in PM",
        pm_sig_tier == "NS"       & !concordant ~ "Discordant, NS in PM",
        TRUE ~ "Background"
      )
    )
  
  key_df <- shared %>% filter(pm_sig_tier != "NS")
  conc_pct <- if (nrow(key_df) > 0) round(mean(key_df$concordant) * 100, 1) else NA_real_
  
  label_df <- key_df %>%
    arrange(pval_PM) %>%
    slice_head(n = 30)
  
  p <- ggplot(shared, aes(x = logFC_PM, y = logFC_miBrain)) +
    geom_hline(yintercept = 0, linetype = "dashed", color = "grey60", linewidth = 0.4) +
    geom_vline(xintercept = 0, linetype = "dashed", color = "grey60", linewidth = 0.4) +
    geom_point(data = shared %>% filter(point_cat == "Background"), color = "grey82", alpha = 0.20, size = 0.5) +
    geom_point(data = shared %>% filter(point_cat == "Concordant, NS in PM"), color = "#2166AC", alpha = 0.55, size = 2.5, shape = 21, fill = "white", stroke = 0.8) +
    geom_point(data = shared %>% filter(point_cat == "Discordant, NS in PM"), color = "#D6604D", alpha = 0.55, size = 2.5, shape = 21, fill = "white", stroke = 0.8) +
    geom_point(data = shared %>% filter(point_cat == "Concordant, sig in PM"), color = "#2166AC", alpha = 0.95, size = 3.0) +
    geom_point(data = shared %>% filter(point_cat == "Discordant, sig in PM"), color = "#D6604D", alpha = 0.95, size = 3.0) +
    geom_smooth(method = "lm", se = TRUE, color = "grey30", linewidth = 0.5) +
    geom_text_repel(
      data = label_df,
      aes(label = symbol),
      size = 3.0,
      fontface = "italic",
      color = "grey20",
      max.overlaps = 40,
      box.padding = 0.35,
      point.padding = 0.25,
      segment.size = 0.3,
      segment.color = "grey50"
    ) +
    annotate("text", x = -Inf, y = Inf,
             label = paste0("Pearson r = ", round(r_pearson, 2),
                            "\nSpearman ρ = ", round(r_spearman, 2),
                            "\nn = ", nrow(shared), " shared genes"),
             hjust = -0.07, vjust = 1.25, size = 3.2, color = "grey20") +
    annotate("text", x = Inf, y = Inf,
             label = paste0(ifelse(is.na(conc_pct), "NA", conc_pct), "% directionally concordant\n",
                            "Filled = sig in PM | Open = NS in PM"),
             hjust = 1.07, vjust = 1.25, size = 3.2, color = "#2166AC") +
    annotate("label", x = 0, y = -Inf,
             label = paste0(state_label, " state"),
             fill = state_color, color = "white", fontface = "bold",
             size = 3.8, label.padding = unit(0.35, "lines"),
             label.r = unit(0.25, "lines"), vjust = -0.3) +
    labs(
      title = paste0(state_label, "-state DEG concordance: PM vs miBrain Pilot"),
      subtitle = paste0(
        pm_model_title, " | ",
        "miBrain Pilot DREAM AA vs GG (SC + NUC)\n",
        "Blue = concordant | Red = discordant | Filled = sig in PM | Open = NS in PM"
      ),
      x = "log2FC — Postmortem",
      y = "log2FC — miBrain Pilot (A/A vs G/G)"
    ) +
    theme_classic(base_size = 12) +
    theme(plot.title = element_text(face = "bold", size = 13),
          plot.subtitle = element_text(size = 8, color = "grey40"),
          legend.position = "none")
  
  save_both(p, file.path(out_dir, paste0("StateScatter_", clean_file_tag(state_label), fname_suffix)), out_dir, width = 7, height = 6)
  p
}

############################################################
## 3. Load PM object
############################################################

pm <- load_object(config$pm_object_name, config$pm_rds_path, label = "PM pseudobulk object")
if (!inherits(pm, "SummarizedExperiment")) {
  stop("PM object must be a SummarizedExperiment-like object.")
}

message("Loaded PM object: ", paste(class(pm), collapse = ", "))
message("Assays found: ", paste(assayNames(pm), collapse = ", "))
message("Assay/class check for first state:")
if (config$states_to_run[1] %in% assayNames(pm)) {
  x <- assay(pm, config$states_to_run[1])
  message("  ", config$states_to_run[1], ": class=", class(x)[1],
          ", dim=", paste(dim(x), collapse = " x "),
          ", range=", paste(range(x, na.rm = TRUE), collapse = " to "))
}

############################################################
## 4. Run GA-included PM state DEG
############################################################

sample_info_ga <- make_sample_info_pm(pm, config, homo_only = FALSE)
form_ga <- pm_formula_full(config)
message("PM additive formula: ", deparse(form_ga))

pm_results_ga <- list()
for (st in config$states_to_run) {
  message("\n=== Running GA-included PM state DEG: ", st, " ===")
  prep <- prepare_state_counts(pm, st, sample_info_ga)
  fit_out <- fit_state_dream(
    counts_state = prep$counts_state,
    si = prep$sample_info,
    formula = form_ga,
    coef_pattern = "rs3732765_num",
    min_genes_after_filter = config$min_genes_after_filter,
    out_dir = config$out_dir_ga,
    state_name = st,
    file_tag = "PM_DREAM_GA_included"
  )
  pm_results_ga[[st]] <- fit_out
  
  if (!is.null(fit_out)) {
    message("Top genes for ", st, ":")
    print(head(fit_out$res[, c("symbol", "logFC", "P.Value", "adj.P.Val")], 15))
    
    # Key genes lookup
    hits <- fit_out$res %>% filter(symbol %in% config$key_genes) %>% arrange(P.Value)
    if (nrow(hits) > 0) {
      message("Key genes found in ", st, ": ")
      print(hits %>% select(symbol, logFC, P.Value, adj.P.Val))
    }
  }
}

############################################################
## 5. Run homozygote-only PM state DEG
############################################################

sample_info_homo <- make_sample_info_pm(pm, config, homo_only = TRUE) %>%
  filter(genotype_binary %in% c("G/G", "A/A"))
form_homo <- pm_formula_homo()
message("PM homozygote formula: ", deparse(form_homo))

pm_results_homo <- list()
for (st in config$states_to_run) {
  message("\n=== Running homozygote-only PM state DEG: ", st, " ===")
  prep <- prepare_state_counts(pm, st, sample_info_homo)
  fit_out <- fit_state_dream(
    counts_state = prep$counts_state,
    si = prep$sample_info,
    formula = form_homo,
    coef_pattern = "genotype_binaryA",
    min_genes_after_filter = config$min_genes_after_filter,
    out_dir = config$out_dir_homo,
    state_name = st,
    file_tag = "PM_DREAM_homo_only"
  )
  pm_results_homo[[st]] <- fit_out
  
  if (!is.null(fit_out)) {
    message("Top genes for ", st, " (homo-only):")
    print(head(fit_out$res[, c("symbol", "logFC", "P.Value", "adj.P.Val")], 15))
    
    hits <- fit_out$res %>% filter(symbol %in% config$key_genes) %>% arrange(P.Value)
    if (nrow(hits) > 0) {
      message("Key genes found in ", st, " (homo-only):")
      print(hits %>% select(symbol, logFC, P.Value, adj.P.Val))
    }
  }
}

############################################################
## 6. Prepare helper for state-specific comparison tables
############################################################

read_mb_state_results <- function(path, label) {
  if (is.null(path) || !file.exists(path)) {
    message("No miBrain file for ", label, "; comparison skipped for that state.")
    return(NULL)
  }
  read_tsv(path, show_col_types = FALSE) %>%
    rename(
      symbol = any_of("symbol"),
      logFC_miBrain = any_of("logFC"),
      pval_miBrain = any_of("P.Value"),
      fdr_miBrain = any_of("adj.P.Val")
    )
}

make_state_comparison <- function(pm_fit_out, mb_df, state_label, out_dir, mode_label) {
  if (is.null(pm_fit_out) || is.null(mb_df)) return(NULL)
  
  pm_df <- pm_fit_out$res %>%
    select(symbol, logFC_PM = logFC, pval_PM = P.Value, fdr_PM = adj.P.Val)
  
  shared <- inner_join(
    mb_df %>% select(symbol, logFC_miBrain, pval_miBrain, fdr_miBrain),
    pm_df,
    by = "symbol"
  )
  
  if (nrow(shared) == 0) return(NULL)
  
  shared <- shared %>%
    mutate(
      concordant = sign(logFC_PM) == sign(logFC_miBrain),
      pm_sig_tier = case_when(
        fdr_PM < 0.05 ~ "FDR",
        pval_PM < 0.05 ~ "Nominal",
        TRUE ~ "NS"
      )
    )
  
  r_p <- cor(shared$logFC_PM, shared$logFC_miBrain, method = "pearson")
  r_s <- cor(shared$logFC_PM, shared$logFC_miBrain, method = "spearman")
  
  message("\n=== ", state_label, " comparison (", mode_label, ") ===")
  message("Shared genes: ", nrow(shared))
  message("Pearson r: ", round(r_p, 3))
  message("Spearman rho: ", round(r_s, 3))
  message("miBrain nominal p<0.05: ", sum(shared$pval_miBrain < 0.05, na.rm = TRUE))
  message("Concordant among miBrain nominal p<0.05: ",
          sum(shared$concordant & shared$pval_miBrain < 0.05, na.rm = TRUE),
          " (",
          round(mean(shared$concordant[shared$pval_miBrain < 0.05], na.rm = TRUE) * 100, 1),
          "%)")
  
  # Save table
  out_tbl <- shared %>% arrange(pval_PM)
  write_csv(out_tbl, file.path(out_dir, paste0("StateCompare_", clean_file_tag(state_label), "_", mode_label, ".csv")))
  
  # Plot
  plot_df <- shared %>%
    mutate(
      is_key = symbol %in% config$key_genes,
      label = ifelse(is_key & pval_PM < config$highlight_fdr_cutoff, symbol, NA_character_),
      point_type = case_when(
        pval_PM < 0.05 & concordant ~ "PM sig, concordant",
        pval_PM < 0.05 & !concordant ~ "PM sig, discordant",
        TRUE ~ "Background"
      )
    )
  
  p <- ggplot(plot_df, aes(x = logFC_PM, y = logFC_miBrain)) +
    geom_hline(yintercept = 0, linetype = "dashed", color = "grey60", linewidth = 0.4) +
    geom_vline(xintercept = 0, linetype = "dashed", color = "grey60", linewidth = 0.4) +
    geom_point(data = filter(plot_df, point_type == "Background"), color = "grey80", alpha = 0.25, size = 0.5) +
    geom_point(data = filter(plot_df, point_type == "PM sig, discordant"), color = "#D6604D", alpha = 0.85, size = 2.0) +
    geom_point(data = filter(plot_df, point_type == "PM sig, concordant"), color = "#2166AC", alpha = 0.90, size = 2.2) +
    geom_smooth(method = "lm", se = TRUE, color = "grey30", linewidth = 0.5) +
    geom_text_repel(
      data = filter(plot_df, !is.na(label)),
      aes(label = label),
      size = 2.8,
      color = "#2166AC",
      max.overlaps = 40,
      box.padding = 0.35,
      point.padding = 0.2,
      segment.color = "grey50",
      segment.size = 0.3
    ) +
    annotate("text", x = -Inf, y = Inf,
             label = paste0("Pearson r = ", round(r_p, 2),
                            "\nSpearman ρ = ", round(r_s, 2),
                            "\nn = ", nrow(shared), " shared genes"),
             hjust = -0.07, vjust = 1.25, size = 3.2, color = "grey20") +
    annotate("text", x = Inf, y = Inf,
             label = paste0(round(mean(shared$concordant & shared$pval_miBrain < 0.05, na.rm = TRUE) * 100, 1),
                            "% concordant among PM sig genes\n",
                            "Filled = PM sig | Blue = concordant"),
             hjust = 1.07, vjust = 1.25, size = 3.2, color = "#2166AC") +
    annotate("label", x = 0, y = -Inf,
             label = paste0(state_label, " state"),
             fill = config$state_colors[[state_label]], color = "white",
             fontface = "bold", size = 3.8,
             label.padding = unit(0.35, "lines"), label.r = unit(0.25, "lines"),
             vjust = -0.3) +
    labs(
      title = paste0(state_label, "-state DEG concordance: PM vs miBrain Pilot"),
      subtitle = paste0("Mode: ", mode_label, " | filled = PM sig | open = direction only"),
      x = "log2FC — Postmortem",
      y = "log2FC — miBrain Pilot (A/A vs G/G)"
    ) +
    theme_classic(base_size = 12) +
    theme(plot.title = element_text(face = "bold", size = 13),
          plot.subtitle = element_text(size = 8, color = "grey40"),
          legend.position = "none")
  
  save_both(
    p,
    file.path(out_dir, paste0("StateScatter_", clean_file_tag(state_label), "_", mode_label)),
    out_dir,
    width = 7,
    height = 6
  )
  
  list(table = out_tbl, plot = p, shared = shared)
}

############################################################
## 7. Load miBrain comparison tables
############################################################

mb_tables <- lapply(config$miBrain_state_paths, read_mb_state_results)

############################################################
## 8. Build comparison plots for each state and model mode
############################################################

comparison_outputs <- list()
for (st in config$states_to_run) {
  mb_df <- mb_tables[[st]]
  if (is.null(mb_df)) next
  
  # GA-included comparisons if PM results exist.
  if (!is.null(pm_results_ga[[st]])) {
    comparison_outputs[[paste0(st, "_GA")]] <- make_state_comparison(
      pm_fit_out = pm_results_ga[[st]],
      mb_df = mb_df,
      state_label = st,
      out_dir = config$out_dir_cmpga,
      mode_label = "GA_included"
    )
  }
  
  # Homo-only comparisons if PM results exist.
  if (!is.null(pm_results_homo[[st]])) {
    comparison_outputs[[paste0(st, "_HOMO")]] <- make_state_comparison(
      pm_fit_out = pm_results_homo[[st]],
      mb_df = mb_df,
      state_label = st,
      out_dir = config$out_dir_cmpho,
      mode_label = "homo_only"
    )
  }
}

############################################################
## 9. Summary tables
############################################################

# Collate PM results into lightweight summaries for easy reading.
make_summary_from_results <- function(res_list, mode_name, out_dir) {
  summary_tbl <- bind_rows(lapply(names(res_list), function(st) {
    x <- res_list[[st]]
    if (is.null(x)) return(NULL)
    x$res %>%
      mutate(state = st, mode = mode_name) %>%
      select(state, mode, symbol, logFC, P.Value, adj.P.Val)
  }))
  
  if (nrow(summary_tbl) > 0) {
    write_csv(summary_tbl, file.path(out_dir, paste0("PM_state_DEG_summary_", mode_name, ".csv")))
  }
  summary_tbl
}

summary_ga   <- make_summary_from_results(pm_results_ga,   "GA_included", config$out_dir_ga)
summary_homo <- make_summary_from_results(pm_results_homo, "homo_only",  config$out_dir_homo)

############################################################
## 10. Final messages
############################################################

message("\nScript 5 complete.")
message("Outputs written to: ", normalizePath(config$out_dir))
message("  - PM state DEG results for GA-included and homo-only models")
message("  - optional PM vs miBrain state comparison plots/tables")
message("  - PM summary tables by mode")
