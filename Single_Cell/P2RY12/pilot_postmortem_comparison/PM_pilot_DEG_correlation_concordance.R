############################################################
## SCRIPT 4
## DEG concordance + pathway analysis
##
## Purpose:
##   - read PM DEG table
##   - read pilot DEG table
##   - merge on symbol
##   - compute beta correlations
##   - make scatter plots
##   - rank concordant genes
##   - run fgsea / pathway enrichment on the ranked list
##   - save pathway plots and tables
##
## This script is separate from abundance analysis.
## It does NOT:
##   - compute state proportions
##   - run limma on abundance tables
##   - rebuild pseudobulk
############################################################

############################################################
## 0. Configuration
############################################################

config <- list(
  out_dir = file.path(getwd(), "results", "script4_deg_concordance_pathway"),
  
  # Saved DEG tables from Script 1 / PM DEG pipeline.
  pm_deg_path    = NULL,
  pilot_deg_path  = NULL,
  
  # Column names in the DEG files
  symbol_col = "symbol",
  beta_col   = "logFC",
  p_col      = "P.Value",
  fdr_col    = "adj.P.Val",
  
  # Dataset labels for plots
  pm_label    = "Postmortem",
  pilot_label = "miBrain Pilot",
  
  # Ranking / filtering
  abs_beta_filter = 0.5,
  top_n_labels = 20,
  label_all_concordant = FALSE,
  
  # Gene sets to label in scatter plots if present
  highlight_genes = c(
    "P2RY12", "TREM2", "APOE", "SPP1", "PCDH9",
    "LAMP1", "LAMP2", "ATG7", "MAP1LC3B", "SQSTM1",
    "RAC1", "RHEB", "MTOR", "AKT1", "AKT2", "CTSD"
  ),
  
  # Pathway analysis settings
  fgsea_min_size = 15,
  fgsea_max_size = 500,
  fgsea_nperm    = 10000,
  fgsea_seed     = 42,
  
  # Gene-set collections to run
  run_hallmark = TRUE,
  run_reactome = TRUE,
  run_kegg     = TRUE,
  run_gobp     = TRUE,
  run_gomf     = TRUE,
  run_gocc     = TRUE,
  
  # File names
  scatter_all_name   = "DEG_beta_concordance_ALL",
  scatter_filt_name  = "DEG_beta_concordance_FILTERED",
  ranked_all_name    = "DEG_beta_concordance_ALL_genes_ranked.csv",
  ranked_conc_name   = "concordant_genes_ranked.csv",
  ranked_disc_name   = "discordant_genes_ranked.csv",
  cor_summary_name   = "DEG_beta_correlation_summary.csv",
  fgsea_dir_name     = "fgsea_results",
  fgsea_plot_dir     = "fgsea_plots"
)

dir.create(config$out_dir, showWarnings = FALSE, recursive = TRUE)
dir.create(file.path(config$out_dir, config$fgsea_dir_name), showWarnings = FALSE, recursive = TRUE)
dir.create(file.path(config$out_dir, config$fgsea_plot_dir), showWarnings = FALSE, recursive = TRUE)

############################################################
## 1. Packages
############################################################

required_pkgs <- c(
  "dplyr",
  "readr",
  "ggplot2",
  "ggrepel",
  "forcats",
  "tidyr",
  "msigdbr",
  "fgsea",
  "stringr",
  "tibble"
)

missing_pkgs <- required_pkgs[!vapply(required_pkgs, requireNamespace, logical(1), quietly = TRUE)]
if (length(missing_pkgs) > 0) {
  stop(
    "Missing packages: ", paste(missing_pkgs, collapse = ", "),
    "\nInstall them first (ideally with renv) before running this script."
  )
}

suppressPackageStartupMessages({
  library(dplyr)
  library(readr)
  library(ggplot2)
  library(ggrepel)
  library(forcats)
  library(tidyr)
  library(msigdbr)
  library(fgsea)
  library(stringr)
  library(tibble)
})

############################################################
## 2. Helper functions
############################################################

`%||%` <- function(x, y) if (!is.null(x)) x else y

load_results_table <- function(path, label = "results table") {
  if (is.null(path) || !file.exists(path)) {
    stop("Could not find ", label, ". Provide a valid path.")
  }
  read_csv(path, show_col_types = FALSE)
}

standardize_deg <- function(df, dataset_label, symbol_col, beta_col, p_col, fdr_col) {
  missing_cols <- setdiff(c(symbol_col, beta_col), colnames(df))
  if (length(missing_cols) > 0) {
    stop(dataset_label, " DEG table is missing required columns: ", paste(missing_cols, collapse = ", "))
  }
  
  out <- df %>%
    rename(
      symbol = all_of(symbol_col),
      beta = all_of(beta_col)
    )
  
  if (p_col %in% colnames(out)) {
    out <- out %>% rename(pval = all_of(p_col))
  } else {
    out$pval <- NA_real_
  }
  
  if (fdr_col %in% colnames(out)) {
    out <- out %>% rename(fdr = all_of(fdr_col))
  } else {
    out$fdr <- NA_real_
  }
  
  out %>%
    mutate(
      symbol = as.character(symbol),
      dataset = dataset_label
    )
}

save_both <- function(plot_obj, file_base, out_dir, width = 6, height = 6, dpi = 300) {
  ggsave(file.path(out_dir, paste0(file_base, ".png")), plot_obj, width = width, height = height, dpi = dpi)
  ggsave(file.path(out_dir, paste0(file_base, ".pdf")), plot_obj, width = width, height = height)
}

compute_correlations <- function(df) {
  d <- df %>% filter(!is.na(beta_pm), !is.na(beta_pilot))
  if (nrow(d) < 2) {
    return(list(n = nrow(d), pearson = NA_real_, spearman = NA_real_))
  }
  list(
    n = nrow(d),
    pearson = cor(d$beta_pm, d$beta_pilot, method = "pearson"),
    spearman = cor(d$beta_pm, d$beta_pilot, method = "spearman")
  )
}

make_deg_scatter <- function(df,
                             title,
                             subtitle,
                             highlight_genes = NULL,
                             label_all = FALSE,
                             top_n = 20,
                             beta_xlim = NULL,
                             beta_ylim = NULL) {
  plot_df <- df %>%
    mutate(
      concordant = sign(beta_pm) == sign(beta_pilot),
      point_type = case_when(
        concordant ~ "Concordant",
        !concordant ~ "Discordant",
        TRUE ~ "Background"
      )
    )
  
  label_df <- if (label_all) {
    plot_df
  } else {
    plot_df %>% arrange(desc(abs(beta_pilot))) %>% slice_head(n = top_n)
  }
  
  if (!is.null(highlight_genes)) {
    extra <- plot_df %>% filter(symbol %in% highlight_genes)
    label_df <- bind_rows(label_df, extra) %>% distinct(symbol, .keep_all = TRUE)
  }
  
  corr <- compute_correlations(plot_df)
  cor_label <- paste0(
    "Pearson r = ", round(corr$pearson, 2),
    "\nSpearman ", "ρ = ", round(corr$spearman, 2),
    "\nn = ", corr$n
  )
  
  concordance_pct <- round(mean(plot_df$concordant, na.rm = TRUE) * 100, 1)
  
  p <- ggplot(plot_df, aes(x = beta_pm, y = beta_pilot)) +
    geom_hline(yintercept = 0, linetype = "dashed", color = "grey70", linewidth = 0.4) +
    geom_vline(xintercept = 0, linetype = "dashed", color = "grey70", linewidth = 0.4) +
    geom_point(aes(color = concordant), alpha = 0.35, size = 0.8) +
    geom_smooth(method = "lm", se = TRUE, color = "black", linewidth = 0.6) +
    geom_text_repel(
      data = label_df,
      aes(label = symbol),
      size = 2.8,
      max.overlaps = if (label_all) Inf else top_n,
      box.padding = 0.3,
      segment.color = "grey50",
      segment.size = 0.3
    ) +
    scale_color_manual(
      values = c("TRUE" = "#2166AC", "FALSE" = "#D6604D"),
      labels = c("TRUE" = "Concordant", "FALSE" = "Discordant"),
      name = NULL
    ) +
    annotate(
      "text",
      x = -Inf, y = Inf,
      label = cor_label,
      hjust = -0.1, vjust = 1.2,
      size = 3.8
    ) +
    annotate(
      "text",
      x = Inf, y = -Inf,
      label = paste0(concordance_pct, "% concordant direction"),
      hjust = 1.1, vjust = -0.5,
      size = 3.5, color = "grey40"
    ) +
    labs(
      x = paste0("Beta — ", config$pm_label, " (A/A vs G/G)"),
      y = paste0("Beta — ", config$pilot_label, " (A/A vs G/G)"),
      title = title,
      subtitle = subtitle
    ) +
    theme_classic(base_size = 12) +
    theme(
      plot.title = element_text(face = "bold"),
      plot.subtitle = element_text(size = 9, color = "grey40"),
      legend.position = "top"
    )
  
  if (!is.null(beta_xlim)) p <- p + coord_cartesian(xlim = beta_xlim)
  if (!is.null(beta_ylim)) p <- p + coord_cartesian(ylim = beta_ylim)
  
  p
}

make_lollipop_plot <- function(df,
                               title,
                               subtitle,
                               top_n = 40,
                               label_all = FALSE) {
  # Keep one row per gene and display both datasets side-by-side.
  plot_df <- df %>%
    mutate(concordant = sign(beta_pm) == sign(beta_pilot))
  
  label_df <- if (label_all) {
    plot_df
  } else {
    plot_df %>%
      arrange(desc(abs(beta_pilot))) %>%
      slice_head(n = top_n)
  }
  
  long_df <- plot_df %>%
    select(symbol, beta_pm, beta_pilot, concordant) %>%
    pivot_longer(
      cols = c(beta_pm, beta_pilot),
      names_to = "dataset",
      values_to = "beta"
    ) %>%
    mutate(
      dataset = recode(dataset,
                       beta_pm = config$pm_label,
                       beta_pilot = config$pilot_label),
      direction = ifelse(beta > 0, "A/A > G/G", "G/G > A/A")
    )
  
  label_genes <- unique(label_df$symbol)
  
  ggplot(long_df, aes(x = beta, y = fct_reorder(symbol, beta, .fun = mean), color = direction, shape = dataset)) +
    geom_vline(xintercept = 0, linetype = "dashed", color = "grey60") +
    geom_line(aes(group = symbol), color = "grey70", linewidth = 0.4) +
    geom_point(size = 2.5, alpha = 0.9) +
    geom_text_repel(
      data = long_df %>% filter(symbol %in% label_genes) %>% distinct(symbol, .keep_all = TRUE),
      aes(label = symbol),
      size = 2.5,
      max.overlaps = if (label_all) Inf else top_n,
      box.padding = 0.3,
      segment.color = "grey50",
      segment.size = 0.3,
      show.legend = FALSE
    ) +
    scale_color_manual(
      values = c("A/A > G/G" = "#2166AC", "G/G > A/A" = "#D6604D"),
      name = "Direction"
    ) +
    scale_shape_manual(
      values = c(
        !!config$pm_label := 16,
        !!config$pilot_label := 17
      ),
      name = "Dataset"
    ) +
    labs(
      title = title,
      subtitle = subtitle,
      x = "Beta (log2FC, A/A vs G/G)",
      y = NULL
    ) +
    theme_classic(base_size = 11) +
    theme(
      plot.title = element_text(face = "bold"),
      plot.subtitle = element_text(size = 8, color = "grey40"),
      axis.text.y = element_text(size = 8)
    )
}

make_ranked_concordant_table <- function(df) {
  df %>%
    mutate(
      concordant = sign(beta_pm) == sign(beta_pilot),
      abs_beta_pm = abs(beta_pm),
      abs_beta_pilot = abs(beta_pilot),
      direction = case_when(
        beta_pm > 0 & beta_pilot > 0 ~ "Both A/A > G/G",
        beta_pm < 0 & beta_pilot < 0 ~ "Both G/G > A/A",
        beta_pm > 0 & beta_pilot < 0 ~ "PM A/A>G/G, Pilot G/G>A/A",
        beta_pm < 0 & beta_pilot > 0 ~ "PM G/G>A/A, Pilot A/A>G/G",
        TRUE ~ "On axis"
      )
    )
}

load_pathways <- function() {
  list(
    Hallmark = msigdbr(species = "Homo sapiens", category = "H") %>% split(x = .$gene_symbol, f = .$gs_name),
    Reactome = msigdbr(species = "Homo sapiens", category = "C2", subcategory = "CP:REACTOME") %>% split(x = .$gene_symbol, f = .$gs_name),
    KEGG = msigdbr(species = "Homo sapiens", category = "C2", subcategory = "CP:KEGG_LEGACY") %>% split(x = .$gene_symbol, f = .$gs_name),
    GOBP = msigdbr(species = "Homo sapiens", category = "C5", subcategory = "GO:BP") %>% split(x = .$gene_symbol, f = .$gs_name),
    GOMF = msigdbr(species = "Homo sapiens", category = "C5", subcategory = "GO:MF") %>% split(x = .$gene_symbol, f = .$gs_name),
    GOCC = msigdbr(species = "Homo sapiens", category = "C5", subcategory = "GO:CC") %>% split(x = .$gene_symbol, f = .$gs_name)
  )
}

run_fgsea_collection <- function(pathways, ranked_stats, collection_name, out_dir,
                                 min_size = 15, max_size = 500, nperm = 10000,
                                 fdr_cutoff = 0.05, top_n = 20) {
  if (length(pathways) == 0) return(NULL)
  
  set.seed(config$fgsea_seed)
  fg <- fgsea(
    pathways = pathways,
    stats = ranked_stats,
    minSize = min_size,
    maxSize = max_size,
    nPermSimple = nperm
  ) %>% arrange(padj)
  
  write_csv(fg, file.path(out_dir, paste0("fgsea_", collection_name, "_concordant.csv")))
  
  plot_df <- fg %>%
    filter(padj < fdr_cutoff) %>%
    arrange(padj) %>%
    slice_head(n = top_n) %>%
    mutate(
      pathway_clean = gsub(paste0("^", collection_name, "_"), "", pathway),
      pathway_clean = gsub("_", " ", pathway_clean),
      pathway_clean = str_to_title(pathway_clean),
      pathway_clean = str_wrap(pathway_clean, width = 45),
      direction = ifelse(NES > 0, "A/A > G/G", "G/G > A/A"),
      sig_label = case_when(
        padj < 0.001 ~ "***",
        padj < 0.01  ~ "**",
        padj < 0.05  ~ "*",
        TRUE         ~ ""
      )
    )
  
  if (nrow(plot_df) == 0) {
    message("No significant ", collection_name, " pathways at FDR < ", fdr_cutoff)
    return(fg)
  }
  
  p <- ggplot(plot_df, aes(x = NES, y = fct_reorder(pathway_clean, NES), fill = direction)) +
    geom_col(width = 0.72, alpha = 0.9) +
    geom_vline(xintercept = 0, linewidth = 0.5, color = "grey30") +
    geom_text(
      aes(label = sig_label,
          x = ifelse(NES < 0, NES - 0.04, NES + 0.04)),
      hjust = ifelse(plot_df$NES < 0, 1, 0),
      size = 3.2,
      color = "grey20"
    ) +
    scale_fill_manual(
      values = c("A/A > G/G" = "#2166AC", "G/G > A/G" = "#D6604D", "G/G > A/A" = "#D6604D"),
      name = "Higher in"
    ) +
    scale_x_continuous(expand = expansion(mult = c(0.18, 0.18))) +
    labs(
      title = paste0(collection_name, " pathway enrichment"),
      subtitle = paste0(
        "Concordant genes across ", config$pm_label, " and ", config$pilot_label,
        " | FDR < ", fdr_cutoff, " | top ", top_n, " shown\n",
        "NES > 0: higher in A/A (protective) | NES < 0: higher in G/G (risk)"
      ),
      x = "Normalized Enrichment Score (NES)",
      y = NULL
    ) +
    theme_classic(base_size = 11) +
    theme(
      plot.title = element_text(face = "bold", size = 13),
      plot.subtitle = element_text(size = 8, color = "grey40"),
      axis.text.y = element_text(size = 9),
      legend.position = "top"
    )
  
  h <- max(4, 0.38 * nrow(plot_df) + 2.5)
  save_both(p, file.path(config$fgsea_plot_dir, paste0("fgsea_", collection_name, "_barplot")), out_dir, width = 9, height = h)
  
  fg
}

############################################################
## 3. Load DEG tables
############################################################

pm_deg_raw <- load_results_table(config$pm_deg_path, label = "PM DEG table")
pilot_deg_raw <- load_results_table(config$pilot_deg_path, label = "pilot DEG table")

pm_deg <- standardize_deg(
  pm_deg_raw,
  dataset_label = config$pm_label,
  symbol_col = config$symbol_col,
  beta_col = config$beta_col,
  p_col = config$p_col,
  fdr_col = config$fdr_col
) %>%
  select(symbol, beta_pm = beta, pval_pm = pval, fdr_pm = fdr)

pilot_deg <- standardize_deg(
  pilot_deg_raw,
  dataset_label = config$pilot_label,
  symbol_col = config$symbol_col,
  beta_col = config$beta_col,
  p_col = config$p_col,
  fdr_col = config$fdr_col
) %>%
  select(symbol, beta_pilot = beta, pval_pilot = pval, fdr_pilot = fdr)

message("PM DEG rows: ", nrow(pm_deg))
message("Pilot DEG rows: ", nrow(pilot_deg))

############################################################
## 4. Merge on shared genes
############################################################

shared <- inner_join(pm_deg, pilot_deg, by = "symbol") %>%
  filter(!is.na(beta_pm), !is.na(beta_pilot))

message("Shared genes: ", nrow(shared))

############################################################
## 5. Correlations
############################################################

r_all <- compute_correlations(shared)
shared_filt <- shared %>% filter(abs(beta_pilot) >= config$abs_beta_filter)
r_filt <- compute_correlations(shared_filt)

cor_summary <- tibble(
  comparison = c("All shared genes", paste0("|beta_pilot| >= ", config$abs_beta_filter)),
  n_genes = c(r_all$n, r_filt$n),
  pearson_r = c(r_all$pearson, r_filt$pearson),
  spearman_rho = c(r_all$spearman, r_filt$spearman)
)

write_csv(cor_summary, file.path(config$out_dir, config$cor_summary_name))
message("Saved correlation summary table.")
print(cor_summary)

############################################################
## 6. Scatter plots
############################################################

scatter_all_title <- "DEG beta correlation: PM vs miBrain Pilot"
scatter_all_subtitle <- paste0(
  "All shared genes | PM and pilot effect sizes compared across shared symbols"
)

p_scatter_all <- make_deg_scatter(
  df = shared,
  title = scatter_all_title,
  subtitle = scatter_all_subtitle,
  highlight_genes = config$highlight_genes,
  label_all = config$label_all_concordant,
  top_n = config$top_n_labels
)

save_both(p_scatter_all, config$scatter_all_name, config$out_dir, width = 6.5, height = 6.0)

scatter_filt_title <- paste0("DEG beta correlation: PM vs miBrain Pilot (|β| ≥ ", config$abs_beta_filter, ")")
scatter_filt_subtitle <- paste0(
  "Subset enriched for larger pilot effects | labeled genes include requested pathway genes"
)

p_scatter_filt <- make_deg_scatter(
  df = shared_filt,
  title = scatter_filt_title,
  subtitle = scatter_filt_subtitle,
  highlight_genes = config$highlight_genes,
  label_all = TRUE,
  top_n = config$top_n_labels
)

save_both(p_scatter_filt, config$scatter_filt_name, config$out_dir, width = 6.5, height = 6.0)

############################################################
## 7. Rank concordant and discordant genes
############################################################

shared_ranked <- make_ranked_concordant_table(shared) %>%
  arrange(desc(abs_beta_pilot))

write_csv(shared_ranked, file.path(config$out_dir, config$ranked_all_name))
message("Saved full ranked DEG table.")

concordant_ranked <- shared_ranked %>%
  filter(concordant == TRUE) %>%
  arrange(fdr_pm, desc(abs_beta_pm), desc(abs_beta_pilot))

discordant_ranked <- shared_ranked %>%
  filter(concordant == FALSE) %>%
  arrange(fdr_pm, desc(abs_beta_pm), desc(abs_beta_pilot))

write_csv(concordant_ranked, file.path(config$out_dir, config$ranked_conc_name))
write_csv(discordant_ranked, file.path(config$out_dir, config$ranked_disc_name))

message("Concordant genes: ", nrow(concordant_ranked))
message("Discordant genes: ", nrow(discordant_ranked))

# A compact summary of the strongest concordant genes.
strong_concordant <- concordant_ranked %>%
  select(symbol, beta_pm, beta_pilot, pval_pm, pval_pilot, fdr_pm, fdr_pilot, direction) %>%
  head(30)

write_csv(strong_concordant, file.path(config$out_dir, "top_concordant_genes_30.csv"))

############################################################
## 8. Optional concordant / significant lollipop-style plot
############################################################

# This is a compact visualization of the strongest concordant genes.
# It is not a replacement for the scatter plots, just a readable summary.
if (nrow(concordant_ranked) > 0) {
  conc_plot_df <- concordant_ranked %>%
    slice_head(n = min(40, nrow(concordant_ranked))) %>%
    mutate(symbol = factor(symbol, levels = rev(unique(symbol)))) %>%
    pivot_longer(
      cols = c(beta_pm, beta_pilot),
      names_to = "dataset",
      values_to = "beta"
    ) %>%
    mutate(
      dataset = recode(dataset,
                       beta_pm = config$pm_label,
                       beta_pilot = config$pilot_label),
      direction_color = ifelse(beta > 0, "A/A > G/G", "G/G > A/A")
    )
  
  p_conc_lollipop <- ggplot(conc_plot_df,
                            aes(x = beta, y = symbol, color = direction_color, shape = dataset)) +
    geom_vline(xintercept = 0, linetype = "dashed", color = "grey60") +
    geom_line(aes(group = symbol), color = "grey70", linewidth = 0.4) +
    geom_point(size = 2.5, alpha = 0.9) +
    scale_color_manual(
      values = c("A/A > G/G" = "#2166AC", "G/G > A/A" = "#D6604D"),
      name = "Direction"
    ) +
    scale_shape_manual(
      values = c(config$pm_label = 16, config$pilot_label = 17),
      name = "Dataset"
    ) +
    labs(
      title = "Concordant DEGs: PM vs miBrain Pilot",
      subtitle = paste0(
        "Top concordant genes ranked by PM evidence; each gene shown in both datasets"
      ),
      x = "Beta (log2FC, A/A vs G/G)",
      y = NULL
    ) +
    theme_classic(base_size = 11) +
    theme(
      plot.title = element_text(face = "bold"),
      plot.subtitle = element_text(size = 8, color = "grey40"),
      axis.text.y = element_text(size = 8)
    )
  
  save_both(p_conc_lollipop, "DEG_concordant_lollipop_top40", config$out_dir, width = 7.2, height = 8.5)
}

############################################################
## 9. GSEA / fgsea on concordant ranked genes
############################################################

# Ranked vector: use PM beta, filtered to concordant genes only.
# This is the core input for pathway enrichment.
ranked_concordant <- concordant_ranked %>%
  arrange(desc(beta_pm)) %>%
  select(symbol, beta_pm) %>%
  deframe()

ranked_pm_all <- shared %>%
  arrange(desc(beta_pm)) %>%
  select(symbol, beta_pm) %>%
  deframe()

write_csv(
  tibble(symbol = names(ranked_concordant), beta_pm = unname(ranked_concordant)),
  file.path(config$out_dir, "ranked_concordant_gene_list.csv")
)
write_csv(
  tibble(symbol = names(ranked_pm_all), beta_pm = unname(ranked_pm_all)),
  file.path(config$out_dir, "ranked_pm_all_gene_list.csv")
)

message("Ranked concordant genes: ", length(ranked_concordant))
message("Ranked PM all genes: ", length(ranked_pm_all))

pathways <- load_pathways()
fgsea_results <- list()

if (config$run_hallmark) fgsea_results$Hallmark <- run_fgsea_collection(pathways$Hallmark, ranked_concordant, "HALLMARK", file.path(config$out_dir, config$fgsea_dir_name),
                                                                        min_size = config$fgsea_min_size, max_size = config$fgsea_max_size, nperm = config$fgsea_nperm)
if (config$run_reactome) fgsea_results$Reactome <- run_fgsea_collection(pathways$Reactome, ranked_concordant, "REACTOME", file.path(config$out_dir, config$fgsea_dir_name),
                                                                        min_size = config$fgsea_min_size, max_size = config$fgsea_max_size, nperm = config$fgsea_nperm)
if (config$run_kegg)     fgsea_results$KEGG     <- run_fgsea_collection(pathways$KEGG, ranked_concordant, "KEGG", file.path(config$out_dir, config$fgsea_dir_name),
                                                                        min_size = config$fgsea_min_size, max_size = config$fgsea_max_size, nperm = config$fgsea_nperm)
if (config$run_gobp)     fgsea_results$GOBP     <- run_fgsea_collection(pathways$GOBP, ranked_concordant, "GOBP", file.path(config$out_dir, config$fgsea_dir_name),
                                                                        min_size = config$fgsea_min_size, max_size = config$fgsea_max_size, nperm = config$fgsea_nperm)
if (config$run_gomf)     fgsea_results$GOMF     <- run_fgsea_collection(pathways$GOMF, ranked_concordant, "GOMF", file.path(config$out_dir, config$fgsea_dir_name),
                                                                        min_size = config$fgsea_min_size, max_size = config$fgsea_max_size, nperm = config$fgsea_nperm)
if (config$run_gocc)     fgsea_results$GOCC     <- run_fgsea_collection(pathways$GOCC, ranked_concordant, "GOCC", file.path(config$out_dir, config$fgsea_dir_name),
                                                                        min_size = config$fgsea_min_size, max_size = config$fgsea_max_size, nperm = config$fgsea_nperm)

############################################################
## 10. Optional summary plot across top pathways
############################################################

# Build a simple combined top-pathway table from all collections.
combined_top <- bind_rows(
  lapply(names(fgsea_results), function(nm) {
    fg <- fgsea_results[[nm]]
    if (is.null(fg)) return(NULL)
    fg %>% arrange(padj) %>% slice_head(n = 5) %>% mutate(collection = nm)
  })
) %>%
  filter(!is.na(pathway))

if (nrow(combined_top) > 0) {
  combined_top <- combined_top %>%
    mutate(
      pathway_clean = str_replace_all(pathway, "_", " "),
      pathway_clean = str_to_title(pathway_clean),
      direction = ifelse(NES > 0, "A/A > G/G", "G/G > A/A"),
      neg_log_padj = -log10(padj + 1e-10)
    )
  
  p_summary <- ggplot(combined_top,
                      aes(x = NES, y = fct_reorder(pathway_clean, NES), fill = direction, size = neg_log_padj)) +
    geom_point(shape = 21, stroke = 0.3) +
    geom_vline(xintercept = 0, linetype = "dashed", color = "grey50", linewidth = 0.4) +
    facet_wrap(~ collection, scales = "free_y", ncol = 2) +
    scale_fill_manual(values = c("A/A > G/G" = "#2166AC", "G/G > A/G" = "#D6604D", "G/G > A/G" = "#D6604D", "G/G > A/G" = "#D6604D"), name = "Direction") +
    scale_size_continuous(range = c(3, 9), name = "-log10(FDR)") +
    labs(
      title = "Convergent pathway signal across gene-set collections",
      subtitle = paste0("Top 5 pathways per collection | concordant genes ranked by PM beta"),
      x = "NES",
      y = NULL
    ) +
    theme_classic(base_size = 10) +
    theme(
      plot.title = element_text(face = "bold", size = 12),
      plot.subtitle = element_text(size = 8, color = "grey40"),
      strip.text = element_text(face = "bold", size = 10),
      strip.background = element_rect(fill = "grey95"),
      axis.text.y = element_text(size = 8),
      legend.position = "right"
    )
  
  save_both(p_summary, "fgsea_summary_4collections", config$out_dir, width = 11, height = 9)
}

############################################################
## 11. Final outputs
############################################################

message("\nScript 4 complete.")
message("Outputs written to: ", normalizePath(config$out_dir))
message("  - DEG correlation summary")
message("  - DEG concordance scatter plots")
message("  - ranked concordant/discordant tables")
message("  - fgsea tables")
message("  - fgsea pathway plots")
message("  - optional concordant lollipop summary")
