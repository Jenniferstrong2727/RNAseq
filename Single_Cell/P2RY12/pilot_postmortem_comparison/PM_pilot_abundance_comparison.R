############################################################
## SCRIPT 3
## PM vs pilot abundance comparison
##
## Purpose:
##   - read saved PM abundance results
##   - read saved pilot abundance results
##   - merge on module/state
##   - make concordance / comparison plots
##   - save comparison tables
##
## This script is a pure reader + plotter.
## It does NOT:
##   - rebuild state proportions
##   - run limma
##   - run edgeR / DREAM
##   - rebuild pseudobulk matrices
############################################################

############################################################
## 0. Configuration
############################################################

config <- list(
  out_dir = file.path(getwd(), "results", "script3_pm_vs_pilot_comparison"),
  
  # Saved abundance result tables from Script 2 and the pilot pipeline.
  pm_results_path    = NULL,
  pilot_results_path = NULL,
  
  # Column names in the saved abundance tables
  module_col = "module",
  delta_col  = "delta_percent",
  p_col      = "P.Value",
  fdr_col    = "adj.P.Val",
  logfc_col  = "logFC",
  
  # Plot labels
  pm_label    = "Postmortem",
  pilot_label = "miBrain Pilot",
  
  # Optional custom titles
  comparison_title = "PM vs miBrain Pilot: differential abundance of scVI states",
  
  # Ordering / colors for the state modules
  state_levels = c(
    "M1: P2RY12",
    "M2: HSP90AA1",
    "M3: SPP1",
    "M4: TMEM163",
    "M5: PCDH9",
    "M6: CD163",
    "M7: FRMD4A",
    "M8: B2M"
  ),
  
  state_colors = c(
    "M1: P2RY12"   = "#F8766D",
    "M2: HSP90AA1" = "#C49A00",
    "M3: SPP1"     = "#53B400",
    "M4: TMEM163"  = "#00C1A2",
    "M5: PCDH9"    = "#00A6FF",
    "M6: CD163"    = "#A58AFF",
    "M7: FRMD4A"   = "#FB61D7",
    "M8: B2M"      = "#FF6F91"
  ),
  
  # Plot tuning
  top_n_labels = 20,
  label_all_sig = FALSE,
  sig_fdr_cutoff = 0.05,
  beta_filter_abs = 0.5,
  
  # Output basenames
  scatter_all_name     = "PM_vs_Pilot_abundance_concordance_scatter",
  scatter_filt_name    = "PM_vs_Pilot_abundance_concordance_scatter_filtered",
  lollipop_name        = "PM_vs_Pilot_abundance_lollipop",
  comparison_table_name = "PM_vs_Pilot_abundance_comparison_table.csv"
)

dir.create(config$out_dir, showWarnings = FALSE, recursive = TRUE)

############################################################
## 1. Packages
############################################################

required_pkgs <- c(
  "dplyr",
  "readr",
  "ggplot2",
  "ggrepel",
  "forcats",
  "tidyr"
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

standardize_results <- function(df, dataset_label, module_col, delta_col, p_col, fdr_col, logfc_col) {
  missing_cols <- setdiff(c(module_col, delta_col), colnames(df))
  if (length(missing_cols) > 0) {
    stop(dataset_label, " table is missing required columns: ", paste(missing_cols, collapse = ", "))
  }
  
  out <- df %>%
    rename(
      module = all_of(module_col),
      delta_percent = all_of(delta_col)
    )
  
  # Optional columns for downstream plots.
  if (p_col %in% colnames(out)) {
    out <- out %>% rename(P.Value = all_of(p_col))
  }
  if (fdr_col %in% colnames(out)) {
    out <- out %>% rename(adj.P.Val = all_of(fdr_col))
  }
  if (logfc_col %in% colnames(out)) {
    out <- out %>% rename(logFC = all_of(logfc_col))
  }
  
  out %>%
    mutate(
      module = as.character(module),
      dataset = dataset_label,
      sig = if ("adj.P.Val" %in% colnames(out)) adj.P.Val < 0.05 else FALSE,
      sig_label = if ("adj.P.Val" %in% colnames(out)) {
        case_when(
          adj.P.Val < 0.001 ~ "***",
          adj.P.Val < 0.01  ~ "**",
          adj.P.Val < 0.05  ~ "*",
          TRUE              ~ ""
        )
      } else {
        ""
      }
    )
}

save_both <- function(plot_obj, file_base, out_dir, width = 7, height = 5, dpi = 300) {
  ggsave(file.path(out_dir, paste0(file_base, ".png")), plot_obj, width = width, height = height, dpi = dpi)
  ggsave(file.path(out_dir, paste0(file_base, ".pdf")), plot_obj, width = width, height = height)
}

compute_correlations <- function(df, x_col = "delta_pm", y_col = "delta_pilot") {
  complete <- df %>% filter(!is.na(.data[[x_col]]), !is.na(.data[[y_col]]))
  if (nrow(complete) < 2) {
    return(list(
      n = nrow(complete),
      pearson = NA_real_,
      spearman = NA_real_
    ))
  }
  
  list(
    n = nrow(complete),
    pearson = cor(complete[[x_col]], complete[[y_col]], method = "pearson"),
    spearman = cor(complete[[x_col]], complete[[y_col]], method = "spearman")
  )
}

make_scatter_plot <- function(df,
                              title,
                              subtitle,
                              x_col = "delta_pm",
                              y_col = "delta_pilot",
                              label_col = "module",
                              color_col = "concordant",
                              module_colors = NULL,
                              label_all = FALSE,
                              top_n = 20,
                              highlight_genes = NULL) {
  df <- df %>%
    mutate(concordant = sign(.data[[x_col]]) == sign(.data[[y_col]]))
  
  plot_df <- df %>%
    filter(!is.na(.data[[x_col]]), !is.na(.data[[y_col]]))
  
  # Label strategy:
  #   - if label_all = TRUE, label everything (useful when the plot is small)
  #   - otherwise, label the top-N by absolute pilot effect and any highlighted genes.
  label_df <- if (label_all) {
    plot_df
  } else {
    plot_df %>% arrange(desc(abs(.data[[y_col]]))) %>% slice_head(n = top_n)
  }
  
  if (!is.null(highlight_genes)) {
    extra <- plot_df %>% filter(.data[[label_col]] %in% highlight_genes)
    label_df <- bind_rows(label_df, extra) %>% distinct(.data[[label_col]], .keep_all = TRUE)
  }
  
  corr <- compute_correlations(plot_df, x_col = x_col, y_col = y_col)
  cor_label <- paste0(
    "Pearson r = ", round(corr$pearson, 2),
    "\nSpearman ", "ρ = ", round(corr$spearman, 2),
    "\nn = ", corr$n
  )
  
  concordance_pct <- round(mean(plot_df$concordant) * 100, 1)
  
  p <- ggplot(plot_df, aes(x = .data[[x_col]], y = .data[[y_col]])) +
    geom_hline(yintercept = 0, linetype = "dashed", color = "grey70", linewidth = 0.4) +
    geom_vline(xintercept = 0, linetype = "dashed", color = "grey70", linewidth = 0.4) +
    geom_point(aes(color = concordant), alpha = 0.35, size = 0.8) +
    geom_smooth(method = "lm", se = TRUE, color = "black", linewidth = 0.6) +
    geom_text_repel(
      data = label_df,
      aes(label = .data[[label_col]]),
      size = 2.8,
      max.overlaps = if (label_all) Inf else 20,
      box.padding = 0.3,
      segment.color = "grey50",
      segment.size = 0.3
    ) +
    scale_color_manual(
      values = c("TRUE" = "#2166AC", "FALSE" = "#D6604D"),
      labels = c("TRUE" = "Concordant", "FALSE" = "Discordant"),
      name = NULL
    ) +
    annotate("text",
             x = -Inf, y = Inf,
             label = cor_label,
             hjust = -0.1, vjust = 1.2,
             size = 3.8
    ) +
    annotate("text",
             x = Inf, y = -Inf,
             label = paste0(concordance_pct, "% concordant direction"),
             hjust = 1.1, vjust = -0.5,
             size = 3.5, color = "grey40"
    ) +
    labs(x = paste0("Delta abundance % (", "PM", " − ", "Pilot", ")"),
         y = paste0("Delta abundance % (", "Pilot", " − ", "PM", ")"),
         title = title,
         subtitle = subtitle) +
    theme_classic(base_size = 12) +
    theme(
      plot.title = element_text(face = "bold"),
      plot.subtitle = element_text(size = 9, color = "grey40"),
      legend.position = "top"
    )
  
  if (!is.null(module_colors)) {
    # If module_colors is provided we can use them for text/point emphasis in later extensions.
    invisible(module_colors)
  }
  
  p
}

make_lollipop_plot <- function(df,
                               title,
                               subtitle,
                               module_levels,
                               module_colors = NULL,
                               pm_label = "Postmortem",
                               pilot_label = "miBrain Pilot") {
  plot_df <- df %>%
    mutate(
      module = factor(module, levels = module_levels),
      sig = ifelse(is.na(adj.P.Val), FALSE, adj.P.Val < 0.05),
      sig_label = ifelse(is.na(sig_label), "", sig_label)
    )
  
  long_df <- plot_df %>%
    select(module, delta_pm, delta_pilot, sig, sig_label) %>%
    pivot_longer(
      cols = c(delta_pm, delta_pilot),
      names_to = "dataset",
      values_to = "delta"
    ) %>%
    mutate(
      dataset = recode(dataset,
                       delta_pm = pm_label,
                       delta_pilot = pilot_label),
      direction_color = ifelse(delta > 0, "A/A > G/G", "G/G > A/A")
    )
  
  if (is.null(module_colors)) {
    module_colors <- setNames(rep("grey40", length(module_levels)), module_levels)
  }
  
  ggplot(long_df, aes(x = delta, y = module, color = direction_color, shape = dataset)) +
    geom_vline(xintercept = 0, linetype = "dashed", color = "grey60") +
    geom_line(aes(group = module), color = "grey70", linewidth = 0.4) +
    geom_point(size = 2.5, alpha = 0.9) +
    scale_color_manual(
      values = c("A/A > G/G" = "#2166AC", "G/G > A/G" = "#D6604D", "G/G > A/A" = "#D6604D"),
      name = "Direction"
    ) +
    scale_shape_manual(
      values = c(pm_label = 16, pilot_label = 17),
      name = "Dataset"
    ) +
    labs(
      title = title,
      subtitle = subtitle,
      x = "Delta abundance % (A/A − G/G)",
      y = NULL
    ) +
    theme_classic(base_size = 11) +
    theme(
      plot.title = element_text(face = "bold"),
      plot.subtitle = element_text(size = 8, color = "grey40"),
      axis.text.y = element_text(size = 8)
    )
}

############################################################
## 3. Load PM and pilot result tables
############################################################

pm_raw <- load_results_table(config$pm_results_path, label = "PM abundance results")
pilot_raw <- load_results_table(config$pilot_results_path, label = "pilot abundance results")

pm_results <- standardize_results(
  pm_raw,
  dataset_label = config$pm_label,
  module_col = config$module_col,
  delta_col = config$delta_col,
  p_col = config$p_col,
  fdr_col = config$fdr_col,
  logfc_col = config$logfc_col
)

pilot_results <- standardize_results(
  pilot_raw,
  dataset_label = config$pilot_label,
  module_col = config$module_col,
  delta_col = config$delta_col,
  p_col = config$p_col,
  fdr_col = config$fdr_col,
  logfc_col = config$logfc_col
)

# Keep only the columns we need for comparison.
pm_comp <- pm_results %>%
  select(module, delta_pm = delta_percent,
         PM_P = any_of("P.Value"),
         PM_FDR = any_of("adj.P.Val"),
         PM_logFC = any_of("logFC"))

pilot_comp <- pilot_results %>%
  select(module, delta_pilot = delta_percent,
         Pilot_P = any_of("P.Value"),
         Pilot_FDR = any_of("adj.P.Val"),
         Pilot_logFC = any_of("logFC"))

comparison <- full_join(pm_comp, pilot_comp, by = "module") %>%
  mutate(
    module = factor(module, levels = config$state_levels)
  ) %>%
  arrange(module)

# Simple concordance metrics on the delta percent effect sizes.
comparison <- comparison %>%
  mutate(concordant = sign(delta_pm) == sign(delta_pilot))

############################################################
## 4. Save comparison table
############################################################

comparison_table_path <- file.path(config$out_dir, config$comparison_table_name)
write_csv(comparison, comparison_table_path)
message("Saved comparison table: ", comparison_table_path)

############################################################
## 5. Concordance scatter plot
############################################################

scatter_title <- config$comparison_title
scatter_subtitle <- paste0(
  config$pm_label, " vs ", config$pilot_label,
  " | delta abundance comparison across scVI states"
)

# Main plot
p_scatter_all <- make_scatter_plot(
  df = comparison,
  title = scatter_title,
  subtitle = scatter_subtitle,
  x_col = "delta_pm",
  y_col = "delta_pilot",
  label_col = "module",
  color_col = "concordant",
  module_colors = config$state_colors,
  label_all = config$label_all_sig,
  top_n = config$top_n_labels
)

save_both(
  p_scatter_all,
  config$scatter_all_name,
  config$out_dir,
  width = 6.5,
  height = 6.0
)

# Filtered scatter using only modules with stronger pilot effect size
comparison_filt <- comparison %>%
  filter(abs(delta_pilot) >= config$beta_filter_abs | abs(delta_pm) >= config$beta_filter_abs)

p_scatter_filt <- make_scatter_plot(
  df = comparison_filt,
  title = paste0(config$comparison_title, " (filtered)"),
  subtitle = paste0(
    config$pm_label, " vs ", config$pilot_label,
    " | only modules with |delta| >= ", config$beta_filter_abs,
    " in either dataset"
  ),
  x_col = "delta_pm",
  y_col = "delta_pilot",
  label_col = "module",
  color_col = "concordant",
  module_colors = config$state_colors,
  label_all = TRUE,
  top_n = config$top_n_labels
)

save_both(
  p_scatter_filt,
  config$scatter_filt_name,
  config$out_dir,
  width = 6.5,
  height = 6.0
)

############################################################
## 6. Lollipop comparison plot
############################################################

lollipop_title <- "PM vs miBrain Pilot: differential abundance by state"
lollipop_subtitle <- paste0(
  "State-level delta abundance comparison | ",
  config$pm_label, " and ", config$pilot_label
)

p_lollipop <- make_lollipop_plot(
  df = comparison,
  title = lollipop_title,
  subtitle = lollipop_subtitle,
  module_levels = config$state_levels,
  module_colors = config$state_colors,
  pm_label = config$pm_label,
  pilot_label = config$pilot_label
)

save_both(
  p_lollipop,
  config$lollipop_name,
  config$out_dir,
  width = 9,
  height = 5.5
)

############################################################
## 7. Small summary tables for downstream use
############################################################

summary_table <- comparison %>%
  mutate(
    abs_delta_pm = abs(delta_pm),
    abs_delta_pilot = abs(delta_pilot)
  ) %>%
  select(module,
         delta_pm, delta_pilot,
         abs_delta_pm, abs_delta_pilot,
         concordant,
         PM_P, PM_FDR,
         Pilot_P, Pilot_FDR,
         PM_logFC, Pilot_logFC)

write_csv(summary_table, file.path(config$out_dir, "PM_vs_Pilot_abundance_summary_table.csv"))
message("Saved summary table.")

# Optional: ranked concordant and discordant subsets
ranked_concordant <- comparison %>%
  filter(concordant == TRUE) %>%
  arrange(desc(abs(delta_pilot)), desc(abs(delta_pm)))

ranked_discordant <- comparison %>%
  filter(concordant == FALSE) %>%
  arrange(desc(abs(delta_pilot)), desc(abs(delta_pm)))

write_csv(ranked_concordant, file.path(config$out_dir, "PM_vs_Pilot_abundance_concordant_ranked.csv"))
write_csv(ranked_discordant, file.path(config$out_dir, "PM_vs_Pilot_abundance_discordant_ranked.csv"))

############################################################
## 8. Final message
############################################################

message("\nScript 3 complete.")
message("Outputs written to: ", normalizePath(config$out_dir))
message("  - comparison table")
message("  - concordance scatter plot")
message("  - filtered concordance scatter plot")
message("  - lollipop comparison plot")
message("  - ranked concordant/discordant tables")
