

############################################################
## SCRIPT 2
## Postmortem abundance analysis of scVI states
##
## Purpose:
##   - load PM SummarizedExperiment / pseudobulk object
##   - compute PM state proportions
##   - run limma with brain_region covariate
##   - save PM abundance results
##   - make PM-only abundance plot
##   - optionally read pilot_results.csv and make PM vs pilot comparison plot
##
## Notes:
##   - This script does NOT recompute pilot results.
##   - This script does NOT run CAMERA, DEG correlation, or GSEA.
############################################################

############################################################
## 0. Configuration
############################################################

config <- list(
  out_dir = file.path(getwd(), "results", "script2_pm_abundance"),
  
  # PM input object can be:
  #   (1) already loaded in the environment under pm_object_name, or
  #   (2) read from pm_rds_path if provided.
  pm_object_name = "pb_MG",
  pm_rds_path    = NULL,
  
  # Optional saved pilot abundance results for comparison plots
  pilot_results_path = NULL,
  
  # Metadata columns in the PM object
  genotype_col    = "rs3732765",
  region_col      = "brain_region",
  
  # Optional metadata columns used in plotting / saving if present
  dataset_col     = "dataset",
  sample_col      = "sample",
  
  # Genotype values to compare
  allowed_genotypes = c("G/G", "A/A"),
  genotype_levels    = c("A/A", "G/G"),
  genotype_model_lvls = c("G/G", "A/A"),
  
  # State ordering must match the names used in the PM object assays
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
  
  # Colors for the abundance plots
  state_colors = c(
    "M1: P2RY12"   = "#F8766D",
    "M2: HSP90AA1" = "#C49A00",
    "M3: SPP1"     = "#53B400",
    "M4: TMEM163"  = "#00C1A2",
    "M5: PCDH9"    = "#00A6FF",
    "M6: CD163"    = "#A58AFF",
    "M7: FRMD4A"   = "#FB61D7",
    "M8: B2M"      = "#FF6F91"
  )
)

dir.create(config$out_dir, showWarnings = FALSE, recursive = TRUE)

############################################################
## 1. Packages
############################################################

required_pkgs <- c(
  "SummarizedExperiment",
  "dplyr",
  "tidyr",
  "ggplot2",
  "limma",
  "readr",
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
  library(dplyr)
  library(tidyr)
  library(ggplot2)
  library(limma)
  library(readr)
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

save_both <- function(plot_obj, file_base, out_dir, width = 7, height = 5, dpi = 300) {
  ggsave(file.path(out_dir, paste0(file_base, ".png")), plot_obj, width = width, height = height, dpi = dpi)
  ggsave(file.path(out_dir, paste0(file_base, ".pdf")), plot_obj, width = width, height = height)
}

extract_pm_state_signals <- function(pb, state_levels) {
  if (!inherits(pb, "SummarizedExperiment")) {
    stop("PM input must be a SummarizedExperiment-like object.")
  }
  
  assays_here <- assayNames(pb)
  if (length(assays_here) == 0) {
    stop("No assays found in the PM object.")
  }
  
  # Keep only assays that correspond to desired states, when possible.
  state_assays <- intersect(state_levels, assays_here)
  if (length(state_assays) == 0) {
    # If assay names do not match state labels, fall back to all assays.
    state_assays <- assays_here
    message("Assay names do not match configured state labels; using all assays found in the object.")
  }
  
  pm_mod <- lapply(state_assays, function(m) {
    mat <- assay(pb, m)
    if (nrow(mat) == 0 || ncol(mat) == 0) {
      stop("Assay '", m, "' is empty.")
    }
    
    sample_signal <- colMeans(mat)
    data.frame(
      sample = names(sample_signal),
      module = m,
      module_signal = as.numeric(sample_signal),
      stringsAsFactors = FALSE
    )
  }) |>
    bind_rows()
  
  pm_mod
}

make_pm_metadata <- function(pb, genotype_col, region_col) {
  meta <- as.data.frame(colData(pb))
  meta$sample <- rownames(meta)
  
  needed <- c(genotype_col, region_col)
  missing <- needed[!needed %in% colnames(meta)]
  if (length(missing) > 0) {
    stop("Missing PM metadata columns: ", paste(missing, collapse = ", "))
  }
  
  meta
}

build_pm_state_tables <- function(pm_mod, pm_meta, state_levels, allowed_genotypes, genotype_levels, genotype_model_lvls,
                                  genotype_col = "rs3732765", region_col = "brain_region") {
  pm_mod_all <- pm_mod %>%
    left_join(pm_meta, by = "sample") %>%
    filter(.data[[genotype_col]] %in% allowed_genotypes) %>%
    mutate(
      genotype_plot  = factor(.data[[genotype_col]], levels = genotype_levels),
      genotype_model = factor(.data[[genotype_col]], levels = genotype_model_lvls),
      brain_region    = factor(.data[[region_col]]),
      module          = factor(module, levels = state_levels)
    ) %>%
    group_by(sample) %>%
    mutate(prop = module_signal / sum(module_signal)) %>%
    ungroup()
  
  if (nrow(pm_mod_all) == 0) {
    stop("No PM samples remained after genotype filtering.")
  }
  
  message("=== PM sample counts by genotype and brain region ===")
  print(
    pm_mod_all %>%
      distinct(sample, genotype_model, brain_region) %>%
      count(genotype_model, brain_region)
  )
  
  pm_state_wide <- pm_mod_all %>%
    select(sample, genotype_plot, genotype_model, brain_region, module, prop) %>%
    pivot_wider(names_from = module, values_from = prop, values_fill = 0)
  
  for (st in state_levels) {
    if (!st %in% colnames(pm_state_wide)) {
      pm_state_wide[[st]] <- 0
    }
  }
  
  pm_state_wide <- pm_state_wide %>%
    select(sample, genotype_plot, genotype_model, brain_region, all_of(state_levels))
  
  pm_bar_table <- pm_mod_all %>%
    group_by(genotype_plot, module) %>%
    summarise(mean_prop = mean(prop), .groups = "drop") %>%
    group_by(genotype_plot) %>%
    mutate(
      frac_bar = mean_prop / sum(mean_prop),
      percent_bar = frac_bar * 100
    ) %>%
    ungroup() %>%
    arrange(genotype_plot, module)
  
  pm_prop_mat <- t(as.matrix(pm_state_wide[, state_levels]))
  colnames(pm_prop_mat) <- pm_state_wide$sample
  rownames(pm_prop_mat) <- state_levels
  
  list(
    pm_mod_all = pm_mod_all,
    pm_state_wide = pm_state_wide,
    pm_bar_table = pm_bar_table,
    pm_prop_mat = pm_prop_mat
  )
}

fit_pm_abundance <- function(pm_state_wide, pm_prop_mat, state_levels) {
  pm_prop_mat_tr <- asin(sqrt(pm_prop_mat))
  
  design_pm <- model.matrix(~ genotype_model + brain_region, data = pm_state_wide)
  colnames(design_pm) <- make.names(colnames(design_pm))
  
  fit_pm <- lmFit(pm_prop_mat_tr, design_pm)
  fit_pm <- eBayes(fit_pm)
  
  coef_pm_name <- grep("genotype_model", colnames(design_pm), value = TRUE)
  coef_pm_name <- setdiff(coef_pm_name, "(Intercept)")
  if (length(coef_pm_name) != 1) {
    stop("Could not identify genotype coefficient. Found: ", paste(coef_pm_name, collapse = ", "))
  }
  
  pm_limma <- topTable(fit_pm, coef = coef_pm_name, number = Inf, sort.by = "none")
  pm_limma$module <- rownames(pm_limma)
  
  pm_means <- pm_state_wide %>%
    select(sample, genotype_model, all_of(state_levels)) %>%
    pivot_longer(cols = all_of(state_levels), names_to = "module", values_to = "prop") %>%
    group_by(module, genotype_model) %>%
    summarise(mean_prop = mean(prop), .groups = "drop") %>%
    pivot_wider(names_from = genotype_model, values_from = mean_prop)
  
  # robustly rename genotype columns regardless of syntactic conversion
  colnames(pm_means) <- gsub("/", ".", colnames(pm_means))
  if ("G.G" %in% colnames(pm_means)) colnames(pm_means)[colnames(pm_means) == "G.G"] <- "mean_prop_GG"
  if ("A.A" %in% colnames(pm_means)) colnames(pm_means)[colnames(pm_means) == "A.A"] <- "mean_prop_AA"
  
  if (!all(c("mean_prop_GG", "mean_prop_AA") %in% colnames(pm_means))) {
    stop("Could not find genotype columns in pm_means. Found: ", paste(colnames(pm_means), collapse = ", "))
  }
  
  pm_means <- pm_means %>%
    mutate(
      delta_prop = mean_prop_AA - mean_prop_GG,
      delta_percent = delta_prop * 100
    )
  
  pm_results <- pm_limma %>%
    left_join(pm_means, by = "module") %>%
    select(module, logFC, AveExpr, t, P.Value, adj.P.Val, B,
           mean_prop_GG, mean_prop_AA, delta_prop, delta_percent) %>%
    arrange(match(module, state_levels))
  
  list(
    pm_limma = pm_limma,
    pm_results = pm_results,
    design_pm = design_pm
  )
}

make_abundance_lollipop <- function(plot_df, title, subtitle, col_AA = "#D55E00", col_GG = "#0072B2") {
  ggplot(plot_df, aes(x = delta_percent, y = module)) +
    geom_vline(xintercept = 0, linetype = "dashed", color = "grey50", linewidth = 0.5) +
    geom_segment(
      aes(x = 0, xend = delta_percent, y = module, yend = module,
          color = ifelse(delta_percent > 0, "A/A enriched", "G/G enriched")),
      linewidth = 0.8
    ) +
    geom_point(
      aes(color = ifelse(delta_percent > 0, "A/A enriched", "G/G enriched"),
          size = sig,
          shape = sig)
    ) +
    geom_text(
      aes(x = delta_percent + ifelse(delta_percent >= 0, 0.3, -0.3), label = sig_label),
      hjust = ifelse(plot_df$delta_percent >= 0, 0, 1),
      size = 4.5,
      color = "grey20"
    ) +
    scale_color_manual(
      values = c("A/A enriched" = col_AA, "G/G enriched" = col_GG),
      name = "Higher in"
    ) +
    scale_size_manual(values = c(`TRUE` = 4, `FALSE` = 2.5), guide = "none") +
    scale_shape_manual(values = c(`TRUE` = 16, `FALSE` = 21), guide = "none") +
    scale_x_continuous(expand = expansion(mult = c(0.15, 0.15))) +
    labs(x = "Delta abundance % (A/A − G/G)", y = NULL, title = title, subtitle = subtitle) +
    theme_classic(base_size = 12) +
    theme(
      plot.title = element_text(face = "bold", size = 13),
      plot.subtitle = element_text(size = 8, color = "grey40"),
      legend.position = "top",
      axis.text.y = element_text(size = 10)
    )
}

make_combined_abundance_plot <- function(pm_results, pilot_results, state_levels) {
  pm_plot_df <- pm_results %>%
    mutate(
      sig = adj.P.Val < 0.05,
      sig_label = case_when(
        adj.P.Val < 0.001 ~ "***",
        adj.P.Val < 0.01 ~ "**",
        adj.P.Val < 0.05 ~ "*",
        TRUE ~ ""
      ),
      dataset = "Postmortem",
      module = factor(module, levels = rev(state_levels))
    )
  
  pilot_plot_df <- pilot_results %>%
    mutate(
      sig = adj.P.Val < 0.05,
      sig_label = case_when(
        adj.P.Val < 0.001 ~ "***",
        adj.P.Val < 0.01 ~ "**",
        adj.P.Val < 0.05 ~ "*",
        TRUE ~ ""
      ),
      dataset = "miBrain Pilot",
      module = factor(module, levels = rev(state_levels))
    )
  
  combined <- bind_rows(pm_plot_df, pilot_plot_df) %>%
    mutate(
      dataset = factor(dataset, levels = c("Postmortem", "miBrain Pilot")),
      y_num = as.numeric(module),
      y_offset = ifelse(dataset == "Postmortem", 0.18, -0.18),
      y_dodge = y_num + y_offset,
      ds_color = ifelse(dataset == "Postmortem", "#00836B", "#AA4499")
    )
  
  ggplot(combined) +
    geom_vline(xintercept = 0, linetype = "dashed", color = "grey50", linewidth = 0.5) +
    geom_segment(
      aes(x = 0, xend = delta_percent, y = y_dodge, yend = y_dodge,
          color = dataset, linetype = dataset),
      linewidth = 1
    ) +
    geom_point(
      aes(x = delta_percent, y = y_dodge,
          color = dataset, shape = dataset, size = sig,
          fill = ifelse(sig, ds_color, NA))
    ) +
    geom_text(
      aes(x = delta_percent + ifelse(delta_percent >= 0, 0.8, -0.8),
          y = y_dodge, label = sig_label, color = dataset),
      hjust = ifelse(combined$delta_percent >= 0, 0, 1),
      size = 4.5,
      show.legend = FALSE
    ) +
    scale_y_continuous(breaks = seq_along(rev(state_levels)), labels = rev(state_levels), expand = expansion(add = 0.6)) +
    scale_color_manual(values = c("Postmortem" = "#00836B", "miBrain Pilot" = "#AA4499"), name = NULL) +
    scale_fill_manual(values = c("Postmortem" = "#00836B", "miBrain Pilot" = "#AA4499"), guide = "none") +
    scale_linetype_manual(values = c("Postmortem" = "solid", "miBrain Pilot" = "dashed"), name = NULL) +
    scale_shape_manual(values = c("Postmortem" = 21, "miBrain Pilot" = 24), name = NULL) +
    scale_size_manual(values = c(`TRUE` = 4, `FALSE` = 2.5), guide = "none") +
    scale_x_continuous(expand = expansion(mult = c(0.15, 0.18))) +
    labs(
      title = "Differential microglial state abundance: A/A vs G/G",
      subtitle = paste0(
        "Teal (circles, solid) = Postmortem | Pink (triangles, dashed) = miBrain Pilot\n",
        "Filled = FDR<0.05 | Open = not significant | * FDR<0.05  ** FDR<0.01  *** FDR<0.001"
      ),
      x = "Delta abundance % (A/A − G/G)",
      y = NULL
    ) +
    theme_classic(base_size = 12) +
    theme(
      plot.title = element_text(face = "bold", size = 13),
      plot.subtitle = element_text(size = 8, color = "grey40"),
      legend.position = "top",
      legend.text = element_text(size = 11),
      legend.key.width = unit(1.5, "cm"),
      axis.text.y = element_text(size = 10)
    )
}

############################################################
## 3. Load PM object
############################################################

pm <- load_object(config$pm_object_name, config$pm_rds_path, label = "PM object")

if (!inherits(pm, "SummarizedExperiment")) {
  stop("PM object must be a SummarizedExperiment-like object.")
}

message("Loaded PM object: ", paste(class(pm), collapse = ", "))
message("Assays found: ", paste(assayNames(pm), collapse = ", "))

############################################################
## 4. Build PM abundance tables
############################################################

pm_meta <- make_pm_metadata(pm, genotype_col = config$genotype_col, region_col = config$region_col)
pm_mod  <- extract_pm_state_signals(pm, config$state_levels)

pm_tables <- build_pm_state_tables(
  pm_mod = pm_mod,
  pm_meta = pm_meta,
  state_levels = config$state_levels,
  allowed_genotypes = config$allowed_genotypes,
  genotype_levels = config$genotype_levels,
  genotype_model_lvls = config$genotype_model_lvls,
  genotype_col = config$genotype_col,
  region_col = config$region_col
)

pm_mod_all    <- pm_tables$pm_mod_all
pm_state_wide <- pm_tables$pm_state_wide
pm_bar_table  <- pm_tables$pm_bar_table
pm_prop_mat   <- pm_tables$pm_prop_mat

write_csv(pm_bar_table, file.path(config$out_dir, "PM_AA_vs_GG_stackedbar_percentages.csv"))
message("Saved PM stacked-bar summary table.")

############################################################
## 5. Fit PM abundance model
############################################################

pm_fit <- fit_pm_abundance(pm_state_wide, pm_prop_mat, config$state_levels)
pm_limma   <- pm_fit$pm_limma
pm_results <- pm_fit$pm_results

glimpse(pm_results)
write_csv(pm_results, file.path(config$out_dir, "PM_abundance_limma_brainRegionAdj_GGvsAA.csv"))
message("Saved PM abundance results.")

############################################################
## 6. PM-only plot
############################################################

pm_plot_df <- pm_results %>%
  mutate(
    sig = adj.P.Val < 0.05,
    sig_label = case_when(
      adj.P.Val < 0.001 ~ "***",
      adj.P.Val < 0.01 ~ "**",
      adj.P.Val < 0.05 ~ "*",
      TRUE ~ ""
    ),
    module = factor(module, levels = rev(config$state_levels))
  )

p_pm_only <- make_abundance_lollipop(
  plot_df = pm_plot_df,
  title = "Postmortem microglial state abundance: A/A vs G/G",
  subtitle = paste0(
    "n = ", nrow(distinct(pm_mod_all, sample)), " samples | ",
    "limma on arcsin-sqrt transformed proportions | adjusted for brain region\n",
    "* FDR<0.05  ** FDR<0.01  *** FDR<0.001 | filled = FDR<0.05"
  )
)

save_both(p_pm_only, "PM_abundance_lollipop_brainRegionAdj", config$out_dir, width = 7, height = 5)

############################################################
## 7. Optional comparison plot using saved pilot results
############################################################

pilot_results <- NULL
if (!is.null(config$pilot_results_path) && file.exists(config$pilot_results_path)) {
  pilot_results <- read_csv(config$pilot_results_path, show_col_types = FALSE)
  required_cols <- c("module", "delta_percent", "P.Value", "adj.P.Val")
  missing_cols <- required_cols[!required_cols %in% colnames(pilot_results)]
  if (length(missing_cols) > 0) {
    stop("pilot_results file is missing required columns: ", paste(missing_cols, collapse = ", "))
  }
  
  p_combined <- make_combined_abundance_plot(pm_results, pilot_results, config$state_levels)
  save_both(p_combined, "Combined_PM_Pilot_abundance_lollipop_FIXED", config$out_dir, width = 9, height = 5.5)
  
  comparison_table <- pm_results %>%
    select(module,
           PM_delta_percent = delta_percent,
           PM_P = P.Value,
           PM_FDR = adj.P.Val) %>%
    left_join(
      pilot_results %>%
        select(module,
               Pilot_delta_percent = delta_percent,
               Pilot_P = P.Value,
               Pilot_FDR = adj.P.Val),
      by = "module"
    ) %>%
    arrange(match(module, config$state_levels))
  
  write_csv(comparison_table, file.path(config$out_dir, "PM_vs_Pilot_abundance_comparison_table.csv"))
  message("Saved PM vs pilot comparison table and plot.")
} else {
  message("pilot_results_path not provided or file not found; skipping PM vs pilot comparison plot.")
}

############################################################
## 8. Final outputs
############################################################

message("\nScript 2 complete.")
message("Outputs written to: ", normalizePath(config$out_dir))
message("  - PM abundance summary table")
message("  - PM limma results")
message("  - PM-only abundance plot")
message("  - optional PM vs pilot comparison plot if pilot_results_path was supplied")
