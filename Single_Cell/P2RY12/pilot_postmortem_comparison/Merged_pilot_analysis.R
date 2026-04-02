############################################################
## SCRIPT 1
## Pilot miBrain pipeline: merge, export, pseudobulk, DE,
## optional scVI-label summaries
##
## Design goal:
##   - keep this script pilot-only
##   - no PM comparisons here
##   - no hardcoded project-specific object names unless needed
##   - save standardized outputs for downstream scripts
############################################################

############################################################
## 0. Configuration
############################################################

config <- list(
  out_dir = file.path(getwd(), "results", "script1_pilot_pipeline"),
  out_prefix = "merged_microglia",
  
  # Input Seurat objects can be either:
  #   (1) already loaded in the environment under these names, or
  #   (2) read from RDS paths if the paths are provided.
  sc_object_name  = "microglia_subset_pilot",
  nuc_object_name = "microglia_subset_nuclei",
  sc_rds_path     = NULL,
  nuc_rds_path    = NULL,
  
  # Optional scVI-labeled objects for stacked-bar plots
  scvi_sc_object_name  = "scvi_cells",
  scvi_nuc_object_name = "scvi_nuclei",
  scvi_sc_rds_path     = NULL,
  scvi_nuc_rds_path    = NULL,
  
  # Metadata columns
  genotype_col = "genotype",
  origin_col    = "origin",
  state_col     = "C_scANVI",
  
  # Values to keep for pilot DE
  allowed_genotypes = c("G/G", "A/A"),
  genotype_map = c("G/G" = "GG", "A/A" = "AA"),
  
  # Labels used in merged object / pseudobulk
  origin_map = c("single_cell" = "SC", "single_nuclei" = "NUC"),
  
  # scVI state ordering (edit once here, not in every plot)
  state_levels = c(
    "M1: P2RY12",
    "M2: HSP90AA1",
    "M3: SPP1",
    "M4: TMEM163",
    "M5: PCDH9",
    "M6: CD163",
    "M7: FRMD4A",
    "M8: B2M"
  )
)

dir.create(config$out_dir, showWarnings = FALSE, recursive = TRUE)

############################################################
## 1. Packages
############################################################

required_pkgs <- c(
  "Seurat",
  "Matrix",
  "data.table",
  "edgeR",
  "dplyr",
  "ggplot2",
  "tidyr",
  "limma",
  "variancePartition",
  "readr"
)

missing_pkgs <- required_pkgs[!vapply(required_pkgs, requireNamespace, logical(1), quietly = TRUE)]
if (length(missing_pkgs) > 0) {
  stop(
    "Missing packages: ", paste(missing_pkgs, collapse = ", "),
    "\nInstall them first (ideally with renv) before running this script."
  )
}

suppressPackageStartupMessages({
  library(Seurat)
  library(Matrix)
  library(data.table)
  library(edgeR)
  library(dplyr)
  library(ggplot2)
  library(tidyr)
  library(limma)
  library(variancePartition)
  library(readr)
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

pick_assay <- function(seu) {
  assays <- Assays(seu)
  if ("RNA" %in% assays) {
    return("RNA")
  }
  assays[[1]]
}

get_counts <- function(seu, assay = NULL) {
  assay <- assay %||% pick_assay(seu)
  counts <- tryCatch(
    GetAssayData(seu, assay = assay, layer = "counts"),
    error = function(e) GetAssayData(seu, assay = assay, slot = "counts")
  )
  if (!inherits(counts, "dgCMatrix")) {
    counts <- as(counts, "dgCMatrix")
  }
  counts
}

make_origin_prefix <- function(cell_names, origin_map = c("single_cell" = "SC", "single_nuclei" = "NUC")) {
  out <- rep(NA_character_, length(cell_names))
  out[startsWith(cell_names, paste0(names(origin_map)[1], "_"))] <- origin_map[[1]]
  if (length(origin_map) > 1) {
    for (nm in names(origin_map)[-1]) {
      out[startsWith(cell_names, paste0(nm, "_"))] <- origin_map[[nm]]
    }
  }
  out
}

merge_seurat_objects <- function(sc, nuc, out_dir, out_prefix, origin_map) {
  assay_sc  <- pick_assay(sc)
  assay_nuc <- pick_assay(nuc)
  
  genes_sc  <- rownames(get_counts(sc, assay_sc))
  genes_nuc <- rownames(get_counts(nuc, assay_nuc))
  
  common_genes <- intersect(genes_sc, genes_nuc)
  if (length(common_genes) < 1) {
    stop("No overlapping genes between SC and NUC. Check gene IDs (symbols vs Ensembl).")
  }
  
  message("Features in SC:  ", length(genes_sc))
  message("Features in NUC: ", length(genes_nuc))
  message("Using ", length(common_genes), " common genes for merging.")
  
  sc  <- subset(sc, features = common_genes)
  nuc <- subset(nuc, features = common_genes)
  
  merged <- merge(
    sc,
    y = nuc,
    add.cell.ids = c(names(origin_map)[1], names(origin_map)[2]),
    project = out_prefix
  )
  
  # Add a simple origin label inferred from the prefix that merge() adds.
  merged$origin <- NA_character_
  merged$origin[startsWith(colnames(merged), paste0(names(origin_map)[1], "_"))] <- names(origin_map)[1]
  merged$origin[startsWith(colnames(merged), paste0(names(origin_map)[2], "_"))] <- names(origin_map)[2]
  merged$origin_prefix <- unname(origin_map[merged$origin])
  
  rds_out <- file.path(out_dir, paste0(out_prefix, ".rds"))
  saveRDS(merged, rds_out)
  message("Saved merged Seurat object: ", rds_out)
  
  counts <- get_counts(merged)
  mtx_path      <- file.path(out_dir, paste0(out_prefix, "_matrix.mtx"))
  features_path <- file.path(out_dir, paste0(out_prefix, "_features.tsv"))
  barcodes_path  <- file.path(out_dir, paste0(out_prefix, "_barcodes.tsv"))
  meta_path      <- file.path(out_dir, paste0(out_prefix, "_metadata.tsv"))
  
  Matrix::writeMM(counts, file = mtx_path)
  
  features_df <- data.frame(
    feature_id   = rownames(counts),
    feature_name = rownames(counts),
    stringsAsFactors = FALSE
  )
  fwrite(features_df, file = features_path, sep = "\t", col.names = FALSE)
  fwrite(data.frame(barcode = colnames(counts)), file = barcodes_path, sep = "\t", col.names = FALSE)
  
  meta_df <- cbind(barcode = rownames(merged@meta.data), merged@meta.data)
  fwrite(meta_df, file = meta_path, sep = "\t", na = "NA")
  
  message("Exported MTX/TSV bundle to ", normalizePath(out_dir))
  message("  features: ", nrow(features_df),
          " | barcodes: ", ncol(counts),
          " | metadata rows: ", nrow(meta_df))
  
  list(
    merged = merged,
    counts = counts,
    common_genes = common_genes,
    files = list(
      rds = rds_out,
      mtx = mtx_path,
      features = features_path,
      barcodes = barcodes_path,
      metadata = meta_path
    )
  )
}

build_pseudobulk_matrix <- function(seu, genotype_col, allowed_genotypes, genotype_map) {
  meta <- as.data.frame(seu@meta.data)
  if (!genotype_col %in% colnames(meta)) {
    stop("Missing genotype column: ", genotype_col)
  }
  if (!"origin_prefix" %in% colnames(meta)) {
    if (!"origin" %in% colnames(meta)) {
      stop("Missing origin metadata; cannot build pseudobulk.")
    }
    meta$origin_prefix <- NA_character_
    meta$origin_prefix[meta$origin == names(genotype_map)[1]] <- unname(c("SC", "NUC")[1])
  }
  
  meta$gen_std <- NA_character_
  meta$gen_std[meta[[genotype_col]] == names(genotype_map)[1]] <- genotype_map[[1]]
  meta$gen_std[meta[[genotype_col]] == names(genotype_map)[2]] <- genotype_map[[2]]
  
  keep_cells <- rownames(meta)[!is.na(meta$gen_std) & meta[[genotype_col]] %in% allowed_genotypes]
  if (length(keep_cells) == 0) {
    stop("No cells found for the requested genotypes: ", paste(allowed_genotypes, collapse = ", "))
  }
  
  seu2 <- subset(seu, cells = keep_cells)
  meta2 <- as.data.frame(seu2@meta.data)
  
  if (!"origin_prefix" %in% colnames(meta2)) {
    meta2$origin_prefix <- NA_character_
    meta2$origin_prefix[startsWith(rownames(meta2), "SC_")]  <- "SC"
    meta2$origin_prefix[startsWith(rownames(meta2), "NUC_")] <- "NUC"
  }
  
  meta2$gen_std <- NA_character_
  meta2$gen_std[meta2[[genotype_col]] == names(genotype_map)[1]] <- genotype_map[[1]]
  meta2$gen_std[meta2[[genotype_col]] == names(genotype_map)[2]] <- genotype_map[[2]]
  
  message("Origin x genotype table:")
  print(table(meta2$origin_prefix, meta2$gen_std, useNA = "ifany"))
  
  counts <- get_counts(seu2)
  pb_samples <- c("SC_GG", "SC_AA", "NUC_GG", "NUC_AA")
  
  pb_mat <- sapply(pb_samples, function(s) {
    parts <- strsplit(s, "_", fixed = TRUE)[[1]]
    opfx <- parts[1]
    gen  <- parts[2]
    
    cells <- rownames(meta2)[meta2$origin_prefix == opfx & meta2$gen_std == gen]
    if (length(cells) == 0) {
      stop("No cells found for pseudobulk sample ", s,
           ". Check the origin/genotype table above.")
    }
    Matrix::rowSums(counts[, cells, drop = FALSE])
  })
  
  rownames(pb_mat) <- rownames(counts)
  colnames(pb_mat) <- pb_samples
  
  message("Pseudobulk matrix: ", nrow(pb_mat), " genes x ", ncol(pb_mat), " samples")
  pb_df <- cbind(gene = rownames(pb_mat), as.data.frame(pb_mat))
  list(pb_mat = pb_mat, pb_df = pb_df, meta = meta2)
}

run_edger_de <- function(pb_mat, out_dir, out_file = "DE_edgeR_AA_vs_GG.tsv") {
  group <- factor(c("GG", "AA", "GG", "AA"), levels = c("GG", "AA"))
  names(group) <- colnames(pb_mat)
  
  dge <- DGEList(counts = pb_mat, group = group)
  keep <- rowSums(cpm(dge) > 1) >= 2
  message("edgeR: keeping ", sum(keep), " genes of ", nrow(dge$counts), " after CPM filtering.")
  
  dge <- dge[keep, , keep.lib.sizes = FALSE]
  dge <- calcNormFactors(dge)
  
  design <- model.matrix(~ group)
  colnames(design) <- gsub("^group", "", colnames(design))
  
  dge <- estimateDisp(dge, design)
  fit <- glmQLFit(dge, design)
  
  coef_name <- colnames(design)[grepl("AA", colnames(design))]
  if (length(coef_name) != 1) {
    stop("Could not uniquely identify AA coefficient. Design columns: ",
         paste(colnames(design), collapse = ", "))
  }
  coef_idx <- which(colnames(design) == coef_name)
  qlf <- glmQLFTest(fit, coef = coef_idx)
  res <- as.data.frame(topTags(qlf, n = nrow(dge$counts), sort.by = "PValue"))
  
  fwrite(res, file = file.path(out_dir, out_file), sep = "\t")
  saveRDS(
    list(dge = dge, design = design, fit = fit, qlf = qlf, res = res, pb_mat = pb_mat),
    file = file.path(out_dir, "edgeR_objects_AA_vs_GG.rds")
  )
  
  message("Saved edgeR results and RDS objects.")
  list(res = res, dge = dge, design = design, fit = fit, qlf = qlf)
}

run_dream_de <- function(pb_mat, out_dir, out_file = "DE_DREAM_AA_vs_GG.tsv") {
  sample_info <- data.frame(
    sample = colnames(pb_mat),
    genotype = factor(c("GG", "AA", "GG", "AA"), levels = c("GG", "AA")),
    modality = factor(c("SC", "SC", "NUC", "NUC"), levels = c("SC", "NUC"))
  )
  rownames(sample_info) <- sample_info$sample
  
  dge <- DGEList(counts = pb_mat)
  keep <- rowSums(cpm(dge) > 1) >= 2
  message("DREAM: keeping ", sum(keep), " genes of ", nrow(dge$counts), " after CPM filtering.")
  
  dge <- dge[keep, , keep.lib.sizes = FALSE]
  dge <- calcNormFactors(dge)
  
  form <- ~ genotype + modality
  v <- voomWithDreamWeights(dge, form, sample_info)
  fit <- dream(v, form, sample_info)
  fit <- eBayes(fit)
  
  coef_name <- "genotypeAA"
  if (!coef_name %in% colnames(coef(fit))) {
    stop("Could not find coefficient '", coef_name, "'. Available: ",
         paste(colnames(coef(fit)), collapse = ", "))
  }
  
  res <- topTable(fit, coef = coef_name, number = Inf, sort.by = "P")
  fwrite(res, file = file.path(out_dir, out_file), sep = "\t")
  
  res_with_symbol <- res
  res_with_symbol$symbol <- rownames(res_with_symbol)
  res_with_symbol <- res_with_symbol[, c("symbol", setdiff(colnames(res_with_symbol), "symbol"))]
  fwrite(
    res_with_symbol,
    file = file.path(out_dir, sub("\\.tsv$", "_withSymbol.tsv", out_file)),
    sep = "\t"
  )
  
  message("Saved DREAM results.")
  list(res = res, res_with_symbol = res_with_symbol, fit = fit, sample_info = sample_info)
}

attach_scvi_metadata <- function(seu, scvi_obj, cols) {
  if (is.null(scvi_obj)) {
    stop("scVI object is NULL.")
  }
  meta <- as.data.frame(scvi_obj@meta.data)
  needed <- cols[!cols %in% colnames(meta)]
  if (length(needed) > 0) {
    stop("Missing scVI metadata columns: ", paste(needed, collapse = ", "))
  }
  if (is.null(rownames(meta))) {
    stop("scVI metadata must have cell barcodes as rownames.")
  }
  
  # Align to Seurat cell order.
  meta <- meta[colnames(seu), cols, drop = FALSE]
  if (anyNA(rownames(meta))) {
    stop("Some Seurat cells were not found in the scVI metadata.")
  }
  AddMetaData(seu, meta)
}

make_state_bar_table <- function(meta, group_col, state_col, group_levels, state_levels) {
  df <- as.data.frame(meta)
  if (!group_col %in% colnames(df)) stop("Missing group column: ", group_col)
  if (!state_col %in% colnames(df)) stop("Missing state column: ", state_col)
  
  df <- df %>%
    filter(!is.na(.data[[state_col]]), .data[[group_col]] %in% group_levels) %>%
    mutate(
      !!group_col := factor(.data[[group_col]], levels = group_levels),
      !!state_col := factor(.data[[state_col]], levels = state_levels)
    ) %>%
    count(.data[[group_col]], .data[[state_col]], name = "n") %>%
    group_by(.data[[group_col]]) %>%
    mutate(
      prop = n / sum(n),
      percent = round(prop * 100, 1)
    ) %>%
    ungroup()
  
  df
}

make_state_stack_plot <- function(bar_df, group_col, state_col, colors, title, xlab, ylab) {
  ggplot(bar_df, aes(x = .data[[group_col]], y = prop, fill = .data[[state_col]])) +
    geom_col() +
    scale_fill_manual(values = colors) +
    labs(x = xlab, y = ylab, fill = "scVI state", title = title) +
    theme_classic(base_size = 12) +
    theme(
      panel.grid = element_blank(),
      axis.text.x = element_text(size = 12),
      legend.position = "right"
    )
}

save_plot_both <- function(plot_obj, file_base, out_dir, width = 5, height = 4, dpi = 300) {
  ggsave(file.path(out_dir, paste0(file_base, ".png")), plot_obj, width = width, height = height, dpi = dpi)
  ggsave(file.path(out_dir, paste0(file_base, ".pdf")), plot_obj, width = width, height = height)
}

############################################################
## 3. Load input objects
############################################################

sc  <- load_object(config$sc_object_name,  config$sc_rds_path,  label = "pilot sc object")
nuc <- load_object(config$nuc_object_name, config$nuc_rds_path, label = "pilot nuclei object")

############################################################
## 4. Merge Seurat objects and export MTX + TSV bundle
############################################################

merged_bundle <- merge_seurat_objects(
  sc = sc,
  nuc = nuc,
  out_dir = config$out_dir,
  out_prefix = config$out_prefix,
  origin_map = config$origin_map
)

merged <- merged_bundle$merged
counts <- merged_bundle$counts

############################################################
## 5. Pseudobulk (SC/NUC x genotype)
############################################################

pb_out <- build_pseudobulk_matrix(
  seu = merged,
  genotype_col = config$genotype_col,
  allowed_genotypes = config$allowed_genotypes,
  genotype_map = config$genotype_map
)

pb_mat <- pb_out$pb_mat
pb_df  <- pb_out$pb_df

fwrite(pb_df, file = file.path(config$out_dir, "pseudobulk_counts_matrix.tsv"), sep = "\t")
message("Saved pseudobulk counts matrix.")

############################################################
## 6. edgeR DE: A/A vs G/G
############################################################

edgeR_out <- run_edger_de(pb_mat, config$out_dir)
message("Top edgeR genes:")
print(head(edgeR_out$res, 10))

############################################################
## 7. DREAM DE: A/A vs G/G with modality covariate
############################################################

dream_out <- run_dream_de(pb_mat, config$out_dir)
message("Top DREAM genes:")
print(head(dream_out$res, 10))

############################################################
## 8. Optional: attach scVI labels to pilot Seurat objects
##    and save labeled RDS files
############################################################

if (exists(config$scvi_sc_object_name, inherits = TRUE) && exists(config$scvi_nuc_object_name, inherits = TRUE)) {
  scvi_sc  <- get(config$scvi_sc_object_name,  inherits = TRUE)
  scvi_nuc <- get(config$scvi_nuc_object_name, inherits = TRUE)
  
  scvi_cols <- c("C_scANVI", "prediction_confidence")
  
  microglia_subset_pilot_labeled <- attach_scvi_metadata(
    seu = sc,
    scvi_obj = scvi_sc,
    cols = scvi_cols
  )
  
  microglia_subset_nuclei_labeled <- attach_scvi_metadata(
    seu = nuc,
    scvi_obj = scvi_nuc,
    cols = scvi_cols
  )
  
  saveRDS(
    microglia_subset_pilot_labeled,
    file.path(config$out_dir, "microglia_subset_pilot_scvi_labeled.rds")
  )
  saveRDS(
    microglia_subset_nuclei_labeled,
    file.path(config$out_dir, "microglia_subset_nuclei_scvi_labeled.rds")
  )
  message("Saved labeled Seurat objects.")
  
  ##########################################################
  ## 9. Pilot-only scVI state composition plots
  ##########################################################
  
  # Single-cell
  df_sc <- make_state_bar_table(
    meta = microglia_subset_pilot_labeled@meta.data,
    group_col = config$genotype_col,
    state_col = config$state_col,
    group_levels = config$allowed_genotypes,
    state_levels = config$state_levels
  )
  
  p_sc <- make_state_stack_plot(
    bar_df = df_sc,
    group_col = config$genotype_col,
    state_col = config$state_col,
    colors = setNames(scvi_cols <- c(
      "#F8766D", "#C49A00", "#53B400", "#00C1A2",
      "#00A6FF", "#A58AFF", "#FB61D7", "#FF6F91"
    ), config$state_levels),
    title = "Pilot Single-Cell scVI state composition",
    xlab = "Pilot rs3732765 genotype",
    ylab = "Fraction of microglia"
  )
  save_plot_both(p_sc, "Pilot_SC_scANVI_stacked_bar", config$out_dir, width = 5, height = 4)
  
  # Single-nuclei
  df_nuc <- make_state_bar_table(
    meta = microglia_subset_nuclei_labeled@meta.data,
    group_col = config$genotype_col,
    state_col = config$state_col,
    group_levels = config$allowed_genotypes,
    state_levels = config$state_levels
  )
  
  p_nuc <- make_state_stack_plot(
    bar_df = df_nuc,
    group_col = config$genotype_col,
    state_col = config$state_col,
    colors = setNames(c(
      "#F8766D", "#C49A00", "#53B400", "#00C1A2",
      "#00A6FF", "#A58AFF", "#FB61D7", "#FF6F91"
    ), config$state_levels),
    title = "Pilot Single-Nuclei scVI state composition",
    xlab = "Pilot rs3732765 genotype",
    ylab = "Fraction of microglia"
  )
  save_plot_both(p_nuc, "Pilot_NUC_scANVI_stacked_bar", config$out_dir, width = 5, height = 4)
  
  # Combined SC + NUC
  combined_meta <- bind_rows(
    as.data.frame(microglia_subset_pilot_labeled@meta.data),
    as.data.frame(microglia_subset_nuclei_labeled@meta.data)
  )
  df_combined <- make_state_bar_table(
    meta = combined_meta,
    group_col = config$genotype_col,
    state_col = config$state_col,
    group_levels = config$allowed_genotypes,
    state_levels = config$state_levels
  )
  
  p_combined <- make_state_stack_plot(
    bar_df = df_combined,
    group_col = config$genotype_col,
    state_col = config$state_col,
    colors = setNames(c(
      "#F8766D", "#C49A00", "#53B400", "#00C1A2",
      "#00A6FF", "#A58AFF", "#FB61D7", "#FF6F91"
    ), config$state_levels),
    title = "Pilot combined SC + NUC scVI state composition",
    xlab = "Pilot rs3732765 genotype",
    ylab = "Fraction of microglia"
  )
  save_plot_both(p_combined, "Pilot_SC_NUC_combined_scANVI_stacked_bar", config$out_dir, width = 5, height = 4)
  
  # Save percentage tables for downstream use
  write_csv(df_sc, file.path(config$out_dir, "Pilot_SC_scANVI_percentages.csv"))
  write_csv(df_nuc, file.path(config$out_dir, "Pilot_NUC_scANVI_percentages.csv"))
  write_csv(df_combined, file.path(config$out_dir, "Pilot_SC_NUC_combined_scANVI_percentages.csv"))
  
  message("Saved scVI composition tables and plots.")
} else {
  message("Skipping scVI label attachment / stacked-bar plots: scVI objects not found in the environment.")
}

############################################################
## 10. Final outputs
############################################################

message("\nScript 1 complete.")
message("Outputs written to: ", normalizePath(config$out_dir))
message("  - merged Seurat object")
message("  - MTX/TSV export bundle")
message("  - pseudobulk matrix")
message("  - edgeR results")
message("  - DREAM results")
message("  - optional scVI-labeled Seurat objects and stacked-bar plots")
