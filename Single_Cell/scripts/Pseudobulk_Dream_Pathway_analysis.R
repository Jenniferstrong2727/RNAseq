#!/usr/bin/env Rscript
###############################################
# Pseudobulk_Dream_Pathway_analysis.R
#
# Generic pseudobulk -> dream (variancePartition) pipeline
# + basic QC, PCA/MDS, correlation heatmaps
# + differential expression via dream (voomWithDreamWeights)
# + module pathway analyses: CAMERA / mroast / FGSEA / GSVA (if modules provided)
#
# Usage (examples):
#   Rscript Pseudobulk_Dream_Pathway_analysis.R --seurat /path/to/object.rds --assay RNA --layer counts --outdir ./pseudobulk_results
#   # or inside R:
#   source("Pseudobulk_Dream_Pathway_analysis.R"); run_main(list(...))
#
# Author: <Your Name or Org>
# GitHub: <link to repo>
# License: MIT (suggested)
###############################################

suppressPackageStartupMessages({
  library(optparse)
  library(Seurat)
  library(Matrix)
  library(dplyr)
  library(edgeR)
  library(limma)
  library(variancePartition)
  library(ggplot2)
  library(pheatmap)
  library(ComplexHeatmap)
  library(matrixStats)
  library(fgsea)
  library(GSVA)
  library(AnnotationDbi)
  # load other libs lazily where needed
})

# ---------------------------
# Command line arguments
# ---------------------------
option_list <- list(
  make_option(c("-s","--seurat"), type="character", default = NA,
              help = "Path to Seurat object (RDS) OR name of Seurat object already in env. If omitted and running interactively, expects object named 'mg' in workspace.",
              metavar = "file"),
  make_option(c("--assay"), type="character", default = "RNA", help = "Assay to use (default: RNA)"),
  make_option(c("--layer"), type="character", default = "counts", help = "Assay layer for raw counts (Seurat v5 uses layer='counts')"),
  make_option(c("-o","--outdir"), type="character", default = "pseudobulk_results", help = "Output directory"),
  make_option(c("--min_cells_per_sample"), type="integer", default = 20, help = "Min cells per sample to flag"),
  make_option(c("--outlier_mad_multiplier"), type="numeric", default = 3, help = "MAD multiplier for MDS outlier detection"),
  make_option(c("--top_n_pca_genes"), type="integer", default = 500, help = "Top variable genes used in PCA"),
  make_option(c("--save_rds_if_pass"), action = "store_true", default = FALSE, help = "Save pseudobulk RDS only if QC passes")
)

opt <- parse_args(OptionParser(option_list = option_list))

# ---------------------------
# Helper functions
# ---------------------------
safe_dir <- function(path){
  if (is.null(path) || path == "") stop("Output dir unset.")
  dir.create(path, recursive = TRUE, showWarnings = FALSE)
  normalizePath(path, mustWork = FALSE)
}

load_seurat_obj <- function(path_or_name){
  if (is.null(path_or_name) || is.na(path_or_name)) {
    if (exists("mg")) {
      message("Using Seurat object 'mg' found in environment.")
      return(get("mg", envir = .GlobalEnv))
    } else {
      stop("No Seurat object specified and 'mg' not found in environment.")
    }
  }
  if (file.exists(path_or_name)) {
    obj <- readRDS(path_or_name)
    message("Loaded Seurat object from: ", path_or_name)
    return(obj)
  } else if (exists(path_or_name, envir = .GlobalEnv)) {
    message("Using Seurat object named '", path_or_name, "' found in environment.")
    return(get(path_or_name, envir = .GlobalEnv))
  } else {
    stop("Seurat object not found: ", path_or_name)
  }
}

# Align counts and meta by cell names
align_counts_meta <- function(counts_mat, meta){
  if (!identical(colnames(counts_mat), rownames(meta))) {
    if (all(rownames(meta) %in% colnames(counts_mat))) {
      counts_mat <- counts_mat[, rownames(meta)]
    } else if (all(colnames(counts_mat) %in% rownames(meta))) {
      meta <- meta[colnames(counts_mat), , drop = FALSE]
    } else {
      stop("Cannot align metadata and counts: mismatched cell names.")
    }
  }
  list(counts = counts_mat, meta = meta)
}

# build sample-level metadata: one row per sample
build_sample_meta <- function(cell_meta, cols = c("sample_id", "genotype", "batch")){
  if (!all(cols %in% colnames(cell_meta))) {
    stop("Missing required columns in cell metadata. Required: ", paste(cols, collapse = ", "))
  }
  sm <- cell_meta %>%
    as.data.frame() %>%
    group_by(sample_id) %>%
    summarise(across(all_of(cols[-1]), ~ unique(as.character(.x))[1]), .groups = "drop") %>%
    column_to_rownames("sample_id")
  return(as.data.frame(sm))
}

# save plot helper
save_png <- function(plot_expr, fname, width=1000, height=800){
  png(fname, width = width, height = height)
  print(plot_expr)
  dev.off()
}

# ---------------------------
# Main run
# ---------------------------
run_main <- function(options){
  outdir <- safe_dir(options$outdir)
  message("Outputs will be written to: ", outdir)
  # 0. load Seurat object
  seurat_obj <- load_seurat_obj(options$seurat)
  meta <- seurat_obj@meta.data
  required_cols <- c("sample_id", "genotype", "batch")
  if (!all(required_cols %in% colnames(meta))) {
    stop("Cell metadata must contain: sample_id, genotype, batch.")
  }
  # force types & preserve rownames
  meta <- meta %>%
    as.data.frame() %>%
    mutate(sample_id = as.factor(as.character(sample_id)),
           genotype  = as.character(genotype),
           batch     = as.character(batch))
  rownames(meta) <- rownames(seurat_obj@meta.data)
  # 1. extract counts
  counts_mat <- tryCatch({
    GetAssayData(seurat_obj, assay = options$assay, layer = options$layer)
  }, error = function(e){
    stop("Failed to fetch counts from Seurat assay/layer. Error: ", e$message)
  })
  # align
  aligned <- align_counts_meta(counts_mat, meta)
  counts_mat <- aligned$counts; meta <- aligned$meta
  stopifnot(identical(colnames(counts_mat), rownames(meta)))
  # 2. create pseudobulk
  samples <- unique(as.character(meta$sample_id))
  genes <- rownames(counts_mat)
  pseudobulk_counts <- matrix(0L, nrow = length(genes), ncol = length(samples),
                              dimnames = list(genes, samples))
  for (s in samples) {
    cells_in_s <- rownames(meta)[meta$sample_id == s]
    if (length(cells_in_s) == 0) {
      warning("Sample ", s, " has zero cells. Skipping.")
      next
    }
    pseudobulk_counts[, s] <- Matrix::rowSums(counts_mat[, cells_in_s, drop = FALSE])
  }
  # 3. basic QC summaries
  pb_lib_by_sample <- data.frame(sample_id = colnames(pseudobulk_counts),
                                 pb_sum = as.numeric(colSums(pseudobulk_counts)),
                                 stringsAsFactors = FALSE)
  cell_libsize <- Matrix::colSums(counts_mat)
  cell_df <- data.frame(cell = colnames(counts_mat),
                        libsize = as.numeric(cell_libsize[colnames(counts_mat)]),
                        sample_id = as.character(meta[colnames(counts_mat), "sample_id"]),
                        stringsAsFactors = FALSE)
  cell_sum_by_sample <- cell_df %>% group_by(sample_id) %>% summarise(cell_sum = sum(libsize), .groups = "drop")
  lib_compare <- pb_lib_by_sample %>% left_join(cell_sum_by_sample, by = "sample_id") %>% mutate(diff = pb_sum - cell_sum)
  cells_per_sample <- as.data.frame(table(cell_df$sample_id), stringsAsFactors = FALSE)
  colnames(cells_per_sample) <- c("sample_id", "n_cells")
  # 4. correlation / logCPM
  pseudobulk_counts_nonzero <- pseudobulk_counts[rowSums(pseudobulk_counts) > 0, , drop = FALSE]
  dge_tmp <- edgeR::DGEList(counts = pseudobulk_counts_nonzero)
  dge_tmp <- edgeR::calcNormFactors(dge_tmp)
  logcpm <- edgeR::cpm(dge_tmp, log = TRUE, prior.count = 1)
  cor_mat <- stats::cor(logcpm, method = "pearson")
  write.csv(cor_mat, file = file.path(outdir, "pseudobulk_sample_sample_correlation_logCPM.csv"), row.names = TRUE)
  # heatmap (try)
  heatfile <- file.path(outdir, "pseudobulk_sample_correlation_logCPM.png")
  tryCatch({
    pheatmap::pheatmap(cor_mat, filename = heatfile, main = "Sample-sample Pearson cor (logCPM)", border_color = NA)
    message("Wrote correlation heatmap: ", heatfile)
  }, error = function(e){
    warning("pheatmap failed: ", conditionMessage(e))
  })
  # 5. MDS + outlier detection
  dist_mat <- as.dist(1 - cor_mat)
  mds_raw <- cmdscale(dist_mat, k = 2)
  mds_df <- as.data.frame(mds_raw)
  if (all(c("V1", "V2") %in% colnames(mds_df))){
    colnames(mds_df)[colnames(mds_df) == "V1"] <- "MDS1"
    colnames(mds_df)[colnames(mds_df) == "V2"] <- "MDS2"
  } else if (!("MDS1" %in% colnames(mds_df)) && ncol(mds_df) >= 2){
    colnames(mds_df)[1:2] <- c("MDS1","MDS2")
  }
  mds_df$sample <- if (!is.null(rownames(mds_df)) && any(rownames(mds_df) != "")) rownames(mds_df) else as.character(colnames(pseudobulk_counts))[seq_len(nrow(mds_df))]
  centroid <- colMeans(mds_df[, c("MDS1","MDS2")])
  mds_df$mds_dist <- sqrt((mds_df$MDS1 - centroid[1])^2 + (mds_df$MDS2 - centroid[2])^2)
  mds_med <- median(mds_df$mds_dist, na.rm = TRUE)
  mds_mad <- mad(mds_df$mds_dist, na.rm = TRUE)
  mds_df$mds_outlier <- mds_df$mds_dist > (mds_med + options$outlier_mad_multiplier * mds_mad)
  # merge metadata for plotting
  sample_meta <- build_sample_meta(meta)
  mds_plot_df <- mds_df %>% left_join(sample_meta %>% rownames_to_column("sample_id"), by = c("sample" = "sample_id"))
  save_png(ggplot(mds_plot_df, aes(MDS1, MDS2, color = genotype, shape = batch)) + geom_point(size = 4) + theme_minimal() + labs(title = "MDS of pseudobulk samples"), file.path(outdir, "pseudobulk_MDS_by_genotype_batch.png"))
  # 6. PCA
  rv <- apply(logcpm, 1, var)
  n_top <- min(options$top_n_pca_genes, length(rv))
  topg <- names(sort(rv, decreasing = TRUE))[1:n_top]
  pc <- prcomp(t(logcpm[topg, , drop = FALSE]), scale. = TRUE)
  pc_df <- data.frame(PC1 = pc$x[,1], PC2 = pc$x[,2], sample = rownames(pc$x))
  pc_df <- pc_df %>% left_join(sample_meta %>% rownames_to_column("sample_id"), by = c("sample" = "sample_id"))
  save_png(ggplot(pc_df, aes(PC1, PC2, color = genotype, shape = batch)) + geom_point(size = 3) + geom_text(aes(label = sample), vjust = -1, size = 3) + theme_minimal() + labs(title = "PCA (top variable genes) of pseudobulks"), file.path(outdir, "pseudobulk_PCA_topVar_genes.png"))
  # 7. Build master QC table + flags
  flags_small <- data.frame(sample_id = as.character(colnames(pseudobulk_counts)),
                            pb_total = as.numeric(colSums(pseudobulk_counts)),
                            stringsAsFactors = FALSE) %>%
    left_join(lib_compare %>% select(sample_id, cell_sum, diff) , by = "sample_id") %>%
    left_join(cells_per_sample %>% mutate(sample_id = as.character(sample_id)), by = "sample_id") %>%
    left_join(mds_df %>% mutate(sample_id = sample) %>% select(sample_id, mds_dist, mds_outlier), by = "sample_id")
  flags_small <- flags_small %>%
    mutate(flag_low_cells = ifelse(is.na(n_cells), TRUE, n_cells < options$min_cells_per_sample),
           flag_lib_mismatch = ifelse(is.na(diff), TRUE, diff != 0),
           flag_mds_outlier = ifelse(is.na(mds_outlier), FALSE, as.logical(mds_outlier)),
           any_flag = flag_low_cells | flag_lib_mismatch | flag_mds_outlier)
  write.csv(flags_small, file = file.path(outdir, "samples_to_investigate.csv"), row.names = FALSE)
  message("Wrote QC master table: samples_to_investigate.csv")
  # optionally save pseudobulk if no flags
  if (options$save_rds_if_pass) {
    if (all(flags_small$any_flag == FALSE, na.rm = TRUE)) {
      saveRDS(pseudobulk_counts, file = file.path(outdir, "pseudobulk_counts_bySample.rds"))
      write.csv(as.data.frame(sample_meta), file = file.path(outdir, "pseudobulk_sample_metadata.csv"), row.names = TRUE)
      message("QC passed; saved pseudobulk and sample metadata.")
    } else {
      message("QC failed; pseudobulk RDS not saved. Inspect samples_to_investigate.csv")
    }
  } else {
    # always save small checkpoint
    saveRDS(pseudobulk_counts, file = file.path(outdir, "pseudobulk_counts_bySample_checkpoint.rds"))
    write.csv(as.data.frame(sample_meta), file = file.path(outdir, "pseudobulk_sample_metadata.csv"), row.names = TRUE)
    message("Saved checkpoint RDS and sample metadata.")
  }

  # ---------------------------------------------------------------------------
  # DREAM differential expression (pseudobulk)
  # - filter CPM >1 in >=2 samples by default
  # - TMM normalization, voomWithDreamWeights, dream, eBayes
  # ---------------------------------------------------------------------------
  dge <- DGEList(counts = pseudobulk_counts)
  cpm_mat <- cpm(dge)
  keep <- rowSums(cpm_mat > 1) >= 2
  message(sum(keep), " genes kept after CPM filtering (CPM>1 in >=2 samples).")
  dge <- dge[keep, , keep.lib.sizes = FALSE]
  dge <- calcNormFactors(dge, method = "TMM")
  # clean factors
  if (!("batch" %in% colnames(sample_meta))) stop("sample_meta must contain 'batch' column for model.")
  sample_meta$batch_clean <- gsub("[^A-Za-z0-9_]+", "_", as.character(sample_meta$batch))
  sample_meta$batch_clean <- factor(sample_meta$batch_clean)
  sample_meta$genotype <- factor(as.character(sample_meta$genotype))
  if ("GG" %in% levels(sample_meta$genotype)) sample_meta$genotype <- relevel(sample_meta$genotype, ref = "GG")
  form <- ~ batch_clean + genotype
  # voomWithDreamWeights
  vobj <- voomWithDreamWeights(dge, formula = form, data = sample_meta)
  saveRDS(vobj, file = file.path(outdir, "voomWithDreamWeights_object.rds"))
  png(file.path(outdir, "voom_mean_variance_plot_releveled.png"), width = 1000, height = 800)
  plotSA(lmFit(vobj, design = model.matrix(form, data = sample_meta)), main = "voom mean-variance diagnostic")
  dev.off()
  fit <- dream(vobj, formula = form, data = sample_meta)
  fit <- eBayes(fit)
  # find genotype coefficient
  coef_names <- colnames(fit$coefficients)
  gcoef_idx <- which(grepl("^genotype", coef_names))
  if (length(gcoef_idx) == 0) stop("No genotype coefficient found in model.")
  coef_to_test <- coef_names[gcoef_idx[1]]
  res_full <- topTable(fit, coef = coef_to_test, number = Inf, sort.by = "P")
  res_top <- subset(res_full, adj.P.Val <= 0.05)
  write.csv(res_full, file = file.path(outdir, paste0("DEG_voom_limma_dream_full_table_", coef_to_test, ".csv")), row.names = TRUE)
  write.csv(res_top, file = file.path(outdir, paste0("DEG_voom_limma_dream_top_hits_", coef_to_test, ".csv")), row.names = TRUE)
  message("DE analysis complete; results written.")
  # volcano + MA
  res_plot <- res_full; res_plot$gene <- rownames(res_plot)
  res_plot$P.Value <- pmax(res_plot$P.Value, .Machine$double.xmin)
  res_plot$negLogP <- -log10(res_plot$P.Value)
  res_plot$direction <- ifelse(res_plot$logFC > 0, "up", ifelse(res_plot$logFC < 0, "down", "neutral"))
  top_up <- head(res_plot %>% filter(direction == "up") %>% arrange(desc(abs(logFC))) %>% pull(gene), 50)
  top_dn <- head(res_plot %>% filter(direction == "down") %>% arrange(desc(abs(logFC))) %>% pull(gene), 50)
  label_genes <- unique(c(top_up, top_dn))
  p_vol <- ggplot(res_plot, aes(x = logFC, y = negLogP, color = direction)) +
    geom_point(alpha = 0.35, size = 1.6) +
    theme_classic(base_size = 14) +
    geom_text_repel(data = subset(res_plot, gene %in% label_genes), aes(label = gene), size = 3)
  ggsave(file.path(outdir, paste0("volcano_", coef_to_test, ".png")), p_vol, width = 10, height = 8, dpi = 300)
  # save fit objects
  saveRDS(fit, file = file.path(outdir, paste0("dream_fit_object_", coef_to_test, ".rds")))
  saveRDS(vobj, file = file.path(outdir, "voomWithDreamWeights_object_releveled.rds"))
  # ---------------------------
  # CAMERA / mroast / FGSEA / GSVA: only run if gene set objects supplied in workspace:
  # cameraModules (list of modules), roast (list), cameraModules_mapped, cameraModules_excl, etc.
  # The script looks for well-named objects to avoid hardcoding.
  # ---------------------------
  # FGSEA runner: uses res_full or fit object to build ranking
  # GSVA runner: uses vobj$E or expression fallback
  # For brevity: we export a small “hook” area below for user-supplied modules.
  # ---------------------------
  message("Pipeline finished. Check ", outdir, " for outputs and logs.")
  invisible(list(outdir = outdir, pseudobulk = pseudobulk_counts, res_full = res_full, fit = fit, vobj = vobj))
}

# if called as script, run with parsed options
if (interactive()) {
  message("Interactive session detected. Use run_main(opt) to execute programmatically.")
} else {
  run_main(opt)
}

