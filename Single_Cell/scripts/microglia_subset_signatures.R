###############################################################################
# microglia_pipeline_cleaned.R
#
# Clean, modular, and generic pipeline to:
#  - identify microglia from score/module columns
#  - subset microglia and re-process (SCT PCA / UMAP / clustering)
#  - export markers, average expression, cluster abundance & genotype stats
#  - compute module scores (LG + SM), z-scores, and cluster-level assignment
#  - produce per-gene FeaturePlot panels and ComplexHeatmap outputs
#  - persist heatmap matrices and module info into Seurat@misc and save RDS
#
# USAGE:
#  - edit the "USER SETTINGS" section below (input file, object names, output dir)
#  - run interactively or in a batch R session
#
# DEPENDENCIES:
#  Seurat (v4/5), dplyr, tidyr, readr, ComplexHeatmap, circlize, cowplot,
#  patchwork, ggplot2, scales, ggrepel, DescTools (optional for Cramer's V)
#
###############################################################################

########################
# USER SETTINGS (edit)
########################
input_rds            <- "integrated_NEW_MARKERS.rds"   # input integrated object
seurat_obj_varname   <- "obj"                          # variable name to assign the read object
microglia_rds_name   <- "Microglia_gene_signatures.rds"
outdir               <- "/Users/jennifer.strong/Documents/Rworkingdirectory/P2RY12/sc_P2RY12/MICROGLIA_subset"
save_rds_compress    <- "xz"
min_pct_findmarkers  <- 0.10
logfc_thresh_find    <- 0.25
hg_ctrl_for_module   <- 100
z_assignment_thresh  <- 0.25
n_top_modules_export <- 3
cluster_dims         <- 1:12
umap_min_dist        <- 0.02
umap_spread          <- 1.5
neighbors_k          <- 20
clustering_resolution <- 3.0
clustering_algorithm  <- 4    # Leiden

# Optional: curated filter genes for downstream exports
filter_genes <- unique(c(
  "P2RY12","TMEM119","SALL1","CX3CR1","CSF1R","APOE","TREM2","SPP1","IFITM3","ISG15"
  # add/modify as desired
))

########################
# LIBRARIES
########################
suppressMessages({
  library(Seurat)
  library(dplyr)
  library(tidyr)
  library(readr)
  library(ggplot2)
  library(cowplot)
  library(patchwork)
  library(scales)
  library(ComplexHeatmap)
  library(circlize)
  library(grid)
  library(ggrepel)
  library(tibble)
})

dir.create(outdir, recursive = TRUE, showWarnings = FALSE)

########################
# UTILITIES
########################
log_msg <- function(...) message("[", Sys.time(), "] ", paste0(...))

per_gene_z <- function(mat) {
  z <- t(scale(t(as.matrix(mat)), center = TRUE, scale = TRUE))
  z[is.na(z)] <- 0
  z
}

ensure_meta_col <- function(obj, col) {
  if (!col %in% colnames(obj@meta.data)) {
    stop("Required metadata column not found: ", col)
  }
}

########################
# 0. LOAD OBJECT (generic)
########################
log_msg("Loading input RDS: ", input_rds)
stopifnot(file.exists(input_rds))
obj <- readRDS(input_rds)
# optional reassign so script refers to generic name 'obj'
assign(seurat_obj_varname, obj)
rm(obj)
obj <- get(seurat_obj_varname)
log_msg("Object loaded. DefaultAssay: ", DefaultAssay(obj))

########################
# 1. CONFIRM MODULE/SCORE COLUMN PRESENCE
########################
meta <- obj@meta.data
scoreV2_cols <- grep("^scoreV2_", colnames(meta), value = TRUE)
log_msg("Found scoreV2 columns: ", paste(scoreV2_cols, collapse = ", "))

if (length(scoreV2_cols) == 0) {
  stop("No scoreV2_* columns found in object metadata. Please compute module scores first.")
}

# create per-cell 'top from scores' if missing
if (!"moduleV2_top_from_scores" %in% colnames(obj@meta.data)) {
  log_msg("Computing moduleV2_top_from_scores (per-cell top score).")
  scores_v2 <- as.matrix(obj@meta.data[, scoreV2_cols, drop = FALSE])
  top_from_scores <- apply(scores_v2, 1, function(x) {
    nm <- names(x)[which.max(x)]
    sub("^scoreV2_", "", nm)
  })
  obj$moduleV2_top_from_scores <- factor(top_from_scores)
}

# optional existing moduleV2_celltype check
if ("moduleV2_celltype" %in% colnames(obj@meta.data)) {
  log_msg("Existing moduleV2_celltype found (will be used in inclusive microglia call).")
} else {
  log_msg("No moduleV2_celltype column found; inclusive microglia call will rely on scores only.")
}

########################
# 2. DERIVE INCLUSIVE MICROGLIA CALL (generic)
########################
# identify microglia score column dynamically (case-insensitive)
micro_col <- grep("Microglia", scoreV2_cols, ignore.case = TRUE, value = TRUE)
if (length(micro_col) != 1) {
  # if ambiguous, try exact match 'Microglia' after strip
  micro_col <- scoreV2_cols[grepl("microglia", scoreV2_cols, ignore.case = TRUE)][1]
  if (is.na(micro_col)) stop("Unable to identify a single scoreV2_*Microglia* column.")
}
micro_label <- sub("^scoreV2_", "", micro_col)

# inclusive call: top-from-scores OR moduleV2_celltype == micro_label (if present)
cells_top_micro <- colnames(obj)[obj$moduleV2_top_from_scores == micro_label]
cells_moduleV2_micro <- character(0)
if ("moduleV2_celltype" %in% colnames(obj@meta.data)) {
  cells_moduleV2_micro <- rownames(obj@meta.data)[obj@meta.data$moduleV2_celltype == micro_label]
}
micro_cells_inclusive <- union(cells_top_micro, cells_moduleV2_micro)
log_msg("Inclusive microglia cells found: ", length(micro_cells_inclusive))

obj$microgliaV2_inclusive_call <- factor(
  ifelse(colnames(obj) %in% micro_cells_inclusive, "Microglia_inclusive", "NonMicroglia"),
  levels = c("Microglia_inclusive", "NonMicroglia")
)

########################
# 3. SUBSET MICROGLIA (inclusive) & SAVE INTERMEDIATE RDS
########################
micro_cells <- colnames(obj)[obj$microgliaV2_inclusive_call == "Microglia_inclusive"]
log_msg("Subsetting to microglia (n = ", length(micro_cells), ")")
mg <- subset(obj, cells = micro_cells)

# optional: save quick subset
saveRDS(mg, file = file.path(outdir, "microglia_subset_v2_inclusive_percell.rds"))
log_msg("Saved microglia subset RDS (inclusive).")

########################
# 4. REPROCESS MICROGLIA (SCT PCA → UMAP → neighbors → clusters)
########################
DefaultAssay(mg) <- "SCT"
# remove existing scale data (safe re-scaling)
mg[["SCT"]]@scale.data <- matrix()
VariableFeatures(mg) <- NULL

log_msg("Finding HVG (nfeatures=3000) on SCT...")
mg <- FindVariableFeatures(mg, assay = "SCT", selection.method = "vst", nfeatures = 3000, verbose = FALSE)

log_msg("Scaling SCT (variable features)...")
mg <- ScaleData(mg, assay = "SCT", features = VariableFeatures(mg), verbose = FALSE)

log_msg("Running PCA (npcs=50)...")
mg <- RunPCA(mg, assay = "SCT", features = VariableFeatures(mg), npcs = 50, verbose = FALSE)

# remove old harmony if present
if ("harmony" %in% names(mg@reductions)) mg[["harmony"]] <- NULL

log_msg("Running UMAP, neighbors, clustering...")
mg <- RunUMAP(mg, reduction = "pca", dims = cluster_dims, min.dist = umap_min_dist, spread = umap_spread, verbose = FALSE)
mg <- FindNeighbors(mg, reduction = "pca", dims = cluster_dims, k.param = neighbors_k, verbose = FALSE)
mg <- FindClusters(mg, resolution = clustering_resolution, algorithm = clustering_algorithm, verbose = FALSE)
mg$mg_cluster <- Idents(mg)

# save reprocessed mg object
fn_reproc <- file.path(outdir, "microglia_subset_v2_inclusive_SCT_reprocessed.rds")
saveRDS(mg, fn_reproc)
log_msg("Saved reprocessed microglia object: ", fn_reproc)

########################
# 5. PLOTS: UMAP + FeaturePanel + QC
########################
# UMAP cluster plot
p_umap <- DimPlot(mg, reduction = "umap", group.by = "mg_cluster", label = TRUE, repel = TRUE) +
  theme_void(base_size = 14) + theme(legend.position = "right")
ggsave(file.path(outdir, "UMAP_microglia_mg_clusters.png"), p_umap, width = 7, height = 6, dpi = 600)
ggsave(file.path(outdir, "UMAP_microglia_mg_clusters.pdf"), p_umap, width = 7, height = 6, device = cairo_pdf)

# small FeaturePlot panel helper (per-gene scaling per-gene)
make_feature_panel <- function(object, genes, out_prefix, pt.size = 0.6, ncol = 3) {
  present <- genes %in% rownames(object)
  if (any(!present)) {
    warning("Missing genes (skipping): ", paste(genes[!present], collapse = ", "))
    genes <- genes[present]
    if (length(genes) == 0) return(NULL)
  }
  expr_all <- as.matrix(GetAssayData(object, layer = "data")[genes, , drop = FALSE])
  gene_limits <- vapply(genes, function(g) {
    mat <- expr_all[g, ]
    up <- as.numeric(quantile(mat, probs = 0.995, na.rm = TRUE))
    if (!is.finite(up) || up <= 0) up <- max(mat, na.rm = TRUE)
    if (!is.finite(up) || up <= 0) up <- 1
    up
  }, numeric(1))
  make_feat <- function(g) {
    up <- gene_limits[g]
    p <- FeaturePlot(object, features = g, reduction = "umap", order = TRUE, pt.size = pt.size)
    if (inherits(p, "list")) p <- p[[1]]
    p <- p + scale_colour_gradient(low = "lightgrey", high = "#2b83ba", limits = c(0, up), oob = scales::squish) +
      theme_void() + ggtitle(g) + theme(plot.title = element_text(hjust = 0.5, size = 10), legend.position = "right")
    p
  }
  plots <- lapply(genes, make_feat)
  combined <- wrap_plots(plotlist = plots, ncol = ncol)
  ggsave(paste0(out_prefix, ".png"), combined, width = 10, height = 9, dpi = 300)
  ggsave(paste0(out_prefix, ".pdf"), combined, width = 10, height = 9, device = cairo_pdf)
  combined
}

genes_to_plot <- c("P2RY12","CX3CR1","TMEM119","APOE","HLA-DRA","IFITM3","S100A11")
make_feature_panel(mg, genes_to_plot, file.path(outdir, "FeaturePlots_all_markers"))

# QC violin per cluster
qc_metrics <- c("nFeature_RNA", "nCount_RNA", "percent.mt")
plist_qc <- lapply(qc_metrics, function(feat) {
  VlnPlot(mg, features = feat, group.by = "mg_cluster", pt.size = 0.5) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
    labs(x = "Microglia cluster")
})
p_qc <- plot_grid(plotlist = plist_qc, ncol = 1, align = "v")
ggsave(file.path(outdir, "microglia_QC_violin_by_cluster.png"), p_qc, width = 9, height = 9, dpi = 400)

########################
# 6. FINDALLMARKERS + AVERAGE EXPRESSION (RNA, slot='data')
########################
DefaultAssay(mg) <- "RNA"
Idents(mg) <- mg$mg_cluster

log_msg("Running FindAllMarkers() on RNA (slot = 'data') ...")
markers_all <- FindAllMarkers(object = mg, assay = "RNA", slot = "data",
                              test.use = "wilcox", only.pos = TRUE,
                              min.pct = min_pct_findmarkers, logfc.threshold = logfc_thresh_find)

fn_full_markers <- file.path(outdir, "microglia_FindAllMarkers_full_byCluster_logFCdesc.csv")
write.csv(markers_all %>% arrange(cluster, desc(avg_log2FC)), fn_full_markers, row.names = FALSE)
log_msg("Wrote: ", fn_full_markers)

# filtered markers to curated list (if present)
genes_present <- filter_genes[filter_genes %in% rownames(mg)]
if (length(setdiff(filter_genes, genes_present)) > 0) {
  write.csv(data.frame(missing_genes = setdiff(filter_genes, genes_present)),
            file.path(outdir, "filter_genes_missing_from_object.csv"), row.names = FALSE)
  log_msg("Some filter_genes missing: wrote missing list.")
}
markers_filtered <- markers_all %>% filter(gene %in% genes_present)
fn_filtered_markers <- file.path(outdir, "microglia_FindAllMarkers_filtered_byCluster_logFCdesc.csv")
write.csv(markers_filtered, fn_filtered_markers, row.names = FALSE)
log_msg("Wrote: ", fn_filtered_markers)

# AverageExpression (wide matrix)
log_msg("Computing AverageExpression() ...")
avg_list <- AverageExpression(object = mg, assays = "RNA", slot = "data", group.by = "mg_cluster", verbose = FALSE)
avg_mat <- avg_list$RNA
# reorder columns by numeric cluster if possible
cluster_cols <- colnames(avg_mat)
cluster_nums <- as.numeric(gsub("[^0-9-]", "", cluster_cols))
if (!any(is.na(cluster_nums))) {
  ord_idx <- order(cluster_nums, na.last = TRUE)
  avg_mat <- avg_mat[, ord_idx, drop = FALSE]
}
fn_avg_full <- file.path(outdir, "microglia_AverageExpression_full_genes_byCluster_ordered.csv")
write.csv(as.data.frame(avg_mat) %>% rownames_to_column("gene"), fn_avg_full, row.names = FALSE)
log_msg("Wrote: ", fn_avg_full)

# filtered average expression (curated genes)
if (length(genes_present) > 0) {
  avg_mat_filtered <- avg_mat[rownames(avg_mat) %in% genes_present, , drop = FALSE]
  fn_avg_filtered <- file.path(outdir, "microglia_AverageExpression_filtered_genes_byCluster_ordered.csv")
  write.csv(as.data.frame(avg_mat_filtered) %>% rownames_to_column("gene"), fn_avg_filtered, row.names = FALSE)
  log_msg("Wrote: ", fn_avg_filtered)
}

# long version
avg_long_df <- as.data.frame(avg_mat) %>%
  rownames_to_column("gene") %>%
  pivot_longer(cols = -gene, names_to = "cluster_raw", values_to = "avg_expression") %>%
  mutate(cluster_num = as.numeric(gsub("[^0-9-]", "", cluster_raw))) %>%
  arrange(cluster_num, desc(avg_expression))
fn_avg_long <- file.path(outdir, "microglia_AverageExpression_full_long_ordered.csv")
write.csv(avg_long_df, fn_avg_long, row.names = FALSE)
log_msg("Wrote: ", fn_avg_long)

########################
# 7. MODULE SCORING (LG / SM), z-scores, cluster-mean z, and auto-assign labels
########################
# ---- define modules (example lists) ----
lg_modules <- list(
  homeostatic = c("P2RY12","TMEM119","SALL1","CX3CR1","CSF1R","GPR34","FCRLS"),
  DAM_MGnD    = c("TREM2","TYROBP","APOE","SPP1","GPNMB","LPL","CD9","ITGAX","C1QA","C1QB","C1QC"),
  IFN         = c("IFIT1","IFIT2","IFITM3","IFI27","ISG15","IRF7","STAT1","OASL","GBP2"),
  inflammatory= c("IL1B","IL6","TNF","NLRP3","CASP1","PYCARD","CXCL10","CCL2"),
  phagocytic  = c("MERTK","AXL","MFGE8","GULP1","CD68","ITGAM","ITGB2","LAMP1","LAMP2","CTSB","CTSD"),
  lysosomal   = c("GBA","TFEB","NPC1","NPC2","SCARB2","LAMP1","LAMP2","CTSB","CTSD"),
  lipid_LDAM  = c("PLIN2","CIDEC","G0S2","FABP5","LPL","APOE")
)

sm_modules <- list(
  homeostatic = c("P2RY12","TMEM119","SALL1","CX3CR1"),
  DAM         = c("TREM2","APOE","SPP1","GPNMB"),
  MGnD        = c("TREM2","APOE","TYROBP","C1QA"),
  IFN         = c("IFIT1","IFIT2","IFITM3","ISG15"),
  inflammatory= c("IL1B","TNF","NLRP3","IL6"),
  phagocytic  = c("MERTK","AXL","CD68","ITGAM"),
  lipid_LDAM  = c("PLIN2","LPL","APOE","CIDEC"),
  lysosomal   = c("LAMP1","CTSD","TFEB","GBA")
)

# ensure RNA assay data exist (normalized)
DefaultAssay(mg) <- "RNA"
if (is.null(mg@assays$RNA@data) || ncol(mg@assays$RNA@data) == 0) {
  stop("RNA assay data is empty. Run NormalizeData(mg, assay='RNA') before running module scoring.")
}

# AddModuleScore for LG and SM
log_msg("Adding LG module scores...")
mg <- AddModuleScore(mg, features = lg_modules, assay = "RNA", name = "lg_module", ctrl = hg_ctrl_for_module, search = TRUE)
lg_score_cols <- grep("^lg_module", colnames(mg@meta.data), value = TRUE)
if (length(lg_score_cols) == length(lg_modules)) {
  new_lg_names <- paste0("lg_module_", names(lg_modules))
  colnames(mg@meta.data)[match(lg_score_cols, colnames(mg@meta.data))] <- new_lg_names
  lg_score_cols <- new_lg_names
} else {
  warning("LG: unexpected number of AddModuleScore columns; leaving names as-is")
  lg_score_cols <- grep("^lg_module", colnames(mg@meta.data), value = TRUE)
}

log_msg("Adding SM module scores...")
mg <- AddModuleScore(mg, features = sm_modules, assay = "RNA", name = "sm_module", ctrl = hg_ctrl_for_module, search = TRUE)
sm_score_cols <- grep("^sm_module", colnames(mg@meta.data), value = TRUE)
if (length(sm_score_cols) == length(sm_modules)) {
  new_sm_names <- paste0("sm_module_", names(sm_modules))
  colnames(mg@meta.data)[match(sm_score_cols, colnames(mg@meta.data))] <- new_sm_names
  sm_score_cols <- new_sm_names
} else {
  warning("SM: unexpected number of AddModuleScore columns; leaving names as-is")
  sm_score_cols <- grep("^sm_module", colnames(mg@meta.data), value = TRUE)
}

# compute per-cell z-scores
make_z <- function(x) {
  if (all(is.na(x)) || sd(x, na.rm = TRUE) == 0) return(rep(0, length(x)))
  (x - mean(x, na.rm = TRUE)) / sd(x, na.rm = TRUE)
}
for (cname in lg_score_cols) mg@meta.data[[sub("^lg_module_", "z_lg_", cname)]] <- make_z(mg@meta.data[[cname]])
for (cname in sm_score_cols) mg@meta.data[[sub("^sm_module_", "z_sm_", cname)]] <- make_z(mg@meta.data[[cname]])

# cluster-level mean z-scores (LG & SM)
lg_z_cols <- grep("^z_lg_", colnames(mg@meta.data), value = TRUE)
sm_z_cols <- grep("^z_sm_", colnames(mg@meta.data), value = TRUE)

lg_cluster_means <- sapply(lg_z_cols, function(cn) tapply(mg@meta.data[[cn]], mg@meta.data$mg_cluster, mean, na.rm = TRUE))
lg_cluster_means <- as.data.frame(lg_cluster_means) %>% rownames_to_column("cluster")
write.csv(lg_cluster_means, file.path(outdir, "lg_cluster_mean_zscores.csv"), row.names = FALSE)
log_msg("Wrote LG cluster mean z-scores.")

sm_cluster_means <- sapply(sm_z_cols, function(cn) tapply(mg@meta.data[[cn]], mg@meta.data$mg_cluster, mean, na.rm = TRUE))
sm_cluster_means <- as.data.frame(sm_cluster_means) %>% rownames_to_column("cluster")
write.csv(sm_cluster_means, file.path(outdir, "sm_cluster_mean_zscores.csv"), row.names = FALSE)
log_msg("Wrote SM cluster mean z-scores.")

# assign per-cluster top module if mean z >= threshold
assign_top_module <- function(cluster_df, z_threshold) {
  df <- cluster_df
  module_cols <- colnames(df)[-1]
  df2 <- df %>%
    rowwise() %>%
    mutate(
      best_module = {
        vals <- c_across(all_of(module_cols))
        best_i <- which.max(vals)
        best_val <- vals[best_i]
        if (best_val >= z_threshold) module_cols[best_i] else "Unassigned"
      },
      best_val = {
        vals <- c_across(all_of(module_cols))
        vals[which.max(vals)]
      }
    ) %>% ungroup()
  df2 %>% select(cluster, best_module, best_val)
}

lg_labels <- assign_top_module(lg_cluster_means, z_assignment_thresh)
sm_labels <- assign_top_module(sm_cluster_means, z_assignment_thresh)
write.csv(lg_labels, file.path(outdir, paste0("lg_cluster_auto_labels_z", z_assignment_thresh, "_TOPONLY.csv")), row.names = FALSE)
write.csv(sm_labels, file.path(outdir, paste0("sm_cluster_auto_labels_z", z_assignment_thresh, "_TOPONLY.csv")), row.names = FALSE)
log_msg("Wrote LG/SM auto labels.")

# annotate mg@meta.data per-cell mapped auto label (cluster -> top module)
lg_label_map <- setNames(lg_labels$best_module, lg_labels$cluster)
sm_label_map <- setNames(sm_labels$best_module, sm_labels$cluster)

mg@meta.data$auto_label_lg_top_z <- ifelse(mg@meta.data$mg_cluster %in% names(lg_label_map), lg_label_map[as.character(mg@meta.data$mg_cluster)], "Unassigned")
mg@meta.data$auto_label_sm_top_z <- ifelse(mg@meta.data$mg_cluster %in% names(sm_label_map), sm_label_map[as.character(mg@meta.data$mg_cluster)], "Unassigned")

# export updated metadata
fn_meta_out <- file.path(outdir, "microglia_meta_with_module_scores_and_labels.csv")
write.csv(mg@meta.data, fn_meta_out, row.names = TRUE)
log_msg("Wrote updated metadata with scores & labels: ", fn_meta_out)

########################
# 8. EXPORT TOP-N MODULES PER CLUSTER (LG + SM)
########################
read_or_compute_cluster_means <- function(prefix, zprefix = NULL) {
  fn <- file.path(outdir, paste0(prefix, "_cluster_mean_zscores.csv"))
  if (file.exists(fn)) {
    df <- read_csv(fn, show_col_types = FALSE)
    if (!"cluster" %in% colnames(df)) colnames(df)[1] <- "cluster"
    df$cluster <- as.character(df$cluster)
    return(df)
  } else {
    if (exists("mg")) {
      if (is.null(zprefix)) zprefix <- paste0("z_", prefix, "_")
      zcols <- grep(paste0("^", zprefix), colnames(mg@meta.data), value = TRUE)
      if (length(zcols) == 0) stop("No z-score columns found for prefix: ", prefix)
      mat <- sapply(zcols, function(cn) tapply(mg@meta.data[[cn]], mg@meta.data$mg_cluster, mean, na.rm = TRUE))
      df <- as.data.frame(mat) %>% rownames_to_column("cluster")
      colnames(df)[-1] <- sub(paste0("^", zprefix), "", colnames(df)[-1])
      return(df)
    }
    stop("Neither CSV nor mg in memory available for prefix: ", prefix)
  }
}

top_n_positive_per_cluster <- function(cluster_df, n = 3, min_pos = 0) {
  module_cols <- colnames(cluster_df)[-1]
  out_rows <- list()
  for (i in seq_len(nrow(cluster_df))) {
    clust <- as.character(cluster_df$cluster[i])
    vals <- as.numeric(cluster_df[i, module_cols])
    names(vals) <- module_cols
    pos_idx <- which(vals > min_pos)
    if (length(pos_idx) == 0) {
      for (r in 1:n) out_rows[[length(out_rows) + 1]] <- tibble(cluster = clust, rank = r, module = NA_character_, mean_z = NA_real_)
    } else {
      ordered <- sort(vals[pos_idx], decreasing = TRUE)
      mods <- names(ordered)
      for (r in 1:n) {
        if (r <= length(mods)) out_rows[[length(out_rows) + 1]] <- tibble(cluster = clust, rank = r, module = mods[r], mean_z = as.numeric(ordered[r]))
        else out_rows[[length(out_rows) + 1]] <- tibble(cluster = clust, rank = r, module = NA_character_, mean_z = NA_real_)
      }
    }
  }
  bind_rows(out_rows)
}

lg_means_df <- read_or_compute_cluster_means("lg", zprefix = "z_lg_")
lg_top3 <- top_n_positive_per_cluster(lg_means_df, n = n_top_modules_export, min_pos = 0) %>% mutate(source = "lg")
write_csv(lg_top3, file.path(outdir, "lg_top3_modules_per_cluster.csv"))
sm_means_df <- read_or_compute_cluster_means("sm", zprefix = "z_sm_")
sm_top3 <- top_n_positive_per_cluster(sm_means_df, n = n_top_modules_export, min_pos = 0) %>% mutate(source = "sm")
write_csv(sm_top3, file.path(outdir, "sm_top3_modules_per_cluster.csv"))
combined_top <- bind_rows(lg_top3, sm_top3) %>% arrange(source, as.numeric(gsub("[^0-9-]", "", cluster)), rank)
write_csv(combined_top, file.path(outdir, "all_top3_modules_per_cluster_combined.csv"))
log_msg("Wrote top-N module CSVs.")

########################
# 9. UMAP SPLIT BY GENOTYPE & ABUNDANCE TABLES
########################
# ensure cluster factor ordering (numeric preferred)
orig_levels <- sort(unique(as.character(mg$mg_cluster)))
suppressWarnings({ numeric_levels <- as.numeric(orig_levels) })
if (!any(is.na(numeric_levels))) {
  new_levels <- as.character(sort(numeric_levels))
} else new_levels <- sort(orig_levels)
mg$mg_cluster <- factor(as.character(mg$mg_cluster), levels = new_levels)

# split UMAP
p_split <- DimPlot(mg, reduction = "umap", group.by = "mg_cluster", split.by = "genotype", label = FALSE, repel = TRUE, pt.size = 0.65) +
  theme_void(base_size = 14) + theme(strip.text = element_text(size = 12))
p_legend <- DimPlot(mg, reduction = "umap", group.by = "mg_cluster", label = FALSE, pt.size = 3) + theme_void(base_size = 12) +
  guides(color = guide_legend(title = "mg_cluster", override.aes = list(size = 4))) + theme(legend.text = element_text(size = 9))
legend_grob <- cowplot::get_legend(p_legend)
combined_umap <- cowplot::plot_grid(p_split + theme(legend.position = "none"), legend_grob, ncol = 2, rel_widths = c(1,0.28))
ggsave(file.path(outdir, "UMAP_microglia_split_by_genotype_mgclusters_orderedLegend.png"), combined_umap, width = 12, height = 5.5, dpi = 300)
ggsave(file.path(outdir, "UMAP_microglia_split_by_genotype_mgclusters_orderedLegend.pdf"), combined_umap, width = 12, height = 5.5, device = cairo_pdf)
log_msg("Saved split UMAPs.")

# abundance summary
md <- as.data.frame(mg@meta.data) %>% filter(genotype %in% c("AA","GG"))
tab <- md %>% count(mg_cluster, genotype) %>% pivot_wider(names_from = genotype, values_from = n, values_fill = 0) %>%
  rename(count_AA = AA, count_GG = GG)
tab$count_AA[is.na(tab$count_AA)] <- 0
tab$count_GG[is.na(tab$count_GG)] <- 0
total_AA <- sum(tab$count_AA); total_GG <- sum(tab$count_GG)
tab <- tab %>% mutate(total = count_AA + count_GG,
                      pct_AA_of_cluster = round(100 * count_AA / pmax(total,1), 2),
                      pct_GG_of_cluster = round(100 * count_GG / pmax(total,1), 2),
                      pct_AA_of_allAA = round(100 * count_AA / pmax(total_AA,1), 2),
                      pct_GG_of_allGG = round(100 * count_GG / pmax(total_GG,1), 2),
                      diff_pct_cluster_AA_minus_GG = round(pct_AA_of_cluster - pct_GG_of_cluster, 2)) %>%
  arrange(as.numeric(as.character(mg_cluster)))
write.csv(tab, file.path(outdir, "microglia_cluster_genotype_abundance_summary.csv"), row.names = FALSE)
log_msg("Wrote abundance summary CSV.")

########################
# 10. STATISTICS: Chi-square + Fisher per-cluster (AA vs GG)
########################
ct <- table(mg$mg_cluster, mg$genotype)
chisq_res <- tryCatch(chisq.test(ct), error = function(e) e)
saveRDS(chisq_res, file.path(outdir, "chi_square_test_result.rds"))
log_msg("Chi-square result saved.")

# per-cluster Fisher tests and FDR correction
total_AA <- sum(mg$genotype == "AA")
total_GG <- sum(mg$genotype == "GG")
cluster_stats <- as.data.frame(mg@meta.data) %>%
  filter(genotype %in% c("AA","GG")) %>%
  count(mg_cluster, genotype) %>%
  pivot_wider(names_from = genotype, values_from = n, values_fill = 0) %>%
  rowwise() %>%
  mutate(
    fisher_obj = list(fisher.test(matrix(c(AA, GG, total_AA - AA, total_GG - GG), nrow = 2))),
    p_value = fisher_obj$p.value,
    odds_ratio = fisher_obj$estimate
  ) %>%
  ungroup() %>%
  mutate(padj = p.adjust(p_value, method = "fdr"),
         enrichment = case_when(odds_ratio > 1 ~ "AA-enriched", odds_ratio < 1 ~ "GG-enriched", TRUE ~ "Neutral")) %>%
  arrange(padj)
write.csv(cluster_stats, file.path(outdir, "microglia_cluster_fisher_tests_per_cluster.csv"), row.names = FALSE)
log_msg("Wrote Fisher test summary.")

########################
# 11. HEATMAPS (ComplexHeatmap) — user marker panels A & B (generic)
########################
# User should edit cluster ordering and marker panels below as needed.
# Example cluster_order_g (g* naming) — edit for your data:
cluster_order_g <- c("g4","g1","g9","g17","g11","g16","g6","g18","g7")
cluster_order_MG <- sub("^g", "MG", cluster_order_g)

# Marker panels (example); edit or replace with your preferred genes
homeostatic_A <- c("TMEM119","SALL1","CX3CR1","CSF1R","GPR34")
activated_A   <- c("HLA-DRA","CD68","CD74","CD40","ITGAM")
DAM_A         <- c("TREM2","APOE","SPP1","GPNMB","TYROBP")
genesA <- c(homeostatic_A, activated_A, DAM_A)
groupA <- rep(c("Homeostatic","Activated","DAM"),
              times = c(length(homeostatic_A), length(activated_A), length(DAM_A)))

IFN_B        <- c("IFIT1","IFIT3","ISG15")
Phago_B      <- c("MERTK","AXL","MFGE8")
Lysosomal_B  <- c("LAMP1","CTSD","TFEB")
MGnD_B       <- c("CH25H","LPL","CD9")
LDAM_B       <- c("PLIN2","CIDEC","FABP5")
Proteo_B     <- c("SQSTM1","HSPB1","BAG3")
Antigen_B    <- c("HLA-DRB1","CIITA","CD86")
genesB <- c(IFN_B, Phago_B, Lysosomal_B, MGnD_B, LDAM_B, Proteo_B, Antigen_B)
groupB <- rep(c("IFN","Phagocytic","Lysosomal","MGnD","LDAM","Proteostasis","Antigen"),
              times = c(length(IFN_B), length(Phago_B), length(Lysosomal_B), length(MGnD_B), length(LDAM_B), length(Proteo_B), length(Antigen_B)))

# Build matA_z and matB_z from avg_mat (if present)
if (!exists("avg_mat")) stop("avg_mat is required for heatmap generation (AverageExpression earlier).")
z_all <- per_gene_z(avg_mat)  # rows genes x cols clusters (clusters = column names)
# try g* or MG* naming matching
present_cols <- colnames(z_all)
if (all(cluster_order_g %in% present_cols)) use_cols <- cluster_order_g
else if (all(cluster_order_MG %in% present_cols)) use_cols <- cluster_order_MG
else {
  # try intersection
  use_cols <- intersect(cluster_order_g, present_cols)
  if (length(use_cols) == 0) use_cols <- intersect(cluster_order_MG, present_cols)
  if (length(use_cols) == 0) use_cols <- present_cols
  warning("Using intersection or default available cluster columns for heatmap: ", paste(head(use_cols,10), collapse = ", "))
}

matA_z <- z_all[intersect(rownames(z_all), genesA), use_cols, drop = FALSE]
matB_z <- z_all[intersect(rownames(z_all), genesB), use_cols, drop = FALSE]
colnames(matA_z) <- sub("^g", "MG", colnames(matA_z))
colnames(matB_z) <- sub("^g", "MG", colnames(matB_z))

# reorder rows to user-specified gene order (keep present genes only)
genesA_present <- genesA[genesA %in% rownames(matA_z)]
matA_z_ord <- matA_z[genesA_present, , drop = FALSE]
groupA_present <- groupA[genesA %in% genesA_present]

genesB_present <- genesB[genesB %in% rownames(matB_z)]
matB_z_ord <- matB_z[genesB_present, , drop = FALSE]
groupB_present <- groupB[genesB %in% genesB_present]
substate_levels <- c("IFN","Phagocytic","Lysosomal","Proteostasis","MGnD","LDAM","Antigen")
groupB_present_factor <- factor(groupB_present, levels = substate_levels)

# color map (blue-white-red)
col_fun <- colorRamp2(c(-2, 0, 2), c("#2166ac", "white", "#b2182b"))

# Heatmap A (broad categorization- homeostatic, activated, or DAM)
ht_A <- Heatmap(matA_z_ord,
                name = "z-score",
                col = col_fun,
                cluster_rows = FALSE,
                cluster_columns = FALSE,
                show_column_names = TRUE,
                show_row_names = TRUE,
                column_names_rot = 45,
                row_split = factor(groupA_present, levels = c("Homeostatic","Activated","DAM")),
                row_title_side = "left",
                row_title_gp = gpar(fontsize = 12, fontface = "bold"),
                gap = unit(6, "mm"),
                heatmap_legend_param = list(title = "z-score"))
png(file.path(outdir, "HeatmapA_major_signatures_MGnames.png"), width = 3000, height = 2400, res = 300)
draw(ht_A, heatmap_legend_side = "right", padding = unit(c(8,8,8,8), "mm"))
dev.off()
pdf(file.path(outdir, "HeatmapA_major_signatures_MGnames.pdf"), width = 11, height = 9)
draw(ht_A, heatmap_legend_side = "right", padding = unit(c(8,8,8,8), "mm"))
dev.off()
log_msg("Saved Heatmap A.")

# Heatmap B (more granular microglia substate signatures)
ord_idx_B <- order(as.integer(groupB_present_factor))
matB_z_ord2 <- matB_z_ord[ord_idx_B, , drop = FALSE]
groupB_ordered <- groupB_present_factor[ord_idx_B]
ht_B <- Heatmap(matB_z_ord2,
                name = "z-score",
                col = col_fun,
                cluster_rows = FALSE,
                cluster_columns = FALSE,
                show_column_names = TRUE,
                show_row_names = TRUE,
                column_names_rot = 45,
                row_split = groupB_ordered,
                row_title_side = "left",
                row_title_gp = gpar(fontsize = 12, fontface = "bold"),
                gap = unit(4, "mm"),
                heatmap_legend_param = list(title = "z-score"))
png(file.path(outdir, "HeatmapB_substate_ordered_labels_MGnames.png"), width = 3600, height = 3000, res = 300)
draw(ht_B, heatmap_legend_side = "right", padding = unit(c(10,8,8,8), "mm"))
dev.off()
pdf(file.path(outdir, "HeatmapB_substate_ordered_labels_MGnames.pdf"), width = 12, height = 10)
draw(ht_B, heatmap_legend_side = "right", padding = unit(c(10,8,8,8), "mm"))
dev.off()
log_msg("Saved Heatmap B.")

########################
# 12. PERSIST heatmap matrices & avg_mat into mg@misc and save canonical RDS
########################
stopifnot(exists("mg"), exists("avg_mat"), exists("matA_z"), exists("matB_z"))
if (is.null(mg@misc)) mg@misc <- list()
mg@misc$gene_signatures <- list(
  avg_mat = avg_mat,
  matA_z  = matA_z,
  matB_z  = matB_z,
  genesA  = genesA,
  genesB  = genesB,
  description = "Microglia gene signatures, module scoring, and heatmap matrices",
  created_on  = Sys.time()
)
str(mg@misc$gene_signatures, max.level = 1)

out_rds <- file.path(outdir, microglia_rds_name)
saveRDS(mg, file = out_rds, compress = save_rds_compress)
log_msg("Saved final microglia object to: ", out_rds)

########################
# 13. MANIFEST & DONE
########################
outputs <- list.files(outdir, full.names = TRUE)
write.csv(data.frame(outputs = outputs), file.path(outdir, "microglia_outputs_manifest.csv"), row.names = FALSE)
log_msg("Pipeline complete. Manifest written. Check: ", outdir)


