###############################################################################
# microglia_subset_signatures.R
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
log_msg("Saved microglia object to: ", out_rds)


########################
#added on to prepare for psuedobulking
########################

                      
#!/usr/bin/env Rscript
# microglia_assign_MG_STATE.R
# Generic script to assign cluster -> MG_STATE (Homeostatic / Activated / DAM)
# using per-gene z (cluster-level) and output canonical CSVs + RDS.
#
# Usage: source() or run in an active R session after setting params at top.

# --------- PARAMETERS (edit) -----------
input_rds <- "Microglia_gene_signatures.rds"   # input Seurat RDS (generic path)
outdir <- "outputs_microglia_MG_state"         # output folder (created if missing)
seurat_obj_varname <- "obj"                    # variable name to assign the Seurat object to
output_rds <- file.path(outdir, "Microglia_MG_state.rds")
# ----------------------------------------

# --------- Libraries ----------
suppressPackageStartupMessages({
  library(Seurat)
  library(dplyr)
  library(tidyr)
  library(readr)
  library(tibble)
})

# --------- small logger ----------
log_msg <- function(...) message(paste0(...))

dir.create(outdir, recursive = TRUE, showWarnings = FALSE)

# --------- 0. Load object (generic) ----------
stopifnot(file.exists(input_rds))
log_msg("Loading input RDS: ", input_rds)
tmp_obj <- readRDS(input_rds)
assign(seurat_obj_varname, tmp_obj)
rm(tmp_obj)
obj <- get(seurat_obj_varname)
stopifnot(inherits(obj, "Seurat"))
log_msg("Object loaded. DefaultAssay: ", DefaultAssay(obj))

# ensure basic metadata fields exist
required_meta <- c("mg_cluster", "sample_id", "genotype", "batch")
missing_meta <- setdiff(required_meta, colnames(obj@meta.data))
if (length(missing_meta) > 0) {
  stop("Missing required metadata columns: ", paste(missing_meta, collapse = ", "))
}

# --------- 1. Marker sets (edit if you want different panels) ----------
MG_HOMEOSTATIC <- c("P2RY12","TMEM119","SALL1","CX3CR1","CSF1R","FCRLS","GPR34","P2RX7")
MG_ACTIVATED   <- c("HLA-DRA","CD68","CD74","ITGAM","ITGB2","CD40","CD86","CD14")
MG_DAM         <- c("TREM2","APOE","SPP1","GPNMB","TYROBP","LPL","CD9","CST7")

state_marker_list <- list(
  Homeostatic = MG_HOMEOSTATIC,
  Activated   = MG_ACTIVATED,
  DAM         = MG_DAM
)

# Check presence of marker genes in RNA assay
DefaultAssay(obj) <- "RNA"
present_genes <- intersect(unique(unlist(state_marker_list)), rownames(obj))
missing_genes <- setdiff(unique(unlist(state_marker_list)), present_genes)
if (length(missing_genes) > 0) {
  log_msg("Warning - missing genes (will be skipped): ", paste(missing_genes, collapse = ", "))
} else {
  log_msg("All marker genes present.")
}

# --------- 2. Compute cluster-level mean expression and per-gene z across clusters ----------
# Use layer = "data" for Seurat v5 compatibility (log-normalized values in RNA@data)
expr_layer <- "data"
log_msg("Reading expression data (layer = '", expr_layer, "') and computing AverageExpression by mg_cluster.")
avg_list <- AverageExpression(obj, assays = "RNA", slot = expr_layer, group.by = "mg_cluster", verbose = FALSE)
avg_mat <- avg_list$RNA   # genes x clusters

cluster_names <- colnames(avg_mat)
log_msg("Clusters found in average expression matrix: ", paste(cluster_names, collapse = ", "))

marker_genes <- intersect(unique(unlist(state_marker_list)), rownames(avg_mat))
if (length(marker_genes) == 0) stop("None of the marker genes are present in avg_mat. Aborting.")

avg_markers_mat <- avg_mat[marker_genes, , drop = FALSE]

per_gene_z <- function(mat) {
  # per-gene z-score across clusters (rows = genes)
  z <- t(scale(t(as.matrix(mat)), center = TRUE, scale = TRUE))
  z[is.na(z)] <- 0
  return(z)
}
mat_z <- per_gene_z(avg_markers_mat)

# Write cluster gene mean expression and z
cluster_gene_mean_expr <- as.data.frame(avg_markers_mat) %>% rownames_to_column("gene")
fn_meanexpr <- file.path(outdir, "cluster_gene_mean_expr_genes_x_clusters.csv")
write_csv(cluster_gene_mean_expr, fn_meanexpr)

cluster_gene_z <- as.data.frame(mat_z) %>% rownames_to_column("gene")
fn_genez <- file.path(outdir, "cluster_gene_z_genes_x_clusters.csv")
write_csv(cluster_gene_z, fn_genez)

log_msg("Wrote: ", fn_meanexpr)
log_msg("Wrote: ", fn_genez)

# --------- 3. Compute cluster × state scores (mean z across marker genes for each state) ----------
clusters <- colnames(mat_z)
cluster_state_scores <- matrix(NA_real_, nrow = length(clusters), ncol = length(state_marker_list),
                               dimnames = list(clusters, names(state_marker_list)))

for (st in names(state_marker_list)) {
  genes_st <- intersect(state_marker_list[[st]], rownames(mat_z))
  if (length(genes_st) == 0) {
    cluster_state_scores[, st] <- 0
  } else {
    cluster_state_scores[, st] <- colMeans(mat_z[genes_st, , drop = FALSE], na.rm = TRUE)
  }
}
cluster_state_scores_df <- as.data.frame(cluster_state_scores) %>% rownames_to_column("cluster")
fn_cluster_state_scores <- file.path(outdir, "cluster_state_scores_meanZ.csv")
write_csv(cluster_state_scores_df, fn_cluster_state_scores)
log_msg("Wrote: ", fn_cluster_state_scores)

# --------- 4. Assign cluster -> MG_STATE by top mean z; compute delta and confidence ----------
# ensure cluster column present
cluster_state_scores_df <- cluster_state_scores_df %>% mutate(cluster = as.character(cluster))

# identify numeric state columns
state_score_cols <- cluster_state_scores_df %>% select(-cluster) %>% select(where(is.numeric)) %>% colnames()
if (length(state_score_cols) < 1) stop("No numeric state-score columns found in cluster_state_scores_df.")

# compute assigned state, deltas, cluster sizes
cluster_to_state_df <- cluster_state_scores_df %>%
  rowwise() %>%
  mutate(
    .vals = list(c_across(all_of(state_score_cols))),
    best_i = which.max(.vals),
    ord = list(order(.vals, decreasing = TRUE)),
    second_i = ifelse(length(.vals) >= 2, ord[[1]][2], NA_integer_),
    assigned_state = state_score_cols[best_i],
    top_mean_z = as.numeric(.vals[best_i]),
    second_mean_z = ifelse(!is.na(second_i), as.numeric(.vals[second_i]), NA_real_),
    delta_z = top_mean_z - second_mean_z,
    n_cells_in_cluster = as.integer(ifelse(as.character(cluster) %in% names(table(obj$mg_cluster)),
                                           table(obj$mg_cluster)[as.character(cluster)],
                                           ifelse(sub("^(?:g|G|MG|mg)?", "", cluster) %in% names(table(sub("^(?:g|G|MG|mg)?", "", obj$mg_cluster))),
                                                  table(sub("^(?:g|G|MG|mg)?", "", obj$mg_cluster))[sub("^(?:g|G|MG|mg)?", "", cluster)],
                                                  0))),
    note_flag = NA_character_
  ) %>%
  ungroup() %>%
  select(-.vals, -best_i, -ord, -second_i)

# set confidence levels (relaxed thresholds — adjust if you want stricter)
cluster_to_state_df <- cluster_to_state_df %>%
  mutate(
    conf_level = case_when(
      is.na(delta_z) ~ "low",
      delta_z >= 0.30 ~ "high",
      delta_z >= 0.10 ~ "moderate",
      TRUE ~ "low"
    ),
    note_flag = ifelse(n_cells_in_cluster < 20, paste0("small_n=", n_cells_in_cluster), note_flag)
  )

# Write cluster->state mapping
fn_cluster_to_state <- file.path(outdir, "cluster_to_state.csv")
write_csv(cluster_to_state_df, fn_cluster_to_state)
log_msg("Wrote: ", fn_cluster_to_state)

# --------- 5. Robust mapping keys and per-cell assignment ----------
cts <- cluster_to_state_df %>%
  mutate(mg_cluster_key = as.character(cluster),
         mg_cluster_key_num = sub("^(?:g|G|MG|mg)?", "", mg_cluster_key))

map_conf_primary <- setNames(cts$conf_level, cts$mg_cluster_key)
map_delta_primary <- setNames(cts$delta_z, cts$mg_cluster_key)
map_state_primary <- setNames(cts$assigned_state, cts$mg_cluster_key)

map_conf_num <- setNames(cts$conf_level, cts$mg_cluster_key_num)
map_delta_num <- setNames(cts$delta_z, cts$mg_cluster_key_num)
map_state_num <- setNames(cts$assigned_state, cts$mg_cluster_key_num)

cell_clusters <- as.character(obj$mg_cluster)
cell_clusters_num <- sub("^(?:g|G|MG|mg)?", "", cell_clusters)

conf_vec <- map_conf_primary[cell_clusters]
delta_vec <- map_delta_primary[cell_clusters]
state_vec <- map_state_primary[cell_clusters]

missing_idx <- which(is.na(conf_vec) | is.na(delta_vec) | is.na(state_vec))
if (length(missing_idx) > 0) {
  conf_vec[missing_idx] <- map_conf_num[cell_clusters_num[missing_idx]]
  delta_vec[missing_idx] <- map_delta_num[cell_clusters_num[missing_idx]]
  state_vec[missing_idx] <- map_state_num[cell_clusters_num[missing_idx]]
}

conf_vec[is.na(conf_vec)] <- NA_character_
delta_vec[is.na(delta_vec)] <- NA_real_
state_vec[is.na(state_vec)] <- NA_character_

stopifnot(length(conf_vec) == nrow(obj@meta.data))
stopifnot(length(delta_vec) == nrow(obj@meta.data))
stopifnot(length(state_vec) == nrow(obj@meta.data))

# assign into metadata (use @meta.data directly to avoid Seurat $ assignment edge cases)
obj@meta.data$MG_STATE_assigned   <- as.character(state_vec)
obj@meta.data$MG_STATE_confidence <- as.character(conf_vec)
obj@meta.data$MG_STATE_delta_z    <- as.numeric(delta_vec)
obj@meta.data$MG_STATE_note       <- ifelse(is.na(obj@meta.data$MG_STATE_assigned), "unassigned", NA_character_)

# sync $ accessor
obj$MG_STATE_assigned <- obj@meta.data$MG_STATE_assigned
obj$MG_STATE_confidence <- obj@meta.data$MG_STATE_confidence
obj$MG_STATE_delta_z <- obj@meta.data$MG_STATE_delta_z
obj$MG_STATE_note <- obj@meta.data$MG_STATE_note

log_msg("Per-cell MG_STATE fields added to object metadata.")

# --------- 6. Build sample-level manifest & sample x cluster/state matrices ----------
meta_df <- obj@meta.data %>% as.data.frame() %>% tibble::rownames_to_column("cell_barcode")

sample_summary <- meta_df %>%
  group_by(sample_id) %>%
  summarise(
    n_microglia_total = n(),
    pct_AA_in_sample = 100 * mean(genotype == "AA"),
    pct_GG_in_sample = 100 * mean(genotype == "GG"),
    .groups = "drop"
  )

sample_x_cluster_counts <- meta_df %>%
  count(sample_id, mg_cluster) %>%
  pivot_wider(names_from = mg_cluster, values_from = n, values_fill = 0) %>%
  arrange(sample_id)

sample_x_state_counts <- meta_df %>%
  mutate(MG_STATE_assigned = ifelse(is.na(MG_STATE_assigned), "Unassigned", MG_STATE_assigned)) %>%
  count(sample_id, MG_STATE_assigned) %>%
  pivot_wider(names_from = MG_STATE_assigned, values_from = n, values_fill = 0) %>%
  arrange(sample_id)

sample_meta_one <- meta_df %>%
  group_by(sample_id) %>%
  summarise(
    batch = first(batch),
    genotype = first(genotype),
    .groups = "drop"
  )

sample_manifest <- sample_meta_one %>%
  left_join(sample_summary, by = "sample_id") %>%
  left_join(sample_x_cluster_counts, by = "sample_id") %>%
  left_join(sample_x_state_counts, by = "sample_id")

fn_sample_manifest <- file.path(outdir, "sample_manifest.csv")
fn_sample_x_cluster <- file.path(outdir, "sample_x_cluster_counts.csv")
fn_sample_x_state <- file.path(outdir, "sample_x_state_counts.csv")

write_csv(sample_manifest, fn_sample_manifest)
write_csv(sample_x_cluster_counts, fn_sample_x_cluster)
write_csv(sample_x_state_counts, fn_sample_x_state)

log_msg("Wrote sample manifest and count matrices:\n - ", fn_sample_manifest, "\n - ", fn_sample_x_cluster, "\n - ", fn_sample_x_state)

# --------- 7. Eligibility decisions for sample×unit ----------
min_eligible <- 10
min_good <- 20

wide_to_long <- function(wide_df, unit_name = "unit") {
  wide_df %>%
    pivot_longer(-sample_id, names_to = unit_name, values_to = "n_cells")
}

sample_cluster_long <- wide_to_long(sample_x_cluster_counts, unit_name = "mg_cluster") %>%
  mutate(unit_type = "cluster",
         eligibility = case_when(
           n_cells >= min_good ~ "good",
           n_cells >= min_eligible ~ "eligible",
           TRUE ~ "too_few"
         ))

sample_state_long <- wide_to_long(sample_x_state_counts, unit_name = "MG_STATE_assigned") %>%
  mutate(unit_type = "state",
         eligibility = case_when(
           n_cells >= min_good ~ "good",
           n_cells >= min_eligible ~ "eligible",
           TRUE ~ "too_few"
         ))

sample_unit_eligibility <- bind_rows(sample_cluster_long, sample_state_long) %>%
  arrange(unit_type, sample_id, desc(n_cells))

fn_sample_unit_elig <- file.path(outdir, "sample_unit_eligibility.csv")
write_csv(sample_unit_eligibility, fn_sample_unit_elig)
log_msg("Wrote sample_unit_eligibility: ", fn_sample_unit_elig)

# --------- 8. Identify low-representation units & sample-driven clusters ----------
min_eligible <- 10

cluster_sample_geno_counts <- meta_df %>%
  count(mg_cluster, sample_id, genotype, name = "n_cells") %>%
  mutate(is_eligible = n_cells >= min_eligible) %>%
  group_by(mg_cluster, genotype) %>%
  summarise(n_samples_eligible = sum(is_eligible), .groups = "drop")

clusters_few_samples <- cluster_sample_geno_counts %>%
  filter(n_samples_eligible < 3) %>%
  arrange(mg_cluster, genotype)

state_sample_geno_counts <- meta_df %>%
  mutate(MG_STATE_assigned = ifelse(is.na(MG_STATE_assigned), "Unassigned", MG_STATE_assigned)) %>%
  count(MG_STATE_assigned, sample_id, genotype, name = "n_cells") %>%
  mutate(is_eligible = n_cells >= min_eligible) %>%
  group_by(MG_STATE_assigned, genotype) %>%
  summarise(n_samples_eligible = sum(is_eligible), .groups = "drop")

states_few_samples <- state_sample_geno_counts %>%
  filter(n_samples_eligible < 3) %>%
  arrange(MG_STATE_assigned, genotype)

cluster_sample_frac <- meta_df %>%
  count(mg_cluster, sample_id, name = "n_cells") %>%
  group_by(mg_cluster) %>%
  mutate(total_cells = sum(n_cells), frac = n_cells / total_cells) %>%
  ungroup()

sample_driven_clusters <- cluster_sample_frac %>%
  filter(frac > 0.80) %>%
  arrange(desc(frac))

fn_clusters_few_samples <- file.path(outdir, "clusters_few_samples_per_genotype.csv")
fn_states_few_samples   <- file.path(outdir, "states_few_samples_per_genotype.csv")
fn_sample_driven        <- file.path(outdir, "clusters_sample_driven.csv")

write_csv(clusters_few_samples, fn_clusters_few_samples)
write_csv(states_few_samples, fn_states_few_samples)
write_csv(sample_driven_clusters, fn_sample_driven)

log_msg("Wrote diagnostics:\n - ", fn_clusters_few_samples, "\n - ", fn_states_few_samples, "\n - ", fn_sample_driven)

# --------- 9. Export per_cell_metadata.csv (canonical) ----------
meta_out <- obj@meta.data %>% as.data.frame() %>% tibble::rownames_to_column(var = "cell_barcode")
scoreV2_cols <- grep("^scoreV2_", colnames(meta_out), value = TRUE)

percell_cols <- c(
  "cell_barcode",
  "sample_id",
  "mg_cluster",
  "MG_STATE_assigned",
  "MG_STATE_confidence",
  "MG_STATE_delta_z",
  "genotype",
  "batch",
  "nCount_RNA",
  "percent.mt",
  scoreV2_cols
)

existing_cols <- intersect(percell_cols, colnames(meta_out))
missing_cols  <- setdiff(percell_cols, colnames(meta_out))
if (length(missing_cols) > 0) {
  log_msg("Note: omitted missing per-cell columns: ", paste(missing_cols, collapse = ", "))
}

per_cell_metadata <- meta_out %>% select(all_of(existing_cols))
fn_percell <- file.path(outdir, "per_cell_metadata.csv")
write_csv(per_cell_metadata, fn_percell)
log_msg("Wrote per_cell_metadata: ", fn_percell)

# --------- 10. Save final Seurat object ----------
fn_rds_final <- output_rds
saveRDS(obj, fn_rds_final, compress = "xz")
log_msg("Saved final Seurat object to: ", fn_rds_final)

# --------- 11. Final summaries ----------
log_msg("\nCluster -> MG_STATE summary (first rows):")
print(head(cluster_to_state_df %>% select(cluster, assigned_state, top_mean_z, delta_z, conf_level, n_cells_in_cluster)))

log_msg("\nPer-cell MG_STATE distribution:")
print(table(obj@meta.data$MG_STATE_assigned, useNA = "ifany"))

log_msg("\nAll outputs written to: ", normalizePath(outdir))

#########################                      
#VISUALS FOR THE MG_STATES
##########################

# Generic plotting utilities for state-assigned single-cell datasets
# Drop into your GitHub script. Requires: Seurat, ggplot2, dplyr, tidyr, readr, patchwork, viridis, cowplot
# Example usage (after loading mg or meta_df):
# plot_umap_states(mg, outdir = outdir)
# plot_sample_state_tile(meta_df, outdir = outdir)
# plot_sample_composition_stacked(meta_df, outdir = outdir)

library(Seurat)
library(ggplot2)
library(dplyr)
library(tidyr)
library(readr)
library(patchwork)
library(viridis)
library(cowplot)

# --------------------------
# Helper: save PNG (high-res) + PDF (vector)
# --------------------------
save_both <- function(plot_obj, fname_base, outdir = ".", width = 6, height = 6, png_dpi = 600) {
  if (!dir.exists(outdir)) dir.create(outdir, recursive = TRUE)
  png_file <- file.path(outdir, paste0(fname_base, ".png"))
  pdf_file <- file.path(outdir, paste0(fname_base, ".pdf"))
  ggsave(png_file, plot = plot_obj, width = width, height = height, dpi = png_dpi, bg = "white")
  ggsave(pdf_file, plot = plot_obj, width = width, height = height, device = cairo_pdf)
  message("Wrote: ", png_file, " and ", pdf_file)
  invisible(list(png = png_file, pdf = pdf_file))
}

# --------------------------
# Default palettes (change as desired when calling)
# --------------------------
default_state_colors <- c(
  "Homeostatic" = "#1f78b4",
  "Activated"   = "#ff7f00",
  "DAM"         = "#6a3d9a",
  "Unassigned"  = "#bdbdbd"
)

default_geno_colors <- c(
  "AA" = "#D55E00",   # your preferred AA color
  "GG" = "#0072B2"    # your preferred GG color
)

# --------------------------
# Theme: journal-like, no gridlines
# --------------------------
theme_journal <- function(base_size = 12) {
  theme_classic(base_size = base_size) +
    theme(
      panel.grid = element_blank(),
      axis.title = element_blank(),
      axis.ticks = element_blank(),
      axis.text = element_blank(),
      legend.key = element_rect(fill = "transparent", colour = NA),
      legend.background = element_rect(fill = "transparent"),
      plot.title = element_text(size = base_size + 2, face = "bold")
    )
}

# --------------------------
# Function: plot_umap_states
# - input: Seurat object OR a data.frame with UMAP_1/UMAP_2 + state + genotype
# - main outputs: p_umap_all, p_umap_split (returns a list)
# --------------------------
plot_umap_states <- function(x,
                             state_col = "MG_STATE_assigned",
                             genotype_col = "genotype",
                             outdir = ".",
                             fname_base_all = "UMAP_states_all",
                             fname_base_split = "UMAP_states_split_by_genotype",
                             state_colors = default_state_colors,
                             genotype_colors = default_geno_colors,
                             point_size = 0.6,
                             point_alpha = 0.6,
                             png_dpi = 600,
                             width_all = 6, height_all = 6,
                             width_split = 12, height_split = 6) {

  # Build a metadata data.frame with UMAP coords
  if (inherits(x, "Seurat")) {
    if (!"umap" %in% names(x@reductions)) stop("Seurat object has no 'umap' reduction.")
    emb <- Embeddings(x, "umap")
    meta <- x@meta.data %>% tibble::rownames_to_column("cell_barcode")
    umap_df <- as.data.frame(emb) %>%
      tibble::rownames_to_column("cell_barcode") %>%
      left_join(meta, by = "cell_barcode")
    # ensure col names
    colnames(umap_df)[which(colnames(umap_df) == colnames(emb)[1])] <- "UMAP_1"
    colnames(umap_df)[which(colnames(umap_df) == colnames(emb)[2])] <- "UMAP_2"
  } else if (is.data.frame(x)) {
    umap_df <- x
    if (!all(c("UMAP_1","UMAP_2") %in% colnames(umap_df))) stop("data.frame must contain UMAP_1 and UMAP_2 columns.")
  } else {
    stop("x must be a Seurat object or a data.frame with UMAP_1/UMAP_2.")
  }

  # Prepare columns
  if (!state_col %in% colnames(umap_df)) {
    stop(paste0("State column '", state_col, "' not found in metadata."))
  }
  if (!genotype_col %in% colnames(umap_df)) {
    umap_df[[genotype_col]] <- NA_character_
  }

  umap_df <- umap_df %>%
    mutate(
      !!state_col := ifelse(is.na(.data[[state_col]]), "Unassigned", as.character(.data[[state_col]])),
      !!genotype_col := ifelse(is.na(.data[[genotype_col]]), "NA", as.character(.data[[genotype_col]]))
    )

  # factor ordering: preserve provided state_colors names order if possible
  present_states <- unique(umap_df[[state_col]])
  ordered_states <- intersect(names(state_colors), present_states)
  other_states <- setdiff(present_states, ordered_states)
  state_levels <- c(ordered_states, sort(other_states))
  umap_df[[state_col]] <- factor(umap_df[[state_col]], levels = state_levels)

  # Full UMAP (all cells)
  p_umap_all <- ggplot(umap_df, aes(x = UMAP_1, y = UMAP_2, color = .data[[state_col]])) +
    geom_point(size = point_size, alpha = point_alpha, stroke = 0) +
    scale_color_manual(values = state_colors, na.value = "#bdbdbd") +
    labs(color = "STATE") +
    theme_journal() +
    theme(legend.position = "right")

  # Split by genotype (only genotypes present will be faceted)
  umap_split_df <- umap_df %>% filter(.data[[genotype_col]] %in% names(genotype_colors))
  if (nrow(umap_split_df) == 0) {
    warning("No cells with genotypes matching genotype_colors; split-by-genotype plot will be empty.")
    p_umap_split <- NULL
  } else {
    p_umap_split <- ggplot(umap_split_df, aes(x = UMAP_1, y = UMAP_2, color = .data[[state_col]])) +
      geom_point(size = point_size, alpha = point_alpha, stroke = 0) +
      facet_wrap(as.formula(paste("~", genotype_col)), nrow = 1) +
      scale_color_manual(values = state_colors, na.value = "#bdbdbd") +
      theme_journal() +
      theme(strip.text = element_text(size = 12), legend.position = "right")
  }

  # Save to disk
  save_both(p_umap_all, fname_base_all, outdir = outdir, width = width_all, height = height_all, png_dpi = png_dpi)
  if (!is.null(p_umap_split)) {
    save_both(p_umap_split, fname_base_split, outdir = outdir, width = width_split, height = height_split, png_dpi = png_dpi)
  }

  invisible(list(umap_all = p_umap_all, umap_split = p_umap_split, umap_df = umap_df))
}

# --------------------------
# Function: plot_sample_state_tile
# - sample_id_col: column name for sample IDs
# - state_col: column name for state assignment
# - fills cells by log10(n+1); annotates raw counts (hide zeros)
# --------------------------
plot_sample_state_tile <- function(meta_df,
                                   sample_id_col = "sample_id",
                                   state_col = "MG_STATE_assigned",
                                   outdir = ".",
                                   fname_base = "sample_by_state_tile_log10_counts",
                                   state_colors = default_state_colors,
                                   palette_option = "magma",
                                   png_dpi = 600,
                                   width = 6.5, height = 6) {

  # validate
  if (!sample_id_col %in% colnames(meta_df)) stop("sample_id_col not found in meta_df")
  if (!state_col %in% colnames(meta_df)) stop("state_col not found in meta_df")

  df <- meta_df %>%
    mutate(!!state_col := ifelse(is.na(.data[[state_col]]), "Unassigned", as.character(.data[[state_col]]))) %>%
    count(.data[[sample_id_col]], .data[[state_col]], name = "n") %>%
    tidyr::complete(!!rlang::sym(sample_id_col), !!rlang::sym(state_col),
                    fill = list(n = 0))

  # compute log10(n+1)
  df <- df %>%
    mutate(log10_n1 = log10(n + 1))

  # order samples grouped by genotype if present (attempt to use genotype column if exists)
  if ("genotype" %in% colnames(meta_df)) {
    sample_order_df <- meta_df %>%
      group_by(.data[[sample_id_col]]) %>%
      summarise(total_cells = n(), genotype = first(genotype), .groups = "drop") %>%
      arrange(genotype, desc(total_cells))
    sample_order <- sample_order_df[[sample_id_col]]
  } else {
    sample_order <- df %>% group_by(.data[[sample_id_col]]) %>% summarise(total = sum(n)) %>% arrange(desc(total)) %>% pull(!!rlang::sym(sample_id_col))
  }

  df[[sample_id_col]] <- factor(df[[sample_id_col]], levels = sample_order)
  state_levels <- unique(df[[state_col]])
  df[[state_col]] <- factor(df[[state_col]], levels = state_levels)

  p_tile <- ggplot(df, aes_string(x = state_col, y = sample_id_col, fill = "log10_n1")) +
    geom_tile(color = NA) +
    geom_text(aes(label = ifelse(n > 0, n, "")), size = 3) +
    scale_fill_viridis_c(option = palette_option, na.value = "grey90", name = "log10(n+1)") +
    labs(x = "", y = "") +
    theme_journal(10) +
    theme(axis.text.y = element_text(size = 9), axis.text.x = element_text(size = 10), legend.position = "right")

  save_both(p_tile, fname_base, outdir = outdir, width = width, height = height, png_dpi = png_dpi)
  invisible(list(plot = p_tile, df = df))
}

# --------------------------
# Function: plot_sample_composition_stacked
# - stacked composition per sample (proportion), with genotype strip above (expects 'genotype' column in meta_df)
# --------------------------
plot_sample_composition_stacked <- function(meta_df,
                                           sample_id_col = "sample_id",
                                           state_col = "MG_STATE_assigned",
                                           genotype_col = "genotype",
                                           outdir = ".",
                                           fname_base = "sample_composition_by_state_stacked",
                                           state_colors = default_state_colors,
                                           genotype_colors = default_geno_colors,
                                           png_dpi = 600,
                                           width = 10, height = 5) {
  # validate
  if (!sample_id_col %in% colnames(meta_df)) stop("sample_id_col not found in meta_df")
  if (!state_col %in% colnames(meta_df)) stop("state_col not found in meta_df")
  if (!genotype_col %in% colnames(meta_df)) stop("genotype_col not found in meta_df")

  df <- meta_df %>%
    mutate(!!state_col := ifelse(is.na(.data[[state_col]]), "Unassigned", as.character(.data[[state_col]])),
           !!genotype_col := ifelse(is.na(.data[[genotype_col]]), "NA", as.character(.data[[genotype_col]]))) %>%
    group_by(.data[[sample_id_col]], .data[[genotype_col]], .data[[state_col]]) %>%
    summarise(n = n(), .groups = "drop") %>%
    group_by(.data[[sample_id_col]]) %>%
    mutate(prop = n / sum(n)) %>%
    ungroup()

  # Order samples by genotype blocks; within genotype by Homeostatic proportion if present
  sample_order_df <- df %>%
    group_by(.data[[sample_id_col]]) %>%
    summarise(total = sum(n), genotype = first(.data[[genotype_col]]), .groups = "drop")

  homeo_df <- df %>% filter(.data[[state_col]] == "Homeostatic") %>% select(.data[[sample_id_col]], prop) %>% rename(homeo_prop = prop)
  sample_order2 <- sample_order_df %>%
    left_join(homeo_df, by = sample_id_col) %>%
    arrange(genotype, desc(homeo_prop)) %>%
    pull(!!rlang::sym(sample_id_col))

  # ensure ordering includes all samples
  sample_order2 <- unique(c(sample_order2, sample_order_df[[sample_id_col]]))

  df[[sample_id_col]] <- factor(df[[sample_id_col]], levels = sample_order2)
  df[[state_col]] <- factor(df[[state_col]], levels = unique(df[[state_col]]))

  # Stacked bar (proportions)
  p_bar <- ggplot(df, aes_string(x = sample_id_col, y = "prop", fill = state_col)) +
    geom_bar(stat = "identity", position = "fill", width = 0.8, color = NA) +
    scale_fill_manual(values = state_colors) +
    scale_y_continuous(labels = scales::percent_format(accuracy = 1)) +
    labs(x = "", y = "Proportion", fill = "STATE") +
    theme_journal(10) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 9), legend.position = "right")

  # genotype strip (top)
  sample_geno_df <- sample_order_df %>%
    select(all_of(sample_id_col), genotype) %>%
    mutate(!!sample_id_col := factor(!!rlang::sym(sample_id_col), levels = sample_order2),
           genotype = factor(genotype, levels = names(genotype_colors)))

  # ensure genotype_colors only contains keys present in sample_geno_df
  present_genos <- intersect(names(genotype_colors), unique(as.character(sample_geno_df$genotype)))
  geno_palette <- genotype_colors[present_genos]

  p_geno_strip <- ggplot(sample_geno_df, aes_string(x = sample_id_col, y = 1, fill = "genotype")) +
    geom_tile() +
    scale_fill_manual(values = geno_palette, na.value = "#999999") +
    theme_void() +
    theme(legend.position = "none")

  combined <- p_geno_strip + p_bar + plot_layout(ncol = 1, heights = c(0.07, 0.93))

  save_both(combined, fname_base, outdir = outdir, width = width, height = height, png_dpi = png_dpi)
  invisible(list(plot = combined, df = df))
}

# --------------------------
# End of utilities
# --------------------------
# Example wrapper (uncomment to run in your script):
# meta_df <- mg@meta.data %>% tibble::rownames_to_column("cell_barcode")
# plot_umap_states(mg, outdir = outdir)
# plot_sample_state_tile(meta_df, outdir = outdir)
# plot_sample_composition_stacked(meta_df, outdir = outdir)
                                            
########################
# 13. MANIFEST & DONE
########################
outputs <- list.files(outdir, full.names = TRUE)
write.csv(data.frame(outputs = outputs), file.path(outdir, "microglia_outputs_manifest.csv"), row.names = FALSE)
log_msg("Pipeline complete. Manifest written. Check: ", outdir)


