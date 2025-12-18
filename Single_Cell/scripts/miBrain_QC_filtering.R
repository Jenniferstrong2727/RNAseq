#######################################################################
# miBrain_QC_filtering.R
#
# QC + integration + marker-based scoring for sc/snRNA data on Minerva:
#   - Per-sample QC (mito/ribo filters, MAD-based nFeature/nCount)
#   - Doublet detection using scDblFinder
#   - Restrict to protein-coding genes (from GTF)
#   - SCT-based integration across samples
#   - Harmony batch correction on PCA
#   - Clustering + UMAP
#   - Cluster markers (FindAllMarkers)
#   - NEW canonical marker sets (v2) + ambiguous list:
#       * filter markers to genes_of_interest
#       * module scores (scoreV2_*)
#       * top-module identity per cluster (moduleV2_celltype)
#       * heatmap + UMAP of module-based cell identities
#
# Usage (example on Minerva):
#   module load R/4.2.0
#   bsub -q premium -n 8 -W 12:00 -R "rusage[mem=64000]" \
#     -o logs/miBrain_QC_filtering.%J.out \
#     -e logs/miBrain_QC_filtering.%J.err \
#     Rscript Single_Cell/scripts/miBrain_QC_filtering.R
#
# Before running, EDIT the "USER INPUT" section below.
#######################################################################

suppressPackageStartupMessages({
  library(Seurat)
  library(dplyr)
  library(readxl)
  library(scater)
  library(scDblFinder)
  library(SingleCellExperiment)
  library(harmony)
  library(future)
  library(ggplot2)
  library(cowplot)
  library(rtracklayer)
  library(tidyr)
  library(tibble)
  library(pheatmap)
})

theme_set(theme_cowplot())

############################################################
## 0. USER INPUT: paths & parameters (EDIT THIS BLOCK)
############################################################

# Root project directory on Minerva (for MiBrain, or your project)
project_dir <- "/sc/arion/projects/ad-omics/Jennifer/scRNA_mibrain_P2RY12"

# Metadata: Excel with at least a 'sample_id' column.
# Optional: run_label, batch, line, genotype, clone, etc.
metadata_file <- file.path(project_dir, "metadata", "metadata_miBrain.xlsx")

# Directory with CellBender-filtered .h5 files
# Expect: <cb_base>/<sample_id>/<sample_id>_cellbender_filtered.h5
cb_base <- file.path(project_dir, "03_cellbender")

# GTF used to define protein-coding genes (gzipped)
gtf_path <- file.path(
  project_dir,
  "reference",
  "Homo_sapiens.GRCh38.109.gtf.gz"
)

# Output base directory for Seurat objects & results
output_base <- file.path(project_dir, "04_seurat_objects")

# Integration / clustering parameters
n_cores                 <- 8
n_pcs                   <- 50
harmony_dims            <- 1:30
umap_dims               <- 1:30
cluster_resolution      <- 0.6   # main resolution for seurat_clusters
integration_nfeatures   <- 3000  # SCT integration features

############################################################
## 1. Derived paths & setup (no need to edit)
############################################################

qc_raw_dir      <- file.path(output_base, "qc", "raw")
qc_singlet_dir  <- file.path(output_base, "qc", "singlets")
qc_plot_dir     <- file.path(output_base, "qc", "plots_per_sample")
integr_rds_dir  <- file.path(output_base, "integration", "rds")
integr_plot_dir <- file.path(output_base, "integration", "plots")
integr_mark_dir <- file.path(output_base, "integration", "markers")
integr_mod_dir  <- file.path(output_base, "integration", "module_scores")

message("Project dir:   ", project_dir)
message("Metadata file: ", metadata_file)
message("CellBender dir:", cb_base)
message("GTF path:      ", gtf_path)
message("Output base:   ", output_base, "\n")

dir.create(qc_raw_dir,      recursive = TRUE, showWarnings = FALSE)
dir.create(qc_singlet_dir,  recursive = TRUE, showWarnings = FALSE)
dir.create(qc_plot_dir,     recursive = TRUE, showWarnings = FALSE)
dir.create(integr_rds_dir,  recursive = TRUE, showWarnings = FALSE)
dir.create(integr_plot_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(integr_mark_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(integr_mod_dir,  recursive = TRUE, showWarnings = FALSE)

# Future plan for parallelization
plan("multicore", workers = n_cores)
options(future.globals.maxSize = 30 * 1024^3)  # 30 GB

set.seed(1234)

############################################################
## 2. Metadata + CellBender presence
############################################################

if (!file.exists(metadata_file)) {
  stop("Metadata file not found at: ", metadata_file)
}

meta <- readxl::read_xlsx(metadata_file)

if (!"sample_id" %in% colnames(meta)) {
  stop("Metadata file must contain a 'sample_id' column.")
}

meta <- meta %>%
  mutate(
    sample_id       = as.character(sample_id),
    cellbender_path = file.path(
      cb_base,
      sample_id,
      paste0(sample_id, "_cellbender_filtered.h5")
    ),
    has_cellbender  = file.exists(cellbender_path)
  )

message("CellBender status by sample:")
cols_to_show <- intersect(
  c("sample_id", "genotype", "batch", "has_cellbender"),
  colnames(meta)
)
print(meta %>% dplyr::select(dplyr::all_of(cols_to_show)))

meta_cb <- meta %>% dplyr::filter(has_cellbender)

if (nrow(meta_cb) == 0) {
  stop("No samples with CellBender outputs found. Check paths/metadata.")
}

message("\nRunning pipeline on ", nrow(meta_cb), " CellBender samples.\n")

############################################################
## 3. Protein-coding gene list from GTF
############################################################

if (!file.exists(gtf_path)) {
  stop("GTF not found at: ", gtf_path,
       "\nPlease copy it there or update gtf_path in the script.")
}

message("Importing GTF for protein-coding gene list...")
gtf <- rtracklayer::import(gtf_path)

protein_genes <- as.character(
  gtf[gtf$type == "gene" & gtf$gene_biotype == "protein_coding", ]$gene_name
)
protein_genes <- unique(protein_genes)
message("Number of protein-coding genes in GTF: ", length(protein_genes), "\n")

############################################################
## 4. Per-sample QC + doublets + protein-coding
############################################################

seurat_qc_list <- list()

for (i in seq_len(nrow(meta_cb))) {
  row <- meta_cb[i, ]
  s   <- row$sample_id

  message("\n==============================")
  message("Processing sample: ", s)
  message("CellBender file:   ", row$cellbender_path)

  if (!file.exists(row$cellbender_path)) {
    warning("CellBender file missing for sample ", s, " at ", row$cellbender_path,
            " — skipping.")
    next
  }

  counts <- Read10X_h5(row$cellbender_path)

  project_label <- if ("run_label" %in% colnames(meta_cb)) {
    as.character(row$run_label)
  } else {
    s
  }

  obj <- CreateSeuratObject(
    counts       = counts,
    project      = project_label,
    min.cells    = 3,
    min.features = 200
  )

  # Save raw object
  saveRDS(
    obj,
    file = file.path(qc_raw_dir, paste0(s, "_raw.rds"))
  )

  # Attach useful metadata if present
  obj$sample_id <- s
  if ("run_label" %in% colnames(meta_cb)) obj$run_label <- as.character(row$run_label)
  if ("batch"     %in% colnames(meta_cb)) obj$batch     <- as.character(row$batch)
  if ("line"      %in% colnames(meta_cb)) obj$line      <- as.character(row$line)
  if ("genotype"  %in% colnames(meta_cb)) obj$genotype  <- as.character(row$genotype)
  if ("clone"     %in% colnames(meta_cb)) obj$clone     <- as.character(row$clone)

  message("Initial cell count: ", ncol(obj))

  # QC metrics
  obj[["percent.mt"]]   <- PercentageFeatureSet(obj, pattern = "^MT-")
  obj[["percent.ribo"]] <- PercentageFeatureSet(obj, pattern = "^RPS|^RPL")

  # QC violin plot (raw)
  vln_raw <- VlnPlot(
    obj,
    features = c("nFeature_RNA", "nCount_RNA", "percent.mt", "percent.ribo"),
    ncol = 4, pt.size = 0
  ) + ggtitle(paste0(s, " - raw QC"))

  ggsave(
    file.path(qc_plot_dir, paste0(s, "_qc_raw.png")),
    vln_raw,
    width = 12, height = 4, dpi = 300
  )

  # Mito/ribo filters (edit thresholds if needed)
  obj <- subset(obj, subset = percent.mt < 10 & percent.ribo < 15)
  message("After mt<10 & ribo<15: ", ncol(obj), " cells")

  # MAD-based nFeature/nCount filtering
  nf_med <- median(obj$nFeature_RNA); nf_mad <- mad(obj$nFeature_RNA)
  nc_med <- median(obj$nCount_RNA);   nc_mad <- mad(obj$nCount_RNA)

  nf_upper <- nf_med + 5 * nf_mad
  nc_upper <- nc_med + 5 * nc_mad

  obj <- subset(
    obj,
    subset = nFeature_RNA >= 200 & nFeature_RNA <= nf_upper &
             nCount_RNA   >= 500 & nCount_RNA   <= nc_upper
  )
  message("After MAD filters: ", ncol(obj), " cells")

  # Doublets via scDblFinder
  sce <- as.SingleCellExperiment(obj)
  sce <- logNormCounts(sce)
  sce <- scDblFinder(sce)

  obj_singlet <- as.Seurat(sce, counts = "counts", data = "logcounts")
  obj_singlet <- subset(obj_singlet, subset = scDblFinder.class == "singlet")
  message("After scDblFinder singlets: ", ncol(obj_singlet), " cells")

  # Protein-coding only
  keep_genes <- intersect(rownames(obj_singlet), protein_genes)
  obj_singlet <- subset(obj_singlet, features = keep_genes)
  message("After protein-coding filter: ", nrow(obj_singlet), " genes")

  # Save per-sample QC object
  saveRDS(
    obj_singlet,
    file = file.path(qc_singlet_dir, paste0(s, "_qc_singlets.rds"))
  )

  seurat_qc_list[[s]] <- obj_singlet

  # QC violin plot (post-QC)
  vln_post <- VlnPlot(
    obj_singlet,
    features = c("nFeature_RNA", "nCount_RNA", "percent.mt", "percent.ribo"),
    ncol = 4, pt.size = 0
  ) + ggtitle(paste0(s, " - post-QC"))

  ggsave(
    file.path(qc_plot_dir, paste0(s, "_qc_post.png")),
    vln_post,
    width = 12, height = 4, dpi = 300
  )
}

saveRDS(
  seurat_qc_list,
  file = file.path(qc_singlet_dir, "all_samples_qc_singlets_list.rds")
)

############################################################
## 5. SCT integration + Harmony + UMAP + clustering
############################################################

message("\nStarting SCT integration with Harmony...\n")

seurat_sct_list <- lapply(seurat_qc_list, function(x) {
  SCTransform(x, vst.flavor = "v2", verbose = FALSE)
})

features <- SelectIntegrationFeatures(
  object.list = seurat_sct_list,
  nfeatures   = integration_nfeatures
)

seurat_sct_list <- PrepSCTIntegration(
  object.list     = seurat_sct_list,
  anchor.features = features
)

anchors <- FindIntegrationAnchors(
  object.list          = seurat_sct_list,
  normalization.method = "SCT",
  anchor.features      = features
)

integrated <- IntegrateData(
  anchorset           = anchors,
  normalization.method = "SCT"
)

saveRDS(
  integrated,
  file = file.path(integr_rds_dir, "integrated_SCT_raw.rds")
)

integrated <- RunPCA(integrated, npcs = n_pcs, verbose = FALSE)

integrated <- RunHarmony(
  object        = integrated,
  group.by.vars = "sample_id",
  reduction     = "pca",
  dims          = 1:n_pcs,
  assay.use     = "SCT",
  verbose       = TRUE
)

# Neighbors + clustering + UMAP
integrated <- FindNeighbors(integrated, reduction = "harmony", dims = harmony_dims)
integrated <- FindClusters(integrated, resolution = cluster_resolution)

# Use this clustering as seurat_clusters
integrated$seurat_clusters <- integrated$seurat_clusters
Idents(integrated) <- "seurat_clusters"

integrated <- RunUMAP(
  integrated,
  reduction = "harmony",
  dims      = umap_dims,
  min.dist  = 0.1
)

saveRDS(
  integrated,
  file = file.path(integr_rds_dir, "integrated_harmony_umap.rds")
)

# Basic UMAP plots
p_clusters <- DimPlot(
  integrated,
  reduction = "umap",
  group.by  = "seurat_clusters",
  label     = TRUE, repel = TRUE
) + ggtitle("Integrated – clusters (seurat_clusters)")

ggsave(
  file.path(integr_plot_dir, "UMAP_clusters.png"),
  p_clusters, width = 6, height = 5, dpi = 300
)

if ("genotype" %in% colnames(integrated@meta.data)) {
  p_genotype <- DimPlot(
    integrated,
    reduction = "umap",
    group.by  = "genotype"
  ) + ggtitle("Integrated – genotype")
  ggsave(
    file.path(integr_plot_dir, "UMAP_genotype.png"),
    p_genotype, width = 6, height = 5, dpi = 300
  )
}

############################################################
## 6. NEW canonical marker sets (v2) + ambiguous markers
############################################################

# Microglia (canonical)
markers_microglia <- c(
  "P2RY12", "TMEM119", "SALL1", "CX3CR1",
  "TREM2", "CD68", "PTPRC", "CD45", "HLA-DRA", "MHCII", "CD40", "CD206",
  "C1QA", "C1QB", "C1QC", "C1Q",
  "ITGAM", "CD11b", "CD9", "ITGAX",
  "GPNMB", "HEXB", "FCRLS", "GPR34",
  "MERTK", "PROS1", "TYRO3",
  "AIF1", "SPP1", "CSF1R", "CD74",
  "C3", "SPI1"
)

# Neurons (canonical)
markers_neuron <- c(
  "MAP2", "B3GAT1",
  "RBFOX3",
  "SYN1", "SYNAPSIN1", "SYP",
  "DCX", "NEUROD1", "ASCL1",
  "CUX1", "BCL11B", "CTIP2", "FOXP2", "TBR1",
  "TUBB3",
  "TH",
  "GFRA1", "SLC17A6",
  "CX3CL1",
  "STMN2", "NEFL", "NEFM"
)

# Astrocytes (canonical)
markers_astro <- c(
  "GFAP", "AQP4", "ALDH1L1",
  "EAAT2",
  "GLAST",
  "SLC1A2", "SLC1A3",
  "S100B", "SOX9", "GLUL",
  "SERPINA3", "CHI3L1",
  "GLT1",
  "GJA1",
  "SPARCL1", "FABP7", "CLU"
)

# Oligodendrocytes / OPCs (canonical)
markers_oligo <- c(
  "MBP", "PLP1", "SOX10", "MYRF",
  "MOG", "MAG", "MOBP",
  "O4",
  "PDGFRA", "NG2", "CSPG4",
  "OLIG1", "OLIG2",
  "CLDN11"
)

# Endothelial (canonical)
markers_endo <- c(
  "CD31", "PECAM1",
  "FLT1", "KDR",
  "TIE2", "TEK",
  "VWF",
  "VE-cadherin", "CDH5",
  "CLDN5",
  "ABCB1",
  "ESAM",
  "RBP7",
  "SLC2A1"
)

# Pericytes (canonical)
markers_pericyte <- c(
  "PDGFRB",
  "RGS5",
  "ACTA2",
  "TAGLN",
  "KCNJ8", "ABCC9",
  "COL4A1",
  "CD146", "MCAM"
)

# Ambiguous markers (used for filtering markers but NOT for module scores)
markers_ambiguous <- c(
  "VIM",
  "CD44",
  "ICAM1", "VCAM1",
  "ITGAV",
  "PAMR1",
  "CKB", "CPE",
  "SLC2A5",
  "SOX2",
  "CXXC5",
  "TMEM163",
  "CD105", "ENG",
  "Lef1", "Fzd3", "Notum",
  "Apcdd1", "Axin2", "dixdc1", "Tnfrsf19"
)

# Canonical marker sets used for module scores (v2)
marker_sets_v2 <- list(
  Microglia        = markers_microglia,
  Neurons          = markers_neuron,
  Astrocytes       = markers_astro,
  Oligodendrocytes = markers_oligo,
  Endothelial      = markers_endo,
  Pericytes        = markers_pericyte
)

# Build marker_table = gene × celltype × category
marker_table <- bind_rows(
  tibble(
    gene      = markers_microglia,
    celltype  = "Microglia",
    category  = "canonical"
  ),
  tibble(
    gene      = markers_neuron,
    celltype  = "Neurons",
    category  = "canonical"
  ),
  tibble(
    gene      = markers_astro,
    celltype  = "Astrocytes",
    category  = "canonical"
  ),
  tibble(
    gene      = markers_oligo,
    celltype  = "Oligodendrocytes",
    category  = "canonical"
  ),
  tibble(
    gene      = markers_endo,
    celltype  = "Endothelial",
    category  = "canonical"
  ),
  tibble(
    gene      = markers_pericyte,
    celltype  = "Pericytes",
    category  = "canonical"
  ),
  tibble(
    gene      = markers_ambiguous,
    celltype  = "Ambiguous",
    category  = "ambiguous"
  )
) %>%
  distinct(gene, .keep_all = TRUE)

genes_of_interest <- unique(marker_table$gene)

############################################################
## 7. Cluster markers (FindAllMarkers) + filtering
############################################################

DefaultAssay(integrated) <- DefaultAssay(integrated)  # no change, just explicit
Idents(integrated)        <- "seurat_clusters"

markers_all <- FindAllMarkers(
  integrated,
  only.pos        = TRUE,
  test.use        = "wilcox",
  logfc.threshold = 0.25,
  min.pct         = 0.1
)

write.csv(
  markers_all,
  file = file.path(integr_mark_dir, "markers_all_clusters.csv"),
  row.names = FALSE
)

# Subset to genes in marker_table (canonical + ambiguous)
markers_filtered <- markers_all %>%
  inner_join(marker_table, by = "gene") %>%
  arrange(as.numeric(as.character(cluster)), desc(avg_log2FC))

write.csv(
  markers_filtered,
  file = file.path(integr_mark_dir, "markers_clusters_filtered_to_interest.csv"),
  row.names = FALSE
)

# Absolute expression matrix (AverageExpression) for genes_of_interest
avg_list <- AverageExpression(
  integrated,
  group.by = "seurat_clusters",
  verbose  = FALSE
)
avg_mat <- avg_list[[1]]  # genes x clusters

avg_df_long <- avg_mat %>%
  as.data.frame() %>%
  rownames_to_column("gene") %>%
  pivot_longer(
    cols      = -gene,
    names_to  = "cluster",
    values_to = "avg_expression"
  )

absolute_expr_markers <- avg_df_long %>%
  inner_join(marker_table, by = "gene") %>%
  mutate(
    cluster_num = as.numeric(gsub("[^0-9-]", "", cluster))
  ) %>%
  arrange(cluster_num, desc(avg_expression)) %>%
  select(
    cluster = cluster_num,
    gene,
    celltype,
    category,
    avg_expression
  )

write.csv(
  absolute_expr_markers,
  file = file.path(
    integr_mark_dir,
    "markers_clusters_absolute_expression_filtered.csv"
  ),
  row.names = FALSE
)

############################################################
## 8. NEW module scores (scoreV2_) + top module per cluster
############################################################

DefaultAssay(integrated) <- DefaultAssay(integrated)
Idents(integrated)        <- "seurat_clusters"

# Add module scores with new prefix (v2)
integrated <- AddModuleScore(
  integrated,
  features = marker_sets_v2,
  name     = "celltypeScore_v2_"
)

meta <- integrated@meta.data

# Rename these new columns to scoreV2_<celltype>
score_cols_raw_v2 <- grep("^celltypeScore_v2_", colnames(meta), value = TRUE)
new_names_v2      <- paste0("scoreV2_", names(marker_sets_v2))

colnames(integrated@meta.data)[
  match(score_cols_raw_v2, colnames(integrated@meta.data))
] <- new_names_v2

meta <- integrated@meta.data
score_cols_v2 <- grep("^scoreV2_", colnames(meta), value = TRUE)

# Mean module scores per cluster
cluster_means_v2 <- meta %>%
  mutate(seurat_clusters = as.character(seurat_clusters)) %>%
  group_by(seurat_clusters) %>%
  summarise(
    across(all_of(score_cols_v2), mean, .names = "{.col}"),
    .groups = "drop"
  ) %>%
  arrange(as.numeric(seurat_clusters))

# Long format, pick top module per cluster
cluster_top_v2 <- cluster_means_v2 %>%
  pivot_longer(
    cols      = all_of(score_cols_v2),
    names_to  = "module",
    values_to = "mean_score"
  ) %>%
  group_by(seurat_clusters) %>%
  slice_max(mean_score, n = 1, with_ties = FALSE) %>%
  ungroup() %>%
  mutate(
    module_clean = sub("^scoreV2_", "", module)
  ) %>%
  arrange(as.numeric(seurat_clusters)) %>%
  select(
    seurat_clusters,
    top_module     = module_clean,
    top_mean_score = mean_score
  )

write.csv(
  cluster_top_v2,
  file = file.path(
    integr_mod_dir,
    "moduleV2_top_module_by_cluster.csv"
  ),
  row.names = FALSE
)

############################################################
## 9. moduleV2_celltype per cell + UMAP
############################################################

meta <- integrated@meta.data

cluster_to_module_v2 <- setNames(
  cluster_top_v2$top_module,
  cluster_top_v2$seurat_clusters
)

cell_clusters <- as.character(meta$seurat_clusters)
module_vec_v2 <- cluster_to_module_v2[cell_clusters]
module_vec_v2 <- unname(module_vec_v2)

integrated$moduleV2_celltype <- factor(
  module_vec_v2,
  levels = sort(unique(cluster_top_v2$top_module))
)

# Sanity table
print(
  table(
    Cluster        = integrated$seurat_clusters,
    ModuleV2_label = integrated$moduleV2_celltype
  )
)

p_umap_module_v2 <- DimPlot(
  integrated,
  reduction = "umap",
  group.by  = "moduleV2_celltype",
  label     = TRUE,
  repel     = TRUE
) + ggtitle("UMAP – Canonical module score identities (v2)")

ggsave(
  file.path(integr_plot_dir, "UMAP_moduleV2_celltype.png"),
  plot   = p_umap_module_v2,
  width  = 7,
  height = 6,
  dpi    = 600
)

ggsave(
  file.path(integr_plot_dir, "UMAP_moduleV2_celltype.pdf"),
  plot   = p_umap_module_v2,
  width  = 7,
  height = 6,
  device = cairo_pdf
)

############################################################
## 10. Heatmap of z-scored module scores (v2)
############################################################

meta <- integrated@meta.data
score_cols_v2 <- grep("^scoreV2_", colnames(meta), value = TRUE)
if (length(score_cols_v2) == 0) {
  stop("No columns starting with 'scoreV2_' found in meta.data.")
}

cluster_means_v2 <- meta %>%
  mutate(seurat_clusters = as.character(seurat_clusters)) %>%
  group_by(seurat_clusters) %>%
  summarise(
    across(all_of(score_cols_v2), mean, .names = "{.col}"),
    .groups = "drop"
  ) %>%
  arrange(as.numeric(seurat_clusters))

mat <- as.matrix(cluster_means_v2[, score_cols_v2, drop = FALSE])
rownames(mat) <- cluster_means_v2$seurat_clusters
storage.mode(mat) <- "numeric"

heatmap_mat_z <- scale(mat)
heatmap_mat_z <- as.matrix(heatmap_mat_z)
rownames(heatmap_mat_z) <- rownames(mat)
storage.mode(heatmap_mat_z) <- "numeric"

colnames(heatmap_mat_z) <- gsub("^scoreV2_", "", colnames(heatmap_mat_z))

heatmap_colors <- colorRampPalette(c("navy", "white", "firebrick3"))(100)

png(
  file.path(integr_mod_dir, "heatmap_module_scores_clusters_V2_noAnno.png"),
  width = 7, height = 6, units = "in", res = 600
)
pheatmap(
  heatmap_mat_z,
  cluster_rows  = TRUE,
  cluster_cols  = TRUE,
  color         = heatmap_colors,
  main          = "Z-normalized Module Scores by Cluster (V2 markers)",
  fontsize      = 12,
  fontsize_row  = 10,
  fontsize_col  = 12,
  border_color  = NA
)
dev.off()

pdf(
  file.path(integr_mod_dir, "heatmap_module_scores_clusters_V2_noAnno.pdf"),
  width = 7, height = 6
)
pheatmap(
  heatmap_mat_z,
  cluster_rows  = TRUE,
  cluster_cols  = TRUE,
  color         = heatmap_colors,
  main          = "Z-normalized Module Scores by Cluster (V2 markers)",
  fontsize      = 12,
  fontsize_row  = 10,
  fontsize_col  = 12,
  border_color  = NA
)
dev.off()

############################################################
## 11. Save final integrated object
############################################################

saveRDS(
  integrated,
  file = file.path(
    integr_rds_dir,
    "integrated_harmony_umap_resX_moduleV2.rds"
  )
)

message("\n==== miBrain_QC_filtering pipeline complete ====\n")
print(sessionInfo())


