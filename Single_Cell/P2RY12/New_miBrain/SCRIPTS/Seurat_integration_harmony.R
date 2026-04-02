###############################################
# 03_seurat_integration_harmony.R
#
# Generic Seurat + Harmony Integration Script
#
# Performs:
#   - Loading QC’d Seurat objects (from CellBender or preprocessing)
#   - SCTransform normalization
#   - Seurat SCT integration
#   - Harmony batch correction
#   - UMAP + clustering
#   - Marker detection (FindAllMarkers)
#   - Optional: module scoring
#
# Intended to run after Cell Ranger + CellBender.
#
# Author: <Your Name>
# GitHub: <Your GitHub Link>
###############################################

suppressPackageStartupMessages({
  library(Seurat)
  library(dplyr)
  library(harmony)
  library(ggplot2)
  library(cowplot)
})

theme_set(theme_cowplot())

################################################
## ---------- USER CONFIGURATION --------------
################################################

# Folder where QC Seurat objects (list) is stored
qc_list_path <- "<PATH_TO_QC_SEURAT_LIST>.rds"
# Example: "04_seurat_objects/all_samples_qc_list.rds"

# Output base directory
output_dir <- "<OUTPUT_DIRECTORY>"
# Example: "06_integration"

# Integration settings
nfeatures_integration <- 3000
npcs_pca              <- 50
harmony_dims          <- 1:50
umap_dims             <- 1:30
cluster_resolution    <- 0.3

################################################
## ---------- SETUP OUTPUT DIRS ---------------
################################################

dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
rds_dir      <- file.path(output_dir, "rds")
plots_dir    <- file.path(output_dir, "plots")
markers_dir  <- file.path(output_dir, "markers")
module_dir   <- file.path(output_dir, "module_scores")

for (d in c(rds_dir, plots_dir, markers_dir, module_dir)) {
  dir.create(d, recursive = TRUE, showWarnings = FALSE)
}

################################################
## ---------- 1. LOAD QC SEURAT LIST ----------
################################################

if (!file.exists(qc_list_path)) {
  stop("QC list not found: ", qc_list_path)
}

seurat_qc_list <- readRDS(qc_list_path)
message("Loaded ", length(seurat_qc_list), " Seurat objects.")

# Ensure sample_id exists
for (nm in names(seurat_qc_list)) {
  seurat_qc_list[[nm]]$sample_id <- nm
}

################################################
## ---------- 2. SCT INTEGRATION --------------
################################################

message("Running SCTransform on all samples...")

seurat_sct_list <- lapply(seurat_qc_list, function(x) {
  SCTransform(x, verbose = FALSE)
})

features <- SelectIntegrationFeatures(
  object.list = seurat_sct_list,
  nfeatures   = nfeatures_integration
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
  anchorset            = anchors,
  normalization.method = "SCT"
)

saveRDS(integrated, file.path(rds_dir, "integrated_SCT_raw.rds"))

################################################
## ---------- 3. PCA + HARMONY + UMAP ----------
################################################

integrated <- RunPCA(integrated, npcs = npcs_pca, verbose = FALSE)

integrated <- RunHarmony(
  object        = integrated,
  group.by.vars = "sample_id",
  reduction.use = "pca",
  dims.use      = harmony_dims,
  verbose       = TRUE
)

integrated <- FindNeighbors(integrated, reduction = "harmony", dims = umap_dims)
integrated <- FindClusters(integrated, resolution = cluster_resolution)

integrated <- RunUMAP(
  integrated,
  reduction = "harmony",
  dims      = umap_dims
)

saveRDS(integrated, file.path(rds_dir, "integrated_harmony_umap.rds"))

################################################
## ---------- 4. UMAP PLOTS ---------------------
################################################

p1 <- DimPlot(integrated, reduction = "umap", group.by = "sample_id")
p2 <- DimPlot(integrated, reduction = "umap", group.by = "seurat_clusters", label = TRUE)

ggsave(file.path(plots_dir, "UMAP_by_sample.png"), p1, width=6, height=5)
ggsave(file.path(plots_dir, "UMAP_by_cluster.png"), p2, width=6, height=5)

################################################
## ---------- 5. MARKER DETECTION --------------
################################################

Idents(integrated) <- "seurat_clusters"

markers <- FindAllMarkers(
  integrated,
  only.pos = TRUE,
  test.use = "wilcox",
  logfc.threshold = 0.25,
  min.pct = 0.1
)

write.csv(markers, file.path(markers_dir, "markers_all_clusters.csv"), row.names=FALSE)

################################################
## ---------- 6. OPTIONAL MODULE SCORES --------
################################################

module_sets <- list(
  Neurons   = c("RBFOX3","MAP2","SYN1","TUBB3"),
  Astrocytes = c("GFAP","ALDH1L1","SLC1A2"),
  Microglia  = c("P2RY12","TMEM119","CX3CR1")
)

module_sets <- lapply(module_sets, function(g) intersect(g, rownames(integrated)))
module_sets <- module_sets[lengths(module_sets) > 0]

if (length(module_sets) > 0) {
  integrated <- AddModuleScore(integrated, features = module_sets, name = names(module_sets))
  write.csv(integrated@meta.data, file.path(module_dir, "module_scores.csv"))
}

message("\n==== Seurat Integration + Harmony COMPLETE ====\n")
sessionInfo()

