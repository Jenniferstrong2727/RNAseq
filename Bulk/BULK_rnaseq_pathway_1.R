#!/usr/bin/env Rscript

# BULK_rnaseq_pathway_1.R
# Generic pathway analysis for bulk RNA-seq DE results
# - GSEA (fgseaMultilevel) using MSigDB via msigdbr (Symbols)
# - ORA (clusterProfiler) for GO BP, KEGG, and Hallmark (ENTREZ)
# Inputs: DESeq2 results RDS from a prior DEG step (res_raw + res_shr), optional dds
# Outputs: CSV, RDS, and PNGs under <outdir>/<project>/PATHWAYS/

suppressPackageStartupMessages({
  library(optparse)
  library(dplyr); library(tidyr); library(stringr); library(readr)
  library(ggplot2)
  library(AnnotationDbi); library(org.Hs.eg.db)
  library(msigdbr); library(fgsea)
  library(clusterProfiler); library(enrichplot)
  library(tibble)
})

# ---------- CLI ----------
opt_list <- list(
  make_option("--res_raw", type="character", help="Path to results_raw RDS (DESeq2::results on your contrast) [required]"),
  make_option("--res_shr", type="character", help="Path to results_shrunken RDS (lfcShrink) [required]"),
  make_option("--dds",     type="character", default=NA, help="(Optional) Path to dds RDS (for provenance)"),
  make_option("--outdir",  type="character", help="Output directory root [required]"),
  make_option("--project", type="character", default="Project", help="Project name tag for filenames"),
  make_option("--species", type="character", default="Homo sapiens", help="Species for msigdbr/orgDb"),
  make_option("--msig_collection", type="character", default="H", help="MSigDB collection code (e.g. H, C2, C5, etc)"),
  make_option("--padj_cutoff", type="double", default=0.05, help="padj threshold for ORA"),
  make_option("--lfc_cutoff",  type="double", default=0.5,  help="|LFC| threshold for ORA significance"),
  make_option("--make_cnet", action="store_true", default=FALSE, help="Try cnetplot (requires enrichplot/ggraph working)")
)
opt <- parse_args(OptionParser(option_list = opt_list))

stopifnot(!is.null(opt$res_raw), !is.null(opt$res_shr), !is.null(opt$outdir))
stopifnot(file.exists(opt$res_raw), file.exists(opt$res_shr))

# ---------- Paths ----------
base_dir  <- file.path(opt$outdir, opt$project)
path_dir  <- file.path(base_dir, "PATHWAYS")
plots_dir <- file.path(path_dir, "plots")
dir.create(path_dir,  recursive = TRUE, showWarnings = FALSE)
dir.create(plots_dir, recursive = TRUE, showWarnings = FALSE)

# ---------- Load inputs ----------
res_raw <- readRDS(opt$res_raw)  # DESeq2::results (raw)
res_shr <- readRDS(opt$res_shr)  # DESeq2::lfcShrink result

# To data.frames with gene_id
res_full <- as.data.frame(res_raw) %>% rownames_to_column("gene_id")
shr_df   <- as.data.frame(res_shr) %>% rownames_to_column("gene_id")

# ---------- ID helpers ----------
strip_ver <- function(x) sub("\\.\\d+$", "", x)
looks_ens <- function(x) grepl("^ENSG\\d+", x)

# Map ENSEMBL -> SYMBOL for a data.frame with gene_id
add_symbols <- function(df) {
  df$ENSEMBL_clean <- ifelse(looks_ens(df$gene_id), strip_ver(df$gene_id), NA_character_)
  df$gene_sym <- ifelse(
    is.na(df$ENSEMBL_clean),
    df$gene_id,
    {
      sy <- suppressMessages(mapIds(org.Hs.eg.db, keys=df$ENSEMBL_clean,
                                    column="SYMBOL", keytype="ENSEMBL", multiVals="first"))
      ifelse(is.na(sy) | sy=="", df$ENSEMBL_clean, sy)
    }
  )
  df
}

res_full <- add_symbols(res_full)
shr_df   <- add_symbols(shr_df)

# ---------- GSEA (fgseaMultilevel by default) ----------
# Build ranking vectors
rank_stat <- res_full %>%
  filter(!is.na(stat), !is.na(gene_sym), gene_sym!="") %>%
  arrange(desc(stat)) %>%
  group_by(gene_sym) %>% slice_head(n=1) %>% ungroup() %>%
  { setNames(.$stat, .$gene_sym) }

rank_lfc <- shr_df %>%
  filter(!is.na(log2FoldChange), !is.na(gene_sym), gene_sym!="") %>%
  arrange(desc(log2FoldChange)) %>%
  group_by(gene_sym) %>% slice_head(n=1) %>% ungroup() %>%
  { setNames(.$log2FoldChange, .$gene_sym) }

# MSigDB gene sets (Symbols) for GSEA
msig_sym <- msigdbr(species = opt$species, collection = opt$msig_collection) %>%
  select(gs_name, gene_symbol) %>%
  split(.$gene_symbol, .$gs_name)

set.seed(1)
gsea_stat <- fgsea(pathways = msig_sym, stats = rank_stat) %>% arrange(padj)
gsea_lfc  <- fgsea(pathways = msig_sym, stats = rank_lfc)  %>% arrange(padj)

# Save GSEA tables
write_csv(gsea_stat, file.path(path_dir, sprintf("%s_fgsea_%s_rankByWaldStat.csv",      opt$project, opt$msig_collection)))
write_csv(gsea_lfc,  file.path(path_dir, sprintf("%s_fgsea_%s_rankByShrunkenLFC.csv",   opt$project, opt$msig_collection)))
saveRDS(gsea_stat,   file.path(path_dir, sprintf("%s_fgsea_%s_rankByWaldStat.rds",      opt$project, opt$msig_collection)))
saveRDS(gsea_lfc,    file.path(path_dir, sprintf("%s_fgsea_%s_rankByShrunkenLFC.rds",   opt$project, opt$msig_collection)))

# GSEA summary dotplots (top 20 by FDR)
mk_gsea_dot <- function(gsea_tbl, fname_title, out_png) {
  df <- gsea_tbl %>%
    mutate(negLog10FDR = -log10(padj),
           direction = ifelse(NES >= 0, "Up", "Down")) %>%
    slice_min(order_by = padj, n = 20) %>%
    arrange(NES) %>%
    mutate(pathway = factor(pathway, levels = pathway))
  p <- ggplot(df, aes(x = negLog10FDR, y = pathway, size = abs(NES), color = direction)) +
    geom_point() +
    labs(x = "-log10(FDR)", y = "", size = "|NES|", color = "", title = fname_title) +
    theme_classic(base_size = 12)
  ggsave(out_png, p, width = 8, height = 6, dpi = 300)
}
mk_gsea_dot(gsea_stat,
            sprintf("GSEA %s (rank: Wald stat) — top 20", opt$msig_collection),
            file.path(plots_dir, sprintf("GSEA_dotplot_%s_rankByWaldStat_top20.png", opt$msig_collection)))
mk_gsea_dot(gsea_lfc,
            sprintf("GSEA %s (rank: shrunken LFC) — top 20", opt$msig_collection),
            file.path(plots_dir, sprintf("GSEA_dotplot_%s_rankByShrunkenLFC_top20.png", opt$msig_collection)))

# A few enrichment curves (top 3 up/down + 4 targeted if present)
plot_es_set <- function(gsea_tbl, stats_vec, set_list, tag) {
  top_up   <- gsea_tbl %>% filter(NES > 0) %>% slice_min(padj, n = 3) %>% pull(pathway)
  top_down <- gsea_tbl %>% filter(NES < 0) %>% slice_min(padj, n = 3) %>% pull(pathway)
  targeted <- c("HALLMARK_INFLAMMATORY_RESPONSE",
                "HALLMARK_INTERFERON_GAMMA_RESPONSE",
                "HALLMARK_OXIDATIVE_PHOSPHORYLATION",
                "HALLMARK_TNFA_SIGNALING_VIA_NFKB")
  to_plot <- unique(c(top_up, top_down, targeted))
  to_plot <- to_plot[to_plot %in% names(set_list)]
  for (pwy in to_plot) {
    png(file.path(plots_dir, paste0("GSEA_ES_",
                                    gsub("[^A-Za-z0-9_]+","_", pwy), "_", tag, ".png")),
        width = 1600, height = 1200, res = 200)
    print(fgsea::plotEnrichment(set_list[[pwy]], stats_vec) + ggtitle(pwy))
    dev.off()
  }
}
plot_es_set(gsea_stat, rank_stat, msig_sym, "WaldStat")
plot_es_set(gsea_lfc,  rank_lfc,  msig_sym, "ShrunkLFC")

# ---------- ORA ----------
# Significant SYMBOLs (from shrunken results) + background
sig_syms <- shr_df %>%
  mutate(padj_safe = ifelse(is.na(padj), 1, padj)) %>%
  filter(padj_safe < opt$padj_cutoff,
         !is.na(log2FoldChange),
         abs(log2FoldChange) >= opt$lfc_cutoff,
         !is.na(gene_sym), gene_sym!="") %>%
  pull(gene_sym) %>% unique()

bg_syms <- res_full %>%
  filter(!is.na(pvalue), !is.na(gene_sym), gene_sym!="") %>%
  pull(gene_sym) %>% unique()

sym2ent_vec <- function(x) {
  if (length(x) == 0) return(character(0))
  xu <- unique(x)
  m  <- suppressWarnings(suppressMessages(
    bitr(xu, fromType="SYMBOL", toType="ENTREZID", OrgDb=org.Hs.eg.db)
  ))
  if (!is.null(m) && nrow(m) > 0) {
    rate <- round(100 * length(unique(m$SYMBOL)) / length(xu), 2)
    message("SYMBOL->ENTREZ mapped: ", rate, "% (", length(unique(m$SYMBOL)), "/", length(xu), ")")
    return(unique(m$ENTREZID))
  } else {
    message("SYMBOL->ENTREZ mapped: 0% (0/", length(xu), ")")
    return(character(0))
  }
}
entrez_sig <- sym2ent_vec(sig_syms)
entrez_bg  <- sym2ent_vec(bg_syms)

if (length(entrez_sig) == 0) {
  message("ORA skipped: no padj<=", opt$padj_cutoff,
          " & |LFC|>=", opt$lfc_cutoff, " genes (or none mapped to ENTREZ).")
} else {
  # GO BP
  ego <- enrichGO(gene          = entrez_sig,
                  universe      = entrez_bg,
                  OrgDb         = org.Hs.eg.db,
                  keyType       = "ENTREZID",
                  ont           = "BP",
                  pAdjustMethod = "BH",
                  qvalueCutoff  = 0.05,
                  readable      = TRUE)
  # KEGG
  ekegg <- enrichKEGG(gene          = entrez_sig,
                      universe      = entrez_bg,
                      organism      = "hsa",
                      pAdjustMethod = "BH",
                      qvalueCutoff  = 0.05)

  # Hallmark ORA via msigdbr ENTREZ directly (no joins)
  msig_h_term2gene <- msigdbr(species = opt$species, collection = "H") %>%
    select(gs_name, entrez_gene) %>%
    filter(!is.na(entrez_gene)) %>%
    transmute(gs_name, ENTREZID = as.character(entrez_gene)) %>%
    distinct()
  eh <- enricher(gene = entrez_sig, universe = entrez_bg, TERM2GENE = msig_h_term2gene)

  # Save tables + RDS
  write_csv(as.data.frame(ego),   file.path(path_dir, sprintf("%s_ORA_GO_BP.csv",   opt$project)))
  write_csv(as.data.frame(ekegg), file.path(path_dir, sprintf("%s_ORA_KEGG.csv",    opt$project)))
  write_csv(as.data.frame(eh),    file.path(path_dir, sprintf("%s_ORA_Hallmark.csv",opt$project)))
  saveRDS(ego,   file.path(path_dir, sprintf("%s_ORA_GO_BP.rds",    opt$project)))
  saveRDS(ekegg, file.path(path_dir, sprintf("%s_ORA_KEGG.rds",     opt$project)))
  saveRDS(eh,    file.path(path_dir, sprintf("%s_ORA_Hallmark.rds", opt$project)))

  # Dotplots
  if (nrow(as.data.frame(ego)) > 0) {
    png(file.path(plots_dir, "ORA_GO_BP_dotplot_top20.png"), width=1600, height=1200, res=200)
    print(dotplot(ego, showCategory=20) + ggtitle("GO BP (ORA) — top 20"))
    dev.off()
  }
  if (nrow(as.data.frame(ekegg)) > 0) {
    png(file.path(plots_dir, "ORA_KEGG_dotplot_top20.png"), width=1600, height=1200, res=200)
    print(dotplot(ekegg, showCategory=20) + ggtitle("KEGG (ORA) — top 20"))
    dev.off()
  }
  if (nrow(as.data.frame(eh)) > 0) {
    png(file.path(plots_dir, "ORA_Hallmark_dotplot_top20.png"), width=1600, height=1200, res=200)
    print(dotplot(eh, showCategory=20) + ggtitle("MSigDB Hallmark (ORA) — top 20"))
    dev.off()
  }

  # Optional cnetplot (nice-to-have)
  if (isTRUE(opt$make_cnet) && nrow(as.data.frame(ego)) > 0) {
    fc_vec <- shr_df %>% filter(!is.na(gene_sym)) %>%
      distinct(gene_sym, .keep_all = TRUE) %>%
      { stats::setNames(.$log2FoldChange, .$gene_sym) }
    tryCatch({
      png(file.path(plots_dir, "ORA_GO_BP_cnetplot_top10.png"), width=2000, height=1400, res=200)
      print(cnetplot(ego, showCategory = 10, foldChange = fc_vec))
      dev.off()
    }, error = function(e) {
      message("cnetplot skipped: ", conditionMessage(e))
    })
  }
}

# ---------- Echo ----------
message("\n[Done] Outputs written to: ", path_dir,
        "\n  - GSEA CSV/RDS: ", opt$msig_collection,
        "\n  - ORA CSV/RDS: GO BP / KEGG / Hallmark",
        "\n  - Plots in: ", plots_dir, "\n")
