#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(optparse)
  library(readr); library(readxl); library(dplyr); library(tidyr)
  library(janitor); library(stringr)
  library(DESeq2); library(limma)
  library(ggplot2); library(ggrepel)
  library(matrixStats)
  library(AnnotationDbi)
  if (!requireNamespace("org.Hs.eg.db", quietly = TRUE)) {
    install.packages("BiocManager"); BiocManager::install("org.Hs.eg.db", ask = FALSE, update = FALSE)
  }
  library(org.Hs.eg.db)
})

set.seed(1)

# ---------- helpers ----------
theme_pub <- function() {
  theme_classic(base_size = 12) +
    theme(panel.grid = element_blank(),
          axis.ticks = element_line(),
          axis.text  = element_text(color = "black"),
          axis.title = element_text(color = "black"),
          legend.key = element_blank())
}

clean_fc_names <- function(x) {
  x <- basename(x)
  x <- sub("\\.bam$", "", x, ignore.case = TRUE)
  x <- sub("\\.sam$", "", x, ignore.case = TRUE)
  x <- sub("_Aligned(\\.sortedByCoord)?\\.out$", "", x, ignore.case = TRUE)
  x <- sub("\\.sorted$", "", x, ignore.case = TRUE)
  x <- sub("\\.featureCounts.*$", "", x, ignore.case = TRUE)
  x <- sub("\\.counts?$", "", x, ignore.case = TRUE)
  x <- gsub("\\s+", "", x)
  x
}

add_symbols <- function(df) {
  df$Geneid <- rownames(df)
  looks_ens <- grepl("^ENSG\\d+", df$Geneid)
  df$ENSEMBL_clean <- ifelse(looks_ens, sub("\\.\\d+$", "", df$Geneid), NA_character_)
  if (any(looks_ens)) {
    df$symbol   <- AnnotationDbi::mapIds(org.Hs.eg.db, keys = df$ENSEMBL_clean,
                                         column = "SYMBOL", keytype = "ENSEMBL", multiVals = "first")
    df$genename <- AnnotationDbi::mapIds(org.Hs.eg.db, keys = df$ENSEMBL_clean,
                                         column = "GENENAME", keytype = "ENSEMBL", multiVals = "first")
    df$gene <- ifelse(is.na(df$symbol) | df$symbol == "", df$ENSEMBL_clean, df$symbol)
  } else {
    df$gene <- df$Geneid
  }
  dplyr::relocate(df, gene, .before = 1)
}

annotate_biotype <- function(resDF, raw_counts_df) {
  looks_ens <- grepl("^ENSG\\d+", resDF$Geneid)
  resDF$ENSEMBL_clean <- ifelse(looks_ens, sub("\\.\\d+$", "", resDF$Geneid), NA_character_)
  biotype_col <- grep("biotype|type", colnames(raw_counts_df), ignore.case = TRUE, value = TRUE)
  if (length(biotype_col) > 0) {
    tmp <- raw_counts_df[, c("Geneid", biotype_col[1])]
    names(tmp) <- c("Geneid", "gene_biotype")
    tmp$ENSEMBL_clean <- ifelse(grepl("^ENSG\\d+", tmp$Geneid), sub("\\.\\d+$","", tmp$Geneid), NA_character_)
    bt <- tmp[, c("ENSEMBL_clean", "gene_biotype")]
    dplyr::left_join(resDF, bt, by = "ENSEMBL_clean")
  } else if (any(looks_ens)) {
    if (!requireNamespace("EnsDb.Hsapiens.v86", quietly = TRUE)) {
      install.packages("BiocManager"); BiocManager::install("EnsDb.Hsapiens.v86", ask = FALSE, update = FALSE)
    }
    suppressPackageStartupMessages(library(EnsDb.Hsapiens.v86))
    edb <- EnsDb.Hsapiens.v86
    bt <- as.data.frame(genes(edb, columns = c("gene_id","gene_biotype")))[, c("gene_id","gene_biotype")]
    names(bt) <- c("ENSEMBL_clean","gene_biotype")
    dplyr::left_join(resDF, bt, by = "ENSEMBL_clean")
  } else {
    resDF$gene_biotype <- NA_character_; resDF
  }
}

get_len_kb <- function(raw_tbl) {
  if ("Length" %in% names(raw_tbl)) {
    len <- setNames(as.numeric(raw_tbl$Length)/1000, raw_tbl$Geneid)
    return(len)
  }
  NULL
}

make_TPM_log2 <- function(cts, len_kb_named) {
  common <- base::intersect(rownames(cts), names(len_kb_named))
  if (length(common) < 1000) return(NULL)
  cts2 <- cts[common, , drop = FALSE]
  len2 <- len_kb_named[common]
  rpk  <- sweep(cts2, 1, len2, "/")
  tpm  <- 1e6 * sweep(rpk, 2, colSums(rpk), "/")
  log2(tpm + 1)
}

maybe_write_csv  <- function(df, path, save_files)  { if (isTRUE(save_files)) readr::write_csv(df, path) }
maybe_ggsave     <- function(plot, path, save_files, width=7, height=6, dpi=300) { if (isTRUE(save_files)) ggsave(path, plot, width=width, height=height, dpi=dpi) }
maybe_dir_create <- function(path, save_files) { if (isTRUE(save_files)) dir.create(path, showWarnings = FALSE, recursive = TRUE) }

summarize_res_tbl <- function(df, alpha = 0.05, lfc_cut = 1) {
  df2 <- df %>% dplyr::mutate(padj_safe = ifelse(is.na(padj), 1, padj))
  tested <- sum(!is.na(df2$padj_safe))
  sig    <- df2 %>% dplyr::filter(padj_safe < alpha)
  up     <- sum(sig$log2FoldChange > 0, na.rm = TRUE)
  down   <- sum(sig$log2FoldChange < 0, na.rm = TRUE)
  n_abs1 <- sig %>% dplyr::filter(abs(log2FoldChange) >= lfc_cut) %>% nrow()
  med_abs_lfc <- sig %>% dplyr::summarize(median_abs_LFC = median(abs(log2FoldChange), na.rm = TRUE)) %>% dplyr::pull()
  tibble::tibble(
    tested_genes      = tested,
    significant_FDR   = nrow(sig),
    upregulated       = up,
    downregulated     = down,
    sig_absLFC_ge_1   = n_abs1,
    median_absLFC_sig = round(med_abs_lfc, 3)
  )
}

# ---------- CLI ----------
option_list <- list(
  make_option(c("-c","--counts"),   type="character", help="Path to featureCounts TSV (must include Geneid)"),
  make_option(c("-m","--metadata"), type="character", help="Path to metadata (xlsx or csv/tsv)"),
  make_option(c("-o","--outdir"),   type="character", help="Root output directory (absolute path)"),
  make_option(c("--project"),       type="character", default="PROJECT", help="Project label for output paths"),
  make_option(c("--genotype-col"),  type="character", default="genotype", help="Metadata column with genotype"),
  make_option(c("--sample-col"),    type="character", default="sample",   help="Metadata column with sample IDs"),
  make_option(c("--batch-col"),     type="character", default="handler",  help="Metadata column used as batch"),
  make_option(c("--case-geno"),     type="character", help="Genotype value for CASE group (e.g. 'A/A')"),
  make_option(c("--ctrl-geno"),     type="character", help="Genotype value for CONTROL group (e.g. 'G/G')"),
  make_option(c("--alpha"),         type="double", default=0.05, help="FDR cutoff"),
  make_option(c("--lfc-cut"),       type="double", default=1.0,  help="abs(LFC) cutoff for summary"),
  make_option(c("--topn"),          type="integer", default=100,  help="Top-N labels per direction on volcano"),
  make_option(c("--no-tpm-pca"),    action="store_true", default=FALSE, help="Disable TPM PCA"),
  make_option(c("--no-varpart"),    action="store_true", default=FALSE, help="Disable variancePartition"),
  make_option(c("--dryrun"),        action="store_true", default=FALSE, help="Don’t write files (preview only)")
)
opt <- parse_args(OptionParser(option_list=option_list))

# sanity checks
required <- c("counts","metadata","outdir","case-geno","ctrl-geno")
miss <- required[!nzchar(sapply(required, function(k) opt[[gsub("-","_",k)]] %||% ""))]
if (length(miss) > 0) stop("Missing required flags: ", paste("--", miss, collapse=" "))

# load metadata (xlsx or csv/tsv)
meta_path <- opt$metadata
if (grepl("\\.xlsx?$", meta_path, ignore.case=TRUE)) {
  meta <- readxl::read_excel(meta_path)
} else {
  meta <- readr::read_delim(meta_path, delim = ifelse(grepl("\\.tsv$", meta_path, ignore.case=TRUE), "\t", ","), show_col_types = FALSE)
}
meta <- meta %>% janitor::clean_names()

# map columns
sample_col <- opt$`sample-col`; geno_col <- opt$`genotype-col`; batch_col <- opt$`batch-col`
stopifnot(sample_col %in% names(meta), geno_col %in% names(meta), batch_col %in% names(meta))

# output dirs
root_output_dir <- normalizePath(opt$outdir, mustWork = FALSE)
run_tag <- format(Sys.time(), "%Y%m%d-%H%M%S")
scratch_dir <- file.path(root_output_dir, paste0("run_", run_tag))
project <- opt$project
qc_dir  <- file.path(scratch_dir, project, "QC")
out_dir <- file.path(scratch_dir, project, "DEG_handler_adjusted")
if (!opt$dryrun) { dir.create(qc_dir, recursive=TRUE, showWarnings=FALSE); dir.create(out_dir, recursive=TRUE, showWarnings=FALSE) }

# load counts
raw_full <- readr::read_tsv(opt$counts, comment = "#", show_col_types = FALSE)

# standardize Geneid col
stopifnot(any(grepl("^geneid$", names(raw_full), ignore.case = TRUE)))
gene_col <- names(raw_full)[grepl("^geneid$", names(raw_full), ignore.case = TRUE)][1]
raw_full <- raw_full %>% dplyr::rename(Geneid = dplyr::all_of(gene_col))

# harmonize sample names in counts vs metadata
raw_sample_cols <- base::setdiff(names(raw_full), c("Geneid","Chr","Start","End","Strand","Length"))
raw_clean <- clean_fc_names(raw_sample_cols) %>% stringr::str_replace_all("_", "-") %>% toupper()
meta[[sample_col]] <- meta[[sample_col]] %>% as.character() %>% stringr::str_trim() %>% stringr::str_replace_all("_", "-") %>% toupper()
names(raw_full)[match(raw_sample_cols, names(raw_full))] <- raw_clean

# factors & group mapping
case_geno  <- opt$`case-geno`; ctrl_geno <- opt$`ctrl-geno`
meta <- meta %>%
  dplyr::mutate(
    sample  = as.character(.data[[sample_col]]),
    handler = factor(.data[[batch_col]]),
    group   = factor(
      ifelse(.data[[geno_col]] == case_geno, "case",
             ifelse(.data[[geno_col]] == ctrl_geno, "control", NA)),
      levels = c("control","case")
    )
  )
stopifnot(!anyDuplicated(meta$sample))
stopifnot(all(meta$sample %in% names(raw_full)))

# counts matrix
raw <- raw_full[, c("Geneid", meta$sample)]
cts <- as.matrix(raw[, -1, drop = FALSE]); mode(cts) <- "integer"
rownames(cts) <- raw$Geneid
stopifnot(nrow(cts) > 1000, ncol(cts) == nrow(meta))

# lib size plot
df_lib <- data.frame(sample = colnames(cts), libs = colSums(cts)) %>%
  dplyr::left_join(meta, by = dplyr::join_by(sample))
p_lib <- ggplot(df_lib, aes(reorder(sample, libs), libs/1e6, fill = group)) +
  geom_col() + coord_flip() + theme_pub() +
  labs(title = paste0(project, ": Library size (millions of reads)"), x = "", y = "Millions")
print(p_lib); maybe_ggsave(p_lib, file.path(qc_dir, "01_library_sizes.png"), !opt$dryrun, width=6, height=5, dpi=300)

# blind VST PCA
coldata0 <- meta %>% tibble::column_to_rownames("sample")
dds_blind <- DESeqDataSetFromMatrix(cts, coldata0, design = ~ 1)
keep0 <- rowSums(counts(dds_blind) >= 10) >= 2
dds_blind <- dds_blind[keep0, ]
vsd_blind <- vst(dds_blind, blind = TRUE)
v_blind   <- SummarizedExperiment::assay(vsd_blind)
pcs_b     <- prcomp(t(v_blind), scale. = TRUE)
varexp_b  <- round(100 * pcs_b$sdev^2 / sum(pcs_b$sdev^2), 1)
pcdf_b    <- data.frame(pcs_b$x[, 1:2], sample = rownames(pcs_b$x),
                        group = coldata0$group, handler = coldata0$handler)
p_pca_blind <- ggplot(pcdf_b, aes(PC1, PC2, color = group, shape = handler, label = sample)) +
  geom_point(size = 3) +
  ggrepel::geom_text_repel(size = 3, box.padding = 0.3, max.overlaps = 50) +
  labs(title = paste0(project, ": PCA (blind VST; BEFORE correction)"),
       x = paste0("PC1 (", varexp_b[1], "%)"),
       y = paste0("PC2 (", varexp_b[2], "%)")) +
  theme_pub()
print(p_pca_blind)
maybe_ggsave(p_pca_blind, file.path(qc_dir, "02_pca_blind_before.png"), !opt$dryrun, width=7, height=6, dpi=600)

# optional TPM PCA
if (!opt$`no-tpm-pca`) {
  logtpm_filt <- NULL
  len_kb <- get_len_kb(raw_full)
  logtpm <- if (!is.null(len_kb)) make_TPM_log2(cts, setNames(len_kb, names(len_kb))) else NULL
  if (!is.null(logtpm)) {
    keep_nonzero <- rowSums(2^logtpm - 1) > 0
    keep_var     <- rowVars(logtpm) > 0
    logtpm_filt  <- logtpm[keep_nonzero & keep_var, , drop = FALSE]
    logtpm_filt  <- logtpm_filt[, rownames(coldata0)]
    pcs_tpm    <- prcomp(t(logtpm_filt), scale. = TRUE)
    varexp_tpm <- round(100 * pcs_tpm$sdev^2 / sum(pcs_tpm$sdev^2), 1)
    pcdf_tpm   <- data.frame(pcs_tpm$x[, 1:2], sample = colnames(logtpm_filt),
                             group = coldata0$group, handler = coldata0$handler)
    p_pca_tpm <- ggplot(pcdf_tpm, aes(PC1, PC2, color = group, shape = handler, label = sample)) +
      geom_point(size = 3) +
      ggrepel::geom_text_repel(size = 3, box.padding = 0.3, max.overlaps = 50) +
      labs(title = paste0(project, ": PCA on log2 TPM (blind; BEFORE correction)"),
           x = paste0("PC1 (", varexp_tpm[1], "%)"),
           y = paste0("PC2 (", varexp_tpm[2], "%)")) +
      theme_pub()
    print(p_pca_tpm)
    maybe_ggsave(p_pca_tpm, file.path(qc_dir, "02b_pca_TPM_blind.png"), !opt$dryrun, width=8, height=7, dpi=600)
  } else {
    message(project, ": No usable Length column found; skipping TPM PCA.")
  }
}

# optional variancePartition
if (!opt$`no-varpart`) {
  if (!requireNamespace("variancePartition", quietly = TRUE)) {
    install.packages("BiocManager"); BiocManager::install("variancePartition", ask = FALSE, update = FALSE)
  }
  library(variancePartition)
  expr_b <- v_blind
  md_b   <- as.data.frame(SummarizedExperiment::colData(dds_blind))
  form_b <- ~ (1|group) + (1|handler)
  vp_b   <- fitExtractVarPartModel(expr_b, form_b, md_b)
  if (!opt$dryrun) {
    png(file.path(qc_dir, "02c_variancePartition_blind.png"), width=2400, height=1800, res=300)
    print(plotVarPart(vp_b)); dev.off()
  }
}

# corrected model & PCA
dds <- DESeqDataSetFromMatrix(cts, coldata0, design = ~ handler + group)
keep <- rowSums(counts(dds) >= 10) >= 2
dds  <- dds[keep, ]
X <- model.matrix(~ handler + group, data = as.data.frame(SummarizedExperiment::colData(dds)))
if (qr(X)$rank < ncol(X)) warning("Design is not full rank — check handler/group balance.")

vsd <- vst(dds, blind = FALSE)
v   <- SummarizedExperiment::assay(vsd)
v_vis <- limma::removeBatchEffect(
  v,
  batch = SummarizedExperiment::colData(dds)$handler,
  covariates = model.matrix(~ group, data = as.data.frame(SummarizedExperiment::colData(dds)))[, -1]
)
pcs <- prcomp(t(v_vis), scale. = TRUE)
varexp <- round(100 * pcs$sdev^2 / sum(pcs$sdev^2), 1)
pc_df <- data.frame(pcs$x[, 1:2], sample = rownames(pcs$x),
                    group = SummarizedExperiment::colData(dds)$group,
                    handler = SummarizedExperiment::colData(dds)$handler)
p_pca_after <- ggplot(pc_df, aes(PC1, PC2, color = group, shape = handler, label = sample)) +
  geom_point(size = 3) +
  ggrepel::geom_text_repel(size = 3, box.padding = 0.3, max.overlaps = 50) +
  labs(title = paste0(project, ": PCA (VST; handler removed for viz, group preserved)"),
       x = paste0("PC1 (", varexp[1], "%)"),
       y = paste0("PC2 (", varexp[2], "%)")) +
  theme_pub()
print(p_pca_after)
maybe_ggsave(p_pca_after, file.path(qc_dir, "03_pca_VST_after.png"), !opt$dryrun, width=8, height=7, dpi=600)

# DESeq2 fit
dds <- DESeq(dds)

if (!opt$dryrun) {
  png(file.path(out_dir, "01_dispersions.png"), width = 1600, height = 1200, res = 200)
  plotDispEsts(dds); dev.off()
}

# results (raw + apeglm-shrunken)
res_raw <- results(dds, contrast = c("group", "case", "control"))
res_raw <- res_raw[order(res_raw$padj, res_raw$pvalue), ]

if (!opt$dryrun) {
  png(file.path(out_dir, "02_pvalue_histogram.png"), width = 1600, height = 1200, res = 200)
  hist(res_raw$pvalue, breaks = 50, main = "P-value histogram", xlab = "p"); dev.off()
}

if (!requireNamespace("apeglm", quietly = TRUE)) {
  install.packages("BiocManager"); BiocManager::install("apeglm", ask = FALSE, update = FALSE)
}
coef_name <- grep("group.*case.*vs.*control", resultsNames(dds), value = TRUE)[1]
stopifnot(!is.na(coef_name))
res_shr <- lfcShrink(dds, coef = coef_name, type = "apeglm")
res_shr <- res_shr[order(res_shr$padj, res_shr$pvalue), ]

if (!opt$dryrun) {
  png(file.path(out_dir, "03_MA_raw.png"), width = 1600, height = 1200, res = 200)
  plotMA(res_raw, main = "MA plot (raw LFC)", ylim = c(-6, 6)); dev.off()
  png(file.path(out_dir, "04_MA_shrunken.png"), width = 1600, height = 1200, res = 200)
  plotMA(res_shr, main = "MA plot (apeglm-shrunken LFC)", ylim = c(-6, 6)); dev.off()
}

# symbols + biotype
resDF_raw <- add_symbols(as.data.frame(res_raw))
resDF_shr <- add_symbols(as.data.frame(res_shr))
resDF_shr <- annotate_biotype(resDF_shr, raw_full)

maybe_write_csv(resDF_raw, file.path(out_dir, "05_DESeq2_results_raw_byPadj.csv"), !opt$dryrun)
maybe_write_csv(resDF_raw %>% dplyr::arrange(dplyr::desc(log2FoldChange)),
                file.path(out_dir, "05_DESeq2_results_raw_byLFC.csv"), !opt$dryrun)
maybe_write_csv(resDF_shr %>% dplyr::arrange(padj, pvalue),
                file.path(out_dir, "06_DESeq2_results_LFCshrunken_all_byPadj.csv"), !opt$dryrun)
maybe_write_csv(resDF_shr %>% dplyr::arrange(dplyr::desc(log2FoldChange)),
                file.path(out_dir, "06_DESeq2_results_LFCshrunken_all_byLFC.csv"), !opt$dryrun)

# protein-coding + volcano
res_pc_shr <- resDF_shr %>% dplyr::filter(!is.na(gene_biotype) & gene_biotype == "protein_coding")
maybe_write_csv(res_pc_shr %>% dplyr::arrange(padj, pvalue),
                file.path(out_dir, "07_DESeq2_results_LFCshrunken_proteinCoding_byPadj.csv"), !opt$dryrun)
maybe_write_csv(res_pc_shr %>% dplyr::arrange(dplyr::desc(log2FoldChange)),
                file.path(out_dir, "07_DESeq2_results_LFCshrunken_proteinCoding_byLFC.csv"), !opt$dryrun)

volc_df <- res_pc_shr %>% dplyr::mutate(
  padj_plot = ifelse(is.na(padj), 1, padj),
  log10padj = -log10(pmax(padj_plot, .Machine$double.xmin)),
  direction = dplyr::case_when(
    !is.na(padj) & padj < opt$alpha & log2FoldChange > 0 ~ "Up (FDR<0.05)",
    !is.na(padj) & padj < opt$alpha & log2FoldChange < 0 ~ "Down (FDR<0.05)",
    TRUE ~ "NS"
  )
)
sig <- dplyr::filter(volc_df, direction != "NS")
label_topN <- function(N) {
  dplyr::bind_rows(
    sig %>% dplyr::arrange(dplyr::desc(log2FoldChange)) %>% dplyr::slice_head(n = N),
    sig %>% dplyr::arrange(log2FoldChange)              %>% dplyr::slice_head(n = N)
  )
}
labset <- if (nrow(sig) > 0) label_topN( min(opt$topn, max(10, floor(nrow(sig)/10))) ) else sig

p_volc <- ggplot(volc_df, aes(x = log2FoldChange, y = log10padj)) +
  geom_point(aes(color = direction), alpha = 0.7, size = 1.3) +
  scale_color_manual(values = c("Down (FDR<0.05)"="#3B82F6","NS"="grey70","Up (FDR<0.05)"="#EF4444")) +
  ggrepel::geom_text_repel(data = labset, aes(label = gene),
                           size = 2.6, max.overlaps = Inf, box.padding = 0.35, point.padding = 0.15) +
  theme_bw() +
  labs(title = paste0(project, ": Volcano (protein-coding; shrunken)"),
       x = "log2 fold change (case vs control)", y = "-log10(FDR)", color = "")
print(p_volc)
maybe_ggsave(p_volc, file.path(out_dir, "08_volcano_shrunken_proteinCoding.png"), !opt$dryrun, width=7.5, height=6, dpi=600)

sum_all <- summarize_res_tbl(resDF_shr, alpha = opt$alpha, lfc_cut = opt$`lfc-cut`)
res_pc_shr_nullsafe <- if (nrow(res_pc_shr) == 0) resDF_shr[0,] else res_pc_shr
sum_pc  <- summarize_res_tbl(res_pc_shr_nullsafe, alpha = opt$alpha, lfc_cut = opt$`lfc-cut`)
maybe_write_csv(sum_all, file.path(out_dir, "09_SUMMARY_all_shrunken.csv"), !opt$dryrun)
maybe_write_csv(sum_pc,  file.path(out_dir, "09_SUMMARY_proteinCoding_shrunken.csv"), !opt$dryrun)

sig_tbl <- res_pc_shr %>% dplyr::mutate(padj_safe = ifelse(is.na(padj), 1, padj)) %>%
  dplyr::filter(padj_safe < opt$alpha, !is.na(log2FoldChange))
if (nrow(sig_tbl) > 0) {
  top_up <- sig_tbl %>% dplyr::filter(log2FoldChange > 0) %>%
    dplyr::arrange(dplyr::desc(log2FoldChange)) %>% dplyr::slice_head(n = opt$topn) %>%
    dplyr::mutate(direction_set = "UP_byLFC", rank_byLFC = dplyr::row_number()) %>%
    dplyr::select(direction_set, rank_byLFC, gene, Geneid, log2FoldChange, lfcSE, stat, pvalue, padj, baseMean)
  top_dn <- sig_tbl %>% dplyr::filter(log2FoldChange < 0) %>%
    dplyr::arrange(log2FoldChange) %>% dplyr::slice_head(n = opt$topn) %>%
    dplyr::mutate(direction_set = "DOWN_byLFC", rank_byLFC = dplyr::row_number()) %>%
    dplyr::select(direction_set, rank_byLFC, gene, Geneid, log2FoldChange, lfcSE, stat, pvalue, padj, baseMean)
  maybe_write_csv(top_up, file.path(out_dir, sprintf("10_DEG_top%d_UP_byLFC.csv", opt$topn)), !opt$dryrun)
  maybe_write_csv(top_dn, file.path(out_dir, sprintf("10_DEG_top%d_DOWN_byLFC.csv", opt$topn)), !opt$dryrun)
}

norm <- counts(dds, normalized = TRUE)
maybe_write_csv(as.data.frame(norm) %>% tibble::rownames_to_column("Geneid"),
                file.path(out_dir, "11_normalized_counts.csv"), !opt$dryrun)

if (!opt$dryrun) writeLines(capture.output(sessionInfo()), con = file.path(out_dir, "sessionInfo.txt"))

message("\nComplete. Outputs in:\n  ", file.path(scratch_dir, project), "\n")
