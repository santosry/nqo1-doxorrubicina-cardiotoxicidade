#!/usr/bin/env Rscript
# ==============================================================================
# ampk_de_analysis.R — Pipeline UPSTREAM (PASSO 1 + PASSO 2)
# ==============================================================================
# Narrativa oficial do estudo (3 passos):
#   PASSO 1 — Analise exploratoria da via AMPK (KEGG hsa04152) na DCM (GSE116250).
#             Resultado: NAO significativo (GSEA NES = -1.19; padj = 0.419).
#   PASSO 2 — Expressao diferencial do TRANSCRIPTOMA COMPLETO -> selecao do gene
#             NQO1 (downregulado, log2FC = -1.686; padj = 6.86e-13).
#   PASSO 3 — (no script.R) Consulta DGB -> Doxorrubicina (q = 2.21e-28).
#
# Este script executa o PASSO 1 e o PASSO 2 e salva TODOS os outputs em
# output/tables/ampk/:
#   - DEG_full_table.csv      (DE do transcriptoma completo)
#   - NQO1_results.csv        (linha do gene NQO1)
#   - DEG_AMPK_pathway.csv    (DE restrito a via AMPK hsa04152)
#   - GSEA_AMPK_results.csv   (GSEA da via AMPK)
#   - GSVA_AMPK_results.csv   (escore GSVA da via AMPK por amostra)
#
# CORRECOES DE AUDITORIA APLICADAS:
#   * GSE116250 e RNA-seq (NAO microarray). Matriz RPKM (Cufflinks).
#   * NQO1 NAO pertence a via AMPK; foi identificado no transcriptoma completo.
#   * NQO1 esta DOWNREGULADO na DCM (nao upregulado).
#   * Normalizacao: log2(RPKM + 1) + limma (melhoria 6).
#   * Ajuste por covariadas sex + idade (melhoria 7).
#   * ICM excluida do contraste DCM vs CTRL (melhoria 4, documentada).
#   * Rotulo do grupo controle padronizado como CTRL (melhoria 5).
#   * Leitura com data.table::fread para arquivos .gz (melhoria 8).
# ==============================================================================

suppressPackageStartupMessages({
  library(data.table)
  library(readr)
  library(dplyr)
  library(tidyr)
  library(limma)
})

# -----------------------------------------------------------------------------
# Helpers
# -----------------------------------------------------------------------------
`%+%` <- function(a, b) paste0(a, b)

log_dir <- file.path("output", "audit", "logs")
dir.create(log_dir, recursive = TRUE, showWarnings = FALSE)
logf <- function(level, msg) {
  line <- sprintf("[%s] %-5s | %s", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), level, msg)
  message(line)
  cat(line, "\n", file = file.path(log_dir, "ampk_de_analysis.log"), append = TRUE)
}

OUT <- file.path("output", "tables", "ampk")
dir.create(OUT, recursive = TRUE, showWarnings = FALSE)

# -----------------------------------------------------------------------------
# 0. Leitura dos dados brutos
# -----------------------------------------------------------------------------
rpkm_path <- "data/raw/GSE116250_rpkm.txt.gz"
pheno_path <- file.path(OUT, "pheno_data_raw.csv")
ampk_genes_path <- file.path(OUT, "AMPK_pathway_genes_KEGG.csv")

logf("INFO", "Lendo matriz RPKM (RNA-seq): " %+% rpkm_path)
rpkm <- data.table::fread(rpkm_path)

gene_names <- rpkm[["Common_name"]]
expr <- as.matrix(rpkm[, !c("Gene", "Common_name"), with = FALSE])
rownames(expr) <- gene_names
storage.mode(expr) <- "numeric"
logf("INFO", "Matriz: " %+% nrow(expr) %+% " genes x " %+% ncol(expr) %+% " amostras")

# Metadados (derivados da series matrix na analise original)
pheno <- readr::read_csv(pheno_path, show_col_types = FALSE) |>
  dplyr::select(sample = title, disease = `disease:ch1`, sex = `Sex:ch1`, age = `age:ch1`) |>
  dplyr::mutate(
    group = dplyr::case_when(
      disease == "non-failing" ~ "CTRL",            # padronizacao do rotulo controle
      disease == "dilated cardiomyopathy" ~ "DCM",
      disease == "ischemic cardiomyopathy" ~ "ICM",
      TRUE ~ NA_character_
    )
  ) |>
  dplyr::filter(sample %in% colnames(expr), !is.na(group))

logf("INFO", "Amostras com metadado: " %+% nrow(pheno) %+%
       " (CTRL=" %+% sum(pheno$group == "CTRL") %+%
       ", DCM=" %+% sum(pheno$group == "DCM") %+%
       ", ICM=" %+% sum(pheno$group == "ICM") %+% ")")

# -----------------------------------------------------------------------------
# Normalizacao (melhoria 6): log2(RPKM + 1)
# -----------------------------------------------------------------------------
expr_log <- log2(expr + 1)
logf("INFO", "Normalizacao: log2(RPKM + 1) — RNA-seq (nao microarray)")

# -----------------------------------------------------------------------------
# Expressao diferencial do TRANSCRIPTOMA COMPLETO — DCM vs CTRL
# (ICM excluida do contraste; ajuste por sex + idade)
# -----------------------------------------------------------------------------
keep <- pheno$sample[pheno$group %in% c("DCM", "CTRL")]
expr_de <- expr_log[, keep, drop = FALSE]

# Filtro de genes pouco expressos (>= 25% das amostras do contraste)
keep_genes <- rowSums(expr[, keep, drop = FALSE] > 1) >= (0.25 * length(keep))
expr_de <- expr_de[keep_genes, , drop = FALSE]
logf("INFO", "Genes apos filtro de expressao: " %+% nrow(expr_de))

design_df <- pheno |>
  dplyr::filter(group %in% c("DCM", "CTRL")) |>
  dplyr::mutate(group = factor(group, levels = c("CTRL", "DCM"))) |>
  dplyr::arrange(match(sample, colnames(expr_de)))

design <- model.matrix(~ 0 + group + sex + age, data = design_df)
colnames(design) <- sub("^group", "", colnames(design))
logf("INFO", "Design (melhoria 7): " %+% paste(colnames(design), collapse = " + "))

fit <- limma::lmFit(expr_de, design)
fit <- limma::contrasts.fit(fit, limma::makeContrasts(DCM - CTRL, levels = design))
fit <- limma::eBayes(fit, trend = TRUE)

tt <- limma::topTable(fit, number = Inf, sort.by = "none")
deg_full <- data.frame(gene = rownames(tt), tt, row.names = NULL, stringsAsFactors = FALSE) |>
  dplyr::mutate(
    neg_log10_padj = -log10(pmax(adj.P.Val, 1e-300)),
    direction = ifelse(logFC > 0, "up", "down"),
    significant = adj.P.Val < 0.05
  ) |>
  dplyr::select(gene, log2FoldChange = logFC, AveExpr, t, pvalue = P.Value,
                padj = adj.P.Val, B, direction, significant, neg_log10_padj)

readr::write_csv(deg_full, file.path(OUT, "DEG_full_table.csv"))
logf("INFO", "DEG_full_table.csv: " %+% nrow(deg_full) %+% " genes")

# -----------------------------------------------------------------------------
# PASSO 2 — Selecao do gene NQO1 (transcriptoma completo)
# -----------------------------------------------------------------------------
nqo1 <- deg_full |> dplyr::filter(gene == "NQO1")
if (nrow(nqo1) == 0) {
  logf("WARN", "NQO1 nao encontrado no DEG.")
} else {
  logf("INFO", sprintf("PASSO 2 | NQO1: log2FC = %.3f | padj = %.3e | %s",
                       nqo1$log2FoldChange[1], nqo1$padj[1], nqo1$direction[1]))
}
readr::write_csv(nqo1, file.path(OUT, "NQO1_results.csv"))
logf("INFO", "NQO1_results.csv gravado (downregulado na DCM; NAO pertence a via AMPK)")

# -----------------------------------------------------------------------------
# PASSO 1 — Caracterizacao da via AMPK (hsa04152, KEGG)
# -----------------------------------------------------------------------------
ampk_genes <- readr::read_csv(ampk_genes_path, show_col_types = FALSE)
ampk_symbols <- unique(ampk_genes$symbol)
logf("INFO", "PASSO 1 | Via AMPK (hsa04152): " %+% length(ampk_symbols) %+% " genes")

deg_ampk <- deg_full |>
  dplyr::filter(gene %in% ampk_symbols) |>
  dplyr::left_join(
    ampk_genes |> dplyr::select(symbol, description) |> dplyr::distinct(symbol, .keep_all = TRUE),
    by = c("gene" = "symbol")
  )
readr::write_csv(deg_ampk, file.path(OUT, "DEG_AMPK_pathway.csv"))
logf("INFO", "DEG_AMPK_pathway.csv: " %+% nrow(deg_ampk) %+% " genes da via AMPK")

# GSEA da via AMPK (fgsea). NOTA: testada isoladamente, padj (BH, 1 teste) = pval.
# O padj original (0.419) vem do GSEA com os 187 pathways KEGG (GSEA_KEGG_all_pathways.csv).
if (requireNamespace("fgsea", quietly = TRUE)) {
  ranked <- deg_full |>
    dplyr::filter(!is.na(log2FoldChange)) |>
    dplyr::arrange(desc(log2FoldChange))
  stats <- setNames(ranked$log2FoldChange, ranked$gene)

  gsea_ampk <- fgsea::fgseaSimple(pathways = list(AMPK_hsa04152 = ampk_symbols),
                                  stats = stats, nperm = 10000)
  gsea_out <- data.frame(
    pathway = "KEGG_AMPK_SIGNALING_PATHWAY_HSA04152",
    pval = gsea_ampk$pval,
    padj = gsea_ampk$padj,
    log2err = gsea_ampk$log2err,
    ES = gsea_ampk$ES,
    NES = gsea_ampk$NES,
    size = gsea_ampk$size,
    leadingEdge = vapply(gsea_ampk$leadingEdge, paste, collapse = "; ", FUN.VALUE = character(1))
  )
  readr::write_csv(gsea_out, file.path(OUT, "GSEA_AMPK_results.csv"))
  logf("INFO", sprintf("PASSO 1 | GSEA AMPK: NES = %.3f | padj = %.3f (NAO significativo)",
                       gsea_out$NES[1], gsea_out$padj[1]))
} else {
  logf("WARN", "fgsea indisponivel; GSEA_AMPK_results.csv nao gerado. Instale: BiocManager::install('fgsea')")
}

# GSVA da via AMPK (por amostra) — DCM vs CTRL
if (requireNamespace("GSVA", quietly = TRUE)) {
  gsva_res <- GSVA::gsva(expr = expr_log[, keep, drop = FALSE],
                         gset.idx.list = list(AMPK = ampk_symbols),
                         method = "gsva", kcdf = "Gaussian", verbose = FALSE)
  gsva_df <- data.frame(
    sample = colnames(gsva_res),
    AMPK_score = as.numeric(gsva_res["AMPK", ]),
    group = design_df$group[match(colnames(gsva_res), design_df$sample)]
  )
  readr::write_csv(gsva_df, file.path(OUT, "GSVA_AMPK_results.csv"))
  logf("INFO", sprintf("PASSO 1 | GSVA AMPK: %d amostras (DCM=%d, CTRL=%d)",
                       nrow(gsva_df), sum(gsva_df$group == "DCM"), sum(gsva_df$group == "CTRL")))
} else {
  logf("WARN", "GSVA indisponivel; GSVA_AMPK_results.csv nao gerado. Instale: BiocManager::install('GSVA')")
}

logf("INFO", "Pipeline upstream (PASSO 1 + PASSO 2) concluido.")
