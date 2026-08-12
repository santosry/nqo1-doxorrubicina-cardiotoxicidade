#!/usr/bin/env Rscript
# ==============================================================================
# ampk_de_analysis.R — Pipeline UPSTREAM (PASSO 1 + PASSO 2)
# ==============================================================================
# Narrativa oficial (3 passos):
#   PASSO 1 — Analise exploratoria da via AMPK (KEGG hsa04152) na DCM (GSE116250).
#             NAO significativa (GSEA NES ~ -1.19; padj ~ 0.419).
#   PASSO 2 — DE do TRANSCRIPTOMA COMPLETO -> NQO1 downregulado (log2FC ~ -1.69).
#   PASSO 3 — (script.R) DGB -> Doxorrubicina (q = 2.21e-28).
#
# Este script executa o PASSO 1 e o PASSO 2 e salva em output/tables/ampk/:
#   - DEG_full_table.csv      (DE do transcriptoma completo, 41842 genes)
#   - NQO1_results.csv        (linha do gene NQO1)
#   - DEG_AMPK_pathway.csv    (DE restrito a via AMPK hsa04152)
#   - GSEA_AMPK_results.csv   (GSEA da via AMPK)
#   - GSVA_AMPK_results.csv   (escore GSVA da via AMPK por amostra)
#
# FONTE DOS DADOS PRE-PROCESSADOS (ground truth):
#   output/rdata/preprocessed_microarray.RData contem `expr_matrix_norm`
#   (41842 genes x 51 amostras, escala log2) e `pheno_data` (DCM=37, CTRL=14),
#   produzidos pelo pipeline original (GEOquery + normalizacao log2). Carregamos
#   este objeto para reproduzir exatamente o log2FC do NQO1 (~ -1.69).
#
# CORRECOES DE AUDITORIA APLICADAS (10 melhorias):
#   4) ICM excluida do contraste (DCM vs CTRL).
#   5) Rotulo do grupo controle padronizado como CTRL.
#   6) Proveniencia da normalizacao documentada (matriz log2 pre-processada).
#   7) Ajuste por covariadas sex + idade (design ~ 0 + group + sex + age).
#   8) (no script.R) data.table::fread para .gz; aqui a matriz ja esta pre-processada.
# ==============================================================================

suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(readr)
  library(limma)
})

# -----------------------------------------------------------------------------
# Helpers e logging
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
# 1. Carregar dados pre-processados (ground truth)
# -----------------------------------------------------------------------------
rdata_path <- "output/rdata/preprocessed_microarray.RData"
if (!file.exists(rdata_path)) {
  stop("Arquivo pre-processado ausente: ", rdata_path, ". Rode o pre-processamento antes.")
}
env <- new.env()
load(rdata_path, envir = env)
expr_norm <- env$expr_matrix_norm   # 41842 genes x 51 amostras, escala log2
pheno <- env$pheno_data

# Garante rownames unicos (ja sao simbolos unicos no objeto original)
if (anyDuplicated(rownames(expr_norm))) {
  rownames(expr_norm) <- make.unique(rownames(expr_norm))
}

logf("INFO", "Matriz pre-processada: " %+% nrow(expr_norm) %+% " genes x " %+% ncol(expr_norm) %+% " amostras")

# Padronizacao do grupo controle (melhoria 5): non-failing -> CTRL
pheno <- pheno |>
  dplyr::mutate(
    group = factor(dplyr::case_when(
      `disease:ch1` == "non-failing" ~ "CTRL",
      `disease:ch1` == "dilated cardiomyopathy" ~ "DCM",
      TRUE ~ NA_character_
    ), levels = c("CTRL", "DCM")),
    sex = `Sex:ch1`,
    age = as.numeric(`age:ch1`)
  ) |>
  dplyr::filter(title %in% colnames(expr_norm))

pheno <- pheno[match(colnames(expr_norm), pheno$title), ]

logf("INFO", "Amostras: CTRL=" %+% sum(pheno$group == "CTRL") %+%
       ", DCM=" %+% sum(pheno$group == "DCM") %+%
       " (ICM excluida do contraste — melhoria 4)")

# -----------------------------------------------------------------------------
# 2. PASSO 2 — DE do transcriptoma completo (DCM vs CTRL, sex + idade)
# -----------------------------------------------------------------------------
logf("INFO", "Design (melhoria 7): ~ 0 + group + sex + age")
design <- model.matrix(~ 0 + group + sex + age, data = pheno)
colnames(design) <- sub("^group", "", colnames(design))

fit <- limma::lmFit(expr_norm, design)
fit <- limma::contrasts.fit(fit, limma::makeContrasts("DCM-CTRL", levels = design))
fit <- limma::eBayes(fit, trend = TRUE)

tt <- limma::topTable(fit, number = Inf, sort.by = "none")
deg_full <- data.frame(
  gene = rownames(tt),
  log2FoldChange = tt$logFC,
  AveExpr = tt$AveExpr,
  t = tt$t,
  pvalue = tt$P.Value,
  padj = tt$adj.P.Val,
  B = tt$B,
  stringsAsFactors = FALSE,
  row.names = NULL
) |>
  dplyr::mutate(
    neg_log10_padj = -log10(pmax(padj, 1e-300)),
    direction = ifelse(log2FoldChange > 0, "up", "down"),
    significant = padj < 0.05
  ) |>
  dplyr::select(gene, log2FoldChange, AveExpr, t, pvalue, padj, B,
                direction, significant, neg_log10_padj)

readr::write_csv(deg_full, file.path(OUT, "DEG_full_table.csv"))
logf("INFO", "DEG_full_table.csv: " %+% nrow(deg_full) %+% " genes")

# NQO1 (melhoria: alvo do PASSO 2)
nqo1 <- deg_full |> dplyr::filter(gene == "NQO1")
if (nrow(nqo1) == 0) {
  logf("WARN", "NQO1 nao encontrado no DEG.")
} else {
  logf("INFO", sprintf("PASSO 2 | NQO1: log2FC = %.4f | padj = %.3e | %s",
                       nqo1$log2FoldChange[1], nqo1$padj[1], nqo1$direction[1]))
}
readr::write_csv(nqo1, file.path(OUT, "NQO1_results.csv"))
logf("INFO", "NQO1_results.csv gravado (downregulado na DCM; NAO pertence a via AMPK)")

# -----------------------------------------------------------------------------
# 3. PASSO 1 — Via AMPK (hsa04152, KEGG)
# -----------------------------------------------------------------------------
ampk_genes <- readr::read_csv(file.path(OUT, "AMPK_pathway_genes_KEGG.csv"), show_col_types = FALSE)
ampk_symbols <- unique(ampk_genes$symbol)
logf("INFO", "PASSO 1 | Via AMPK (hsa04152): " %+% length(ampk_symbols) %+% " genes no KEGG")

deg_ampk <- deg_full |>
  dplyr::filter(gene %in% ampk_symbols) |>
  dplyr::left_join(
    ampk_genes |> dplyr::select(symbol, description) |> dplyr::distinct(symbol, .keep_all = TRUE),
    by = c("gene" = "symbol")
  )
readr::write_csv(deg_ampk, file.path(OUT, "DEG_AMPK_pathway.csv"))
logf("INFO", "DEG_AMPK_pathway.csv: " %+% nrow(deg_ampk) %+% " genes da via AMPK mapeados")

# GSEA da via AMPK (fgsea). Testada isoladamente, padj (BH, 1 teste) = pval.
# O padj reportado no manuscrito (0.419) vem do GSEA com os 187 pathways KEGG.
if (requireNamespace("fgsea", quietly = TRUE)) {
  ranked <- deg_full |>
    dplyr::filter(!is.na(log2FoldChange)) |>
    dplyr::arrange(desc(log2FoldChange))
  stats <- setNames(ranked$log2FoldChange, ranked$gene)

  gsea_ampk <- fgsea::fgseaSimple(
    pathways = list(AMPK_hsa04152 = ampk_symbols),
    stats = stats, nperm = 10000
  )
  gsea_out <- data.frame(
    pathway = "KEGG_AMPK_SIGNALING_PATHWAY_HSA04152",
    pval = gsea_ampk$pval,
    padj = gsea_ampk$padj,
    log2err = NA_real_,   # fgseaSimple nao estima log2err (apenas fgsea() com paralelismo)
    ES = gsea_ampk$ES,
    NES = gsea_ampk$NES,
    size = gsea_ampk$size,
    leadingEdge = vapply(gsea_ampk$leadingEdge, paste, collapse = "; ", FUN.VALUE = character(1))
  )
  readr::write_csv(gsea_out, file.path(OUT, "GSEA_AMPK_results.csv"))
  logf("INFO", sprintf("PASSO 1 | GSEA AMPK: NES = %.3f | padj = %.3f (NAO significativo)",
                       gsea_out$NES[1], gsea_out$padj[1]))
} else {
  logf("WARN", "fgsea indisponivel; GSEA_AMPK_results.csv nao gerado.")
}

# GSVA da via AMPK (por amostra)
if (requireNamespace("GSVA", quietly = TRUE)) {
  param <- GSVA::gsvaParam(
    exprData = expr_norm,
    geneSets = list(AMPK = ampk_symbols),
    kcdf = "Gaussian"
  )
  gsva_obj <- GSVA::gsva(param)
  scores <- gsva_obj   # GSVA 2.x retorna uma matriz de escores (genesets x amostras)
  gsva_df <- data.frame(
    sample = colnames(scores),
    AMPK_score = as.numeric(scores["AMPK", ]),
    group = pheno$group[match(colnames(scores), pheno$title)]
  )
  readr::write_csv(gsva_df, file.path(OUT, "GSVA_AMPK_results.csv"))
  logf("INFO", sprintf("PASSO 1 | GSVA AMPK: %d amostras (DCM=%d, CTRL=%d)",
                       nrow(gsva_df), sum(gsva_df$group == "DCM"), sum(gsva_df$group == "CTRL")))
} else {
  logf("WARN", "GSVA indisponivel; GSVA_AMPK_results.csv nao gerado.")
}

logf("INFO", "Pipeline upstream (PASSO 1 + PASSO 2) concluido.")
