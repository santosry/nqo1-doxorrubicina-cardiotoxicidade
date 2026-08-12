#!/usr/bin/env Rscript
# ==============================================================================
# ampk_de_analysis.R — Reconstrucao do pipeline upstream (DE + AMPK + NQO1 + GSEA/GSVA)
# ==============================================================================
# Objetivo: reproduzir, a partir dos DADOS BRUTOS (data/raw/), as etapas que antes
# eram executadas em scripts externos nao versionados:
#
#   1. Leitura do GSE116250 (RNA-seq, RPKM) + metadados (pheno).
#   2. Expressao diferencial (limma, DCM vs. nao-falha [CTRL], ajustado por sexo+idade).
#   3. Identificacao do gene NQO1 (downregulado na DCM) e extracao da via AMPK (hsa04152).
#   4. GSEA (fgsea) e GSVA da via AMPK.
#   5. Gravacao das tabelas em output/tables/ampk/ (mesmos nomes ja usados).
#
# NOTA IMPORTANTE (correcoes de auditoria):
#   * GSE116250 e RNA-seq (nao microarray). O nome "preprocessed_microarray.RData"
#     e um termo historico incorreto; a matriz e RPKM (Cufflinks).
#   * NQO1 NAO pertence a via AMPK (hsa04152). Ele foi identificado na analise do
#     TRANSCRIPTOMA COMPLETO (DEG_full_table), nao no filtro da via AMPK.
#   * NQO1 esta DOWNREGULADO na DCM (log2FC ~ -1.69), nao upregulado.
#   * A via AMPK NAO apresentou enriquecimento significativo (GSEA padj ~ 0.419).
# ==============================================================================

suppressPackageStartupMessages({
  library(readr)
  library(dplyr)
  library(tidyr)
  library(limma)
  library(GSVA)
  library(fgsea)
})

# -----------------------------------------------------------------------------
# Logging simples (espelha o padrao do script.R)
# -----------------------------------------------------------------------------
log_dir <- file.path("output", "audit", "logs")
dir.create(log_dir, recursive = TRUE, showWarnings = FALSE)
logf <- function(level, msg) {
  line <- sprintf("[%s] %-5s | %s", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), level, msg)
  message(line)
  cat(line, "\n", file = file.path(log_dir, "ampk_de_analysis.log"), append = TRUE)
}

req <- function(pkg) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    stop("Pacote obrigatorio ausente: ", pkg, " (Bioconductor). Instale antes de rodar.")
  }
}

# Helper de concatenacao (evita repetir paste0 nos logs)
`%+%` <- function(a, b) paste0(a, b)

# -----------------------------------------------------------------------------
# 1. Leitura dos dados brutos
# -----------------------------------------------------------------------------
rpkm_path <- "data/raw/GSE116250_rpkm.txt.gz"
pheno_path <- "output/tables/ampk/pheno_data_raw.csv"
ampk_genes_path <- "output/tables/ampk/AMPK_pathway_genes_KEGG.csv"

logf("INFO", "Lendo matriz RPKM (RNA-seq): " %+% rpkm_path)
rpkm <- readr::read_tsv(rpkm_path, show_col_types = FALSE)

# Matriz de expressao: colunas 1 (Gene/Ensembl) e 2 (Common_name) + amostras
expr_mat <- as.matrix(rpkm[, -(1:2)])
rownames(expr_mat) <- rpkm$Common_name
storage.mode(expr_mat) <- "numeric"

logf("INFO", "Dimensoes da matriz: " %+% paste(dim(expr_mat), collapse = " x "))

# Metadados (derivados da series matrix via GEOquery na analise original)
pheno <- readr::read_csv(pheno_path, show_col_types = FALSE) |>
  dplyr::select(sample = title, disease = `disease:ch1`, sex = `Sex:ch1`, age = `age:ch1`) |>
  dplyr::mutate(
    # Padronizacao de rotulo do grupo controle: non-failing/NF -> CTRL
    group = dplyr::case_when(
      disease == "non-failing" ~ "CTRL",
      disease == "dilated cardiomyopathy" ~ "DCM",
      disease == "ischemic cardiomyopathy" ~ "ICM",
      TRUE ~ NA_character_
    )
  ) |>
  dplyr::filter(sample %in% colnames(expr_mat))

logf("INFO", "Amostras com metadado: " %+% nrow(pheno) %+% " (de " %+% ncol(expr_mat) %+% " na matriz)")

# -----------------------------------------------------------------------------
# 2. Expressao diferencial — DCM vs CTRL (ajustado por sexo + idade)
#    ICM e excluida deste contraste (foco DCM), como no GSVA original.
# -----------------------------------------------------------------------------
logf("INFO", "Normalizacao: log2(RPKM + 1) e filtro de genes pouco expressos")
expr_log <- log2(expr_mat + 1)

# Filtra genes expressos em pelo menos 25% das amostras usadas
keep_samples <- pheno$sample[pheno$group %in% c("DCM", "CTRL")]
keep_genes <- rowSums(expr_mat[, keep_samples, drop = FALSE] > 1) >= (0.25 * length(keep_samples))
expr_filt <- expr_log[keep_genes, keep_samples, drop = FALSE]

design_df <- pheno |>
  dplyr::filter(group %in% c("DCM", "CTRL")) |>
  dplyr::mutate(group = factor(group, levels = c("CTRL", "DCM"))) |>
  dplyr::arrange(match(sample, colnames(expr_filt)))

design <- model.matrix(~ 0 + group + sex + age, data = design_df)
colnames(design) <- gsub("^group", "", colnames(design))

logf("INFO", "Design (covariaveis): " %+% paste(colnames(design), collapse = " + "))
logf("INFO", "Contraste: DCM - CTRL")

fit <- limma::lmFit(expr_filt, design)
contr <- limma::makeContrasts(DCM - CTRL, levels = design)
fit2 <- limma::contrasts.fit(fit, contr)
fit2 <- limma::eBayes(fit2, trend = TRUE)

deg_full <- limma::topTable(fit2, number = Inf, sort.by = "none") |>
  tibble::rownames_to_column("gene") |>
  dplyr::mutate(
    neg_log10_padj = -log10(pmax(adj.P.Val, 1e-300)),
    direction = ifelse(logFC > 0, "up", "down"),
    significant = adj.P.Val < 0.05
  ) |>
  dplyr::select(gene, log2FoldChange = logFC, AveExpr, t, pvalue = P.Value, padj = adj.P.Val,
                B, direction, significant, neg_log10_padj)

dir.create("output/tables/ampk", recursive = TRUE, showWarnings = FALSE)
readr::write_csv(deg_full, "output/tables/ampk/DEG_full_table.csv")
logf("INFO", "DEG_full_table.csv gravado: " %+% nrow(deg_full) %+% " genes")

# -----------------------------------------------------------------------------
# 3. NQO1 (alvo farmacogenomico)
# -----------------------------------------------------------------------------
nqo1 <- deg_full[deg_full$gene == "NQO1", ]
if (nrow(nqo1) == 0) {
  logf("WARN", "NQO1 nao encontrado no DEG.")
} else {
  logf("INFO", sprintf("NQO1: log2FC = %.3f | padj = %.3e | direcao = %s",
                       nqo1$log2FoldChange, nqo1$padj, nqo1$direction))
}

# -----------------------------------------------------------------------------
# 4. Via AMPK (hsa04152, KEGG) — filtro dos genes da via + GSEA + GSVA
# -----------------------------------------------------------------------------
logf("INFO", "Carregando genes da via AMPK (hsa04152) de " %+% ampk_genes_path)
ampk_genes <- readr::read_csv(ampk_genes_path, show_col_types = FALSE)
ampk_symbols <- unique(ampk_genes$symbol)

deg_ampk <- deg_full |>
  dplyr::filter(gene %in% ampk_symbols) |>
  dplyr::left_join(
    ampk_genes |> dplyr::select(symbol, description) |> dplyr::distinct(symbol, .keep_all = TRUE),
    by = c("gene" = "symbol")
  )

readr::write_csv(deg_ampk, "output/tables/ampk/DEG_AMPK_pathway.csv")
logf("INFO", "DEG_AMPK_pathway.csv: " %+% nrow(deg_ampk) %+% " genes da via AMPK")
logf("INFO", "NQO1 pertence a via AMPK (hsa04152)? " %+%
       ifelse("NQO1" %in% ampk_symbols, "SIM", "NAO (correto: foi identificado no transcriptoma completo)"))

# GSEA da via AMPK (fgsea) — usando o ranking de todos os genes
# NOTA: aqui a via AMPK e testada isoladamente, logo padj (BH, 1 teste) = pval.
# O padj original reportado no manuscrito (0.419) vem do GSEA com os 187 conjuntos
# KEGG (ver GSEA_KEGG_all_pathways.csv), aplicando BH sobre todos os pathways.
ranked <- deg_full |>
  dplyr::filter(!is.na(log2FoldChange)) |>
  dplyr::arrange(desc(log2FoldChange))
stats <- setNames(ranked$log2FoldChange, ranked$gene)

gsea_ampk <- fgsea::fgseaSimple(pathways = list(AMPK_hsa04152 = ampk_symbols), stats = stats, nperm = 10000)
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
readr::write_csv(gsea_out, "output/tables/ampk/GSEA_AMPK_results.csv")
logf("INFO", sprintf("GSEA AMPK: NES = %.3f | padj = %.3f", gsea_out$NES, gsea_out$padj))

# GSVA da via AMPK (por amostra) — DCM vs CTRL
gsva_res <- GSVA::gsva(expr = expr_log[, keep_samples, drop = FALSE],
                       gset.idx.list = list(AMPK = ampk_symbols),
                       method = "gsva", kcdf = "Gaussian", verbose = FALSE)
gsva_df <- data.frame(
  sample = colnames(gsva_res),
  AMPK_score = as.numeric(gsva_res["AMPK", ]),
  group = design_df$group[match(colnames(gsva_res), design_df$sample)]
)
readr::write_csv(gsva_df, "output/tables/ampk/GSVA_AMPK_scores_per_sample.csv")
logf("INFO", sprintf("GSVA AMPK gravado: %d amostras (DCM=%d, CTRL=%d)",
                     nrow(gsva_df),
                     sum(gsva_df$group == "DCM"),
                     sum(gsva_df$group == "CTRL")))

logf("INFO", "Pipeline upstream reconstruido concluido.")
