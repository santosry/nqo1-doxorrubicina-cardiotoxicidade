#!/usr/bin/env Rscript

# ==============================================================================
# CARDIO - script unico de analise
# ==============================================================================
# Projeto de analise de expressao diferencial cardiovascular, com foco atual nos
# farmacos que modulam a resposta do gene NQO1.
#
# Estrutura organizada:
#   data/raw/
#     Dados brutos baixados de GEO/STRING. Nao editar manualmente.
#
#   output/tables/ampk/
#     Tabelas da analise de expressao diferencial/AMPK.
#
#   output/tables/pubchem/
#     Resultados de associacao gene-composto via PubChem.
#
#   output/tables/nqo1/
#     Fonte Excel, tabelas limpas, ranking e top 20 de farmacos NQO1.
#
#   output/figures/ampk/
#     Figuras da analise AMPK que ainda sao uteis. PCA foi removido.
#
#   output/figures/pubchem/
#     Graficos derivados do PubChem.
#
#   output/figures/nqo1_drug_regulation/
#     Graficos finais de regulacao de NQO1 por farmacos.
#
#   output/figures/ppi/
#     Rede PPI da via AMPK gerada a partir dos arquivos STRING.
#
#   output/audit/
#     Logs estruturados e relatorio de auditoria.
#
#   output/archive/
#     Resultados antigos ou fora do foco atual, preservados para rastreio.
#
# Como rodar no PowerShell, a partir da raiz do projeto:
#   & 'C:\Program Files\R\R-4.6.0\bin\Rscript.exe' script.R
#
# Consulta PubChem:
#   Por padrao, o script usa os resultados PubChem ja salvos em output/tables/pubchem.
#   Para consultar o PubChem/NCBI novamente, rode com:
#     $env:RUN_PUBCHEM_QUERY='true'
#     $env:NCBI_API_KEY='sua_chave_ncbi_opcional'
#   A chave NCBI deve ficar no ambiente, nunca escrita dentro do script.
#
# Resultados principais:
#   output/tables/nqo1/doxorubicin_NQO1_summary.csv
#   output/figures/nqo1_drug_regulation/volcano_NQO1_drugs.png
#   output/figures/nqo1_drug_regulation/dose_response_NQO1_drugs.png
#   output/figures/pubchem/pubchem_top_genes_compound_counts.png
#   output/figures/ppi/PPI_AMPK_network.png
#   output/audit/AUDIT_NQO1_PUBCHEM_PPI.md

options(stringsAsFactors = FALSE)
set.seed(42)

get_project_root <- function() {
  args <- commandArgs(trailingOnly = FALSE)
  file_arg <- grep("^--file=", args, value = TRUE)
  if (length(file_arg) > 0) {
    return(normalizePath(dirname(sub("^--file=", "", file_arg[1])), winslash = "/", mustWork = TRUE))
  }
  normalizePath(getwd(), winslash = "/", mustWork = TRUE)
}

ROOT <- get_project_root()
setwd(ROOT)

DIRS <- list(
  data_raw = file.path(ROOT, "data", "raw"),
  tables_ampk = file.path(ROOT, "output", "tables", "ampk"),
  tables_pubchem = file.path(ROOT, "output", "tables", "pubchem"),
  tables_nqo1 = file.path(ROOT, "output", "tables", "nqo1"),
  tables_ppi = file.path(ROOT, "output", "tables", "ppi"),
  figures_ampk = file.path(ROOT, "output", "figures", "ampk"),
  figures_pubchem = file.path(ROOT, "output", "figures", "pubchem"),
  figures_nqo1 = file.path(ROOT, "output", "figures", "nqo1_drug_regulation"),
  figures_ppi = file.path(ROOT, "output", "figures", "ppi"),
  audit = file.path(ROOT, "output", "audit"),
  logs = file.path(ROOT, "output", "audit", "logs")
)
invisible(lapply(DIRS, dir.create, recursive = TRUE, showWarnings = FALSE))

AUDIT_EVENTS <- data.frame(
  time = character(),
  level = character(),
  step = character(),
  message = character()
)

audit_log <- function(level, step, message) {
  stamp <- format(Sys.time(), "%Y-%m-%d %H:%M:%S")
  line <- sprintf("[%s] %-5s | %-18s | %s", stamp, level, step, message)
  message(line)
  # Espelha a linha em um arquivo .log diario (timestamp + nivel + mensagem)
  log_file <- file.path(DIRS$logs, sprintf("run_%s.log", format(Sys.Date(), "%Y%m%d")))
  cat(line, "\n", file = log_file, append = TRUE)
  AUDIT_EVENTS <<- rbind(
    AUDIT_EVENTS,
    data.frame(time = stamp, level = level, step = step, message = message)
  )
}

require_pkg <- function(pkg) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    stop("Pacote obrigatorio ausente: ", pkg, ". Instale antes de rodar este script.")
  }
  suppressPackageStartupMessages(library(pkg, character.only = TRUE))
}

require_pkg("dplyr")
require_pkg("ggplot2")
require_pkg("stringr")
require_pkg("openxlsx")
require_pkg("igraph")
require_pkg("tidyr")
require_pkg("data.table")

find_first_existing <- function(paths, required = TRUE, label = "arquivo") {
  existing <- paths[file.exists(paths)]
  if (length(existing) > 0) {
    normalizePath(existing[1], winslash = "/", mustWork = TRUE)
  } else if (required) {
    stop("Nao encontrei ", label, ". Caminhos testados: ", paste(paths, collapse = " | "))
  } else {
    NA_character_
  }
}

save_plot <- function(plot, filename, width = 9, height = 6) {
  ggplot2::ggsave(
    filename = filename,
    plot = plot,
    width = width,
    height = height,
    dpi = 300,
    bg = "white"
  )
  audit_log("INFO", "plot", paste("Figura salva:", normalizePath(filename, winslash = "/", mustWork = FALSE)))
}

read_csv_flexible <- function(path) {
  audit_log("INFO", "read", paste("Lendo", path))
  df <- utils::read.csv(path, check.names = FALSE, stringsAsFactors = FALSE)
  names(df) <- trimws(names(df))
  df
}

standardize_names <- function(x) {
  x |>
    stringr::str_replace_all("\\s+", "_") |>
    stringr::str_replace_all("-", "_") |>
    stringr::str_replace_all("\\.", "_")
}

cleanup_unneeded_outputs <- function() {
  unneeded <- c(
    file.path(DIRS$figures_nqo1, "boxplot_top_drugs_NQO1.png"),
    file.path(DIRS$figures_nqo1, "top20_drugs_NQO1.png"),
    file.path(DIRS$figures_ampk, "PCA_final_publication.png"),
    file.path(DIRS$figures_ampk, "PCA_microarray.png"),
    file.path(DIRS$figures_ampk, "PCA_scree_plot.png"),
    file.path(DIRS$tables_nqo1, "DGB_results_NQO1.csv"),
    file.path(DIRS$tables_nqo1, "DGB_results_NQO1_clean.csv")
  )
  for (path in unneeded[file.exists(unneeded)]) {
    file.remove(path)
    audit_log("INFO", "cleanup", paste("Removido por nao ser mais necessario:", path))
  }
}

decimal_to_numeric <- function(x) {
  x <- as.character(x)
  has_comma <- stringr::str_detect(x, ",")
  x <- ifelse(
    has_comma,
    stringr::str_replace(stringr::str_replace_all(x, "\\.", ""), ",", "."),
    x
  )
  suppressWarnings(as.numeric(x))
}

q_to_neglog10 <- function(q_value, cap = 50) {
  y <- -log10(pmax(q_value, 10^-cap, na.rm = TRUE))
  pmin(y, cap)
}

extract_dose_um <- function(dose) {
  dose <- as.character(dose)
  dose <- stringr::str_replace_all(dose, ",", ".")
  value <- suppressWarnings(as.numeric(stringr::str_extract(dose, "\\d+(?:\\.\\d+)?")))
  unit <- stringr::str_to_lower(stringr::str_extract(dose, "nm|µm|μm|um|mm"))
  dplyr::case_when(
    is.na(value) ~ NA_real_,
    unit == "nm" ~ value / 1000,
    unit == "mm" ~ value * 1000,
    TRUE ~ value
  )
}

read_nqo1_xlsx <- function(path) {
  sheets <- openxlsx::getSheetNames(path)
  audit_log("INFO", "xlsx", paste("Abas encontradas:", paste(sheets, collapse = ", ")))

  out <- lapply(sheets, function(sheet) {
    df <- openxlsx::read.xlsx(path, sheet = sheet)
    # Remove coluna-indice residual (1a coluna sem nome, lida como X1 pela openxlsx).
    if (ncol(df) > 0 && grepl("^X\\d+$", names(df)[1])) {
      df <- df[, -1, drop = FALSE]
    }
    names(df) <- standardize_names(names(df))

    if (!"Drug_Name" %in% names(df)) stop("Aba sem coluna Drug Name: ", sheet)
    if (!"p_value" %in% names(df)) stop("Aba sem coluna p-value: ", sheet)
    if (!"q_value" %in% names(df)) stop("Aba sem coluna q-value: ", sheet)
    if (!"Fold_Change" %in% names(df)) stop("Aba sem coluna Fold Change: ", sheet)

    source <- dplyr::case_when(
      stringr::str_detect(sheet, "cmap") ~ "CMAP",
      stringr::str_detect(sheet, "l1000") ~ "LINCS L1000",
      stringr::str_detect(sheet, "creeds") ~ "CREEDS",
      TRUE ~ sheet
    )
    expected_direction <- ifelse(stringr::str_detect(sheet, "_up$"), "up", "down")

    df |>
      dplyr::mutate(
        dplyr::across(dplyr::any_of(c("p_value", "q_value", "Fold_Change")), as.character),
        source_database = source,
        source_sheet = sheet,
        expected_direction = expected_direction
      )
  })

  combined <- dplyr::bind_rows(out)

  combined |>
    dplyr::filter(!is.na(Drug_Name), Drug_Name != "") |>
    dplyr::mutate(
      Drug_Name = stringr::str_squish(as.character(Drug_Name)),
      Drug_Name_Norm = stringr::str_to_lower(stringr::str_squish(as.character(Drug_Name))),
      p_value = decimal_to_numeric(p_value),
      q_value = decimal_to_numeric(q_value),
      Fold_Change = decimal_to_numeric(Fold_Change),
      # Fold_Change do DGB e um ESCORE SINALIZADO (sinal = direcao), nao uma razao
      # de fold-change classica. Mantemos o nome por compatibilidade e adicionamos
      # semantica explicita (effect_score) + uma derivada em escala log2.
      effect_score = Fold_Change,
      log2FC_derived = sign(Fold_Change) * log2(abs(Fold_Change) + 1),
      q_value_for_plot = pmax(q_value, 1e-50),
      neg_log10_q_capped = q_to_neglog10(q_value, cap = 50),
      q_was_zero = !is.na(q_value) & q_value == 0,
      q_value_raw = ifelse(!is.na(q_value) & q_value == 0, "<1e-300", as.character(q_value)),
      regulation = dplyr::case_when(
        Fold_Change > 0 ~ "Upregula NQO1",
        Fold_Change < 0 ~ "Downregula NQO1",
        TRUE ~ "Sem direcao"
      ),
      score = neg_log10_q_capped * abs(Fold_Change)
    ) |>
    dplyr::filter(!is.na(p_value), !is.na(q_value), !is.na(Fold_Change))
}

rank_nqo1_drugs <- function(df) {
  # Normaliza nomes de farmacos para case-insensitive (ex: Doxorubicin = doxorubicin)
  # Mantem o nome original mais frequente como display name
  df |>
    dplyr::group_by(Drug_Name_Norm) |>
    dplyr::summarise(
      Drug_Name = names(sort(table(Drug_Name), decreasing = TRUE))[1],
      mean_fc = mean(Fold_Change, na.rm = TRUE),
      median_fc = median(Fold_Change, na.rm = TRUE),
      max_abs_fc = max(abs(Fold_Change), na.rm = TRUE),
      min_p_value = min(p_value, na.rm = TRUE),
      min_q_value = min(q_value, na.rm = TRUE),
      n_tests = dplyr::n(),
      n_sources = dplyr::n_distinct(source_database),
      sources = paste(sort(unique(source_database)), collapse = "; "),
      score = max(score, na.rm = TRUE),
      direction = dplyr::case_when(
        median_fc > 0 ~ "Upregula NQO1",
        median_fc < 0 ~ "Downregula NQO1",
        TRUE ~ "Sem direcao"
      ),
      .groups = "drop"
    ) |>
    dplyr::arrange(min_q_value, dplyr::desc(score), dplyr::desc(max_abs_fc))
}

plot_volcano <- function(df) {
  top_labels <- df |>
    dplyr::arrange(q_value, dplyr::desc(abs(Fold_Change))) |>
    dplyr::slice_head(n = 12)

  ggplot2::ggplot(df, ggplot2::aes(x = Fold_Change, y = neg_log10_q_capped)) +
    ggplot2::geom_point(ggplot2::aes(color = regulation, shape = q_was_zero), alpha = 0.74, size = 2.1) +
    ggplot2::geom_vline(xintercept = 0, linetype = "dashed", color = "grey45") +
    ggplot2::geom_text(
      data = top_labels,
      ggplot2::aes(label = Drug_Name),
      size = 2.8,
      check_overlap = TRUE,
      vjust = -0.65,
      color = "grey20"
    ) +
    ggplot2::scale_color_manual(values = c(
      "Upregula NQO1" = "#1B9E77",
      "Downregula NQO1" = "#D95F02",
      "Sem direcao" = "grey60"
    )) +
    ggplot2::scale_shape_manual(values = c("FALSE" = 16, "TRUE" = 17)) +
    ggplot2::labs(
      title = "Volcano de farmacos moduladores de NQO1",
      subtitle = "Vistoria: q-values iguais a zero foram capados para evitar escala distorcida",
      x = "Fold change",
      y = "-log10(q-value), capado em 50",
      color = "Regulacao",
      shape = "q-value = 0"
    ) +
    ggplot2::theme_minimal(base_size = 12) +
    ggplot2::theme(legend.position = "bottom")
}

plot_dose_response <- function(df, ranking, max_drugs = 9) {
  if (!"Dose" %in% names(df)) df$Dose <- "Dose nao informada"
  if (!"Time" %in% names(df)) df$Time <- NA_character_

  dose_df_all <- df |>
    dplyr::mutate(
      Dose = as.character(Dose),
      Time = as.character(Time),
      Dose = dplyr::if_else(is.na(Dose) | Dose == "", "Dose nao informada", Dose),
      dose_um = extract_dose_um(Dose),
      time_label = dplyr::if_else(is.na(Time) | Time == "", "Tempo nao informado", Time),
      drug_key = Drug_Name_Norm,
      Drug_Display = dplyr::if_else(
        stringr::str_detect(drug_key, "doxo|doxorubic|doxorubicin"),
        "Doxorubicin",
        Drug_Name
      )
    ) |>
    dplyr::filter(!is.na(dose_um), dose_um > 0)

  if (nrow(dose_df_all) == 0) {
    stop("Nao ha doses numericas na planilha para montar dose-response.")
  }

  dose_drugs_main <- dose_df_all |>
    dplyr::group_by(Drug_Display) |>
    dplyr::summarise(
      n_doses = dplyr::n_distinct(dose_um),
      n_rows = dplyr::n(),
      min_q_value = min(q_value, na.rm = TRUE),
      score = max(score, na.rm = TRUE),
      .groups = "drop"
    ) |>
    dplyr::filter(n_doses >= 2) |>
    dplyr::arrange(min_q_value, dplyr::desc(score)) |>
    dplyr::slice_head(n = max_drugs)

  doxo_drug <- dose_df_all |>
    dplyr::filter(Drug_Display == "Doxorubicin") |>
    dplyr::group_by(Drug_Display) |>
    dplyr::summarise(
      n_doses = dplyr::n_distinct(dose_um),
      n_rows = dplyr::n(),
      min_q_value = min(q_value, na.rm = TRUE),
      score = max(score, na.rm = TRUE),
      .groups = "drop"
    )

  dose_drugs <- dplyr::bind_rows(dose_drugs_main, doxo_drug) |>
    dplyr::distinct(Drug_Display, .keep_all = TRUE) |>
    dplyr::arrange(min_q_value, dplyr::desc(score))

  if (nrow(dose_drugs) == 0) {
    dose_drugs <- dose_df_all |>
      dplyr::group_by(Drug_Display) |>
      dplyr::summarise(
        n_rows = dplyr::n(),
        min_q_value = min(q_value, na.rm = TRUE),
        score = max(score, na.rm = TRUE),
        .groups = "drop"
      ) |>
      dplyr::arrange(min_q_value, dplyr::desc(score)) |>
      dplyr::slice_head(n = max_drugs)
  }

  dose_df <- dose_df_all |>
    dplyr::filter(Drug_Display %in% dose_drugs$Drug_Display) |>
    dplyr::mutate(
      Drug_Display = factor(Drug_Display, levels = rev(dose_drugs$Drug_Display)),
      series = paste(source_database, Cell_Line, time_label, sep = " | ")
    )

  ggplot2::ggplot(dose_df, ggplot2::aes(x = dose_um, y = Fold_Change)) +
    ggplot2::geom_hline(yintercept = 0, linetype = "dashed", color = "grey55", linewidth = 0.35) +
    ggplot2::geom_line(
      ggplot2::aes(group = series, color = source_database),
      alpha = 0.35,
      linewidth = 0.45
    ) +
    ggplot2::geom_point(
      ggplot2::aes(size = neg_log10_q_capped, fill = regulation),
      shape = 21,
      color = "grey25",
      stroke = 0.25,
      alpha = 0.88
    ) +
    ggplot2::facet_wrap(~ Drug_Display, scales = "free_y", ncol = 3) +
    ggplot2::scale_x_log10() +
    ggplot2::scale_color_manual(values = c("CMAP" = "#4C78A8", "LINCS L1000" = "#7F7F7F", "CREEDS" = "#B279A2")) +
    ggplot2::scale_fill_manual(values = c(
      "Upregula NQO1" = "#1B9E77",
      "Downregula NQO1" = "#D95F02",
      "Sem direcao" = "grey70"
    )) +
    ggplot2::scale_size_continuous(range = c(1.5, 5.5)) +
    ggplot2::labs(
      title = "Dose-response dos farmacos que modulam NQO1",
      subtitle = "Somente registros com dose numerica; Doxorubicin incluida a partir dos registros LINCS com 10 µM",
      x = "Dose (µM, escala log10)",
      y = "Fold change de NQO1",
      size = "-log10(q) capado",
      color = "Base",
      fill = "Regulacao"
    ) +
    ggplot2::theme_minimal(base_size = 11) +
    ggplot2::theme(
      legend.position = "bottom",
      panel.grid.minor = ggplot2::element_blank(),
      strip.text = ggplot2::element_text(face = "bold")
    )
}

load_pubchem <- function() {
  path <- find_first_existing(c(
    file.path(DIRS$tables_pubchem, "pubchem_results_all.csv"),
    file.path(ROOT, "output", "pubchem_results_all.csv"),
    file.path(ROOT, "pubchem_results_all.csv")
  ), required = FALSE, label = "resultado PubChem")

  if (is.na(path)) {
    audit_log("WARN", "pubchem", "Resultado PubChem nao encontrado; grafico PubChem nao sera recriado.")
    return(data.frame())
  }
  read_csv_flexible(path)
}

plot_pubchem_counts <- function(pubchem_df) {
  if (nrow(pubchem_df) == 0 || !"Gene" %in% names(pubchem_df)) return(NULL)

  pubchem_df |>
    dplyr::count(Gene, Direction, name = "n_compounds") |>
    dplyr::arrange(dplyr::desc(n_compounds)) |>
    dplyr::slice_head(n = 20) |>
    ggplot2::ggplot(ggplot2::aes(x = n_compounds, y = reorder(Gene, n_compounds), fill = Direction)) +
    ggplot2::geom_col(width = 0.72) +
    ggplot2::scale_fill_manual(values = c("up" = "#1B9E77", "down" = "#D95F02"), na.value = "grey60") +
    ggplot2::labs(
      title = "Top genes AMPK por compostos associados no PubChem",
      x = "Numero de compostos",
      y = NULL,
      fill = "Regulacao do gene"
    ) +
    ggplot2::theme_minimal(base_size = 12) +
    ggplot2::theme(legend.position = "bottom")
}

plot_ppi_network <- function() {
  deg_path <- find_first_existing(c(
    file.path(DIRS$tables_ampk, "DEG_AMPK_pathway.csv"),
    file.path(ROOT, "DEG_AMPK_pathway.csv")
  ), required = FALSE, label = "DEG_AMPK_pathway.csv")
  links_path <- find_first_existing(file.path(DIRS$data_raw, "9606.protein.links.v11.5.min700.txt.gz"),
                                    required = FALSE, label = "STRING links")
  info_path <- find_first_existing(file.path(DIRS$data_raw, "9606.protein.info.v11.5.txt.gz"),
                                   required = FALSE, label = "STRING protein info")

  if (is.na(deg_path) || is.na(links_path) || is.na(info_path)) {
    audit_log("WARN", "ppi", "Arquivos DEG/STRING ausentes; rede PPI nao sera gerada.")
    return(NULL)
  }

  deg <- read_csv_flexible(deg_path) |>
    dplyr::mutate(gene = as.character(gene))
  genes <- unique(deg$gene)

  audit_log("INFO", "ppi", "Lendo STRING protein info.")
  info <- as.data.frame(data.table::fread(info_path, quote = "", stringsAsFactors = FALSE))
  names(info) <- standardize_names(names(info))
  if ("X_string_protein_id" %in% names(info)) {
    names(info)[names(info) == "X_string_protein_id"] <- "string_protein_id"
  }
  info_sub <- info |>
    dplyr::filter(preferred_name %in% genes) |>
    dplyr::select(string_protein_id, preferred_name)

  if (nrow(info_sub) == 0) {
    audit_log("WARN", "ppi", "Nenhum gene AMPK foi mapeado para STRING.")
    return(NULL)
  }

  audit_log("INFO", "ppi", "Lendo STRING links.")
  links <- as.data.frame(data.table::fread(links_path, header = TRUE, stringsAsFactors = FALSE))
  names(links) <- standardize_names(names(links))
  ids <- unique(info_sub$string_protein_id)
  links_sub <- links |>
    dplyr::filter(protein1 %in% ids, protein2 %in% ids) |>
    dplyr::left_join(info_sub, by = c("protein1" = "string_protein_id")) |>
    dplyr::rename(from = preferred_name) |>
    dplyr::left_join(info_sub, by = c("protein2" = "string_protein_id")) |>
    dplyr::rename(to = preferred_name) |>
    dplyr::filter(!is.na(from), !is.na(to), from != to) |>
    dplyr::arrange(dplyr::desc(combined_score))

  if (nrow(links_sub) == 0) {
    audit_log("WARN", "ppi", "Nenhuma interacao STRING entre genes AMPK mapeados.")
    return(NULL)
  }

  utils::write.csv(links_sub, file.path(DIRS$tables_ppi, "PPI_AMPK_STRING_edges.csv"), row.names = FALSE)

  g <- igraph::graph_from_data_frame(
    d = links_sub |> dplyr::select(from, to, combined_score),
    directed = FALSE
  )
  node_df <- data.frame(gene = igraph::V(g)$name) |>
    dplyr::left_join(deg |> dplyr::select(gene, log2FoldChange, padj, direction), by = "gene")
  igraph::V(g)$degree <- igraph::degree(g)
  igraph::V(g)$log2FoldChange <- node_df$log2FoldChange
  igraph::V(g)$direction <- node_df$direction

  layout <- igraph::layout_with_fr(g, weights = igraph::E(g)$combined_score)
  coords <- data.frame(gene = igraph::V(g)$name, x = layout[, 1], y = layout[, 2]) |>
    dplyr::left_join(node_df, by = "gene") |>
    dplyr::mutate(degree = igraph::degree(g)[gene])

  edge_df <- igraph::as_data_frame(g, what = "edges") |>
    dplyr::left_join(coords |> dplyr::select(gene, x, y), by = c("from" = "gene")) |>
    dplyr::rename(x_from = x, y_from = y) |>
    dplyr::left_join(coords |> dplyr::select(gene, x, y), by = c("to" = "gene")) |>
    dplyr::rename(x_to = x, y_to = y)

  utils::write.csv(coords, file.path(DIRS$tables_ppi, "PPI_AMPK_STRING_nodes.csv"), row.names = FALSE)

  p <- ggplot2::ggplot() +
    ggplot2::geom_segment(
      data = edge_df,
      ggplot2::aes(x = x_from, y = y_from, xend = x_to, yend = y_to, linewidth = combined_score),
      color = "grey72",
      alpha = 0.55
    ) +
    ggplot2::geom_point(
      data = coords,
      ggplot2::aes(x = x, y = y, size = degree, fill = log2FoldChange),
      shape = 21,
      color = "grey20",
      stroke = 0.35
    ) +
    ggplot2::geom_text(
      data = coords |> dplyr::arrange(dplyr::desc(degree)) |> dplyr::slice_head(n = 35),
      ggplot2::aes(x = x, y = y, label = gene),
      size = 2.8,
      check_overlap = TRUE,
      vjust = -0.8
    ) +
    ggplot2::scale_fill_gradient2(low = "#D95F02", mid = "white", high = "#1B9E77") +
    ggplot2::scale_linewidth(range = c(0.15, 1.4), guide = "none") +
    ggplot2::labs(
      title = "Rede PPI da via AMPK",
      subtitle = "STRING Homo sapiens, interacoes min700; cor = log2FC em DCM vs controle",
      x = NULL,
      y = NULL,
      fill = "log2FC",
      size = "Grau"
    ) +
    ggplot2::theme_void(base_size = 11) +
    ggplot2::theme(legend.position = "bottom")

  save_plot(p, file.path(DIRS$figures_ppi, "PPI_AMPK_network.png"), width = 10, height = 8)
  audit_log("INFO", "ppi", paste("Rede PPI:", igraph::vcount(g), "nos e", igraph::ecount(g), "arestas."))
  invisible(list(graph = g, nodes = coords, edges = links_sub))
}

main <- function() {
  audit_log("INFO", "start", paste("Raiz do projeto:", ROOT))
  cleanup_unneeded_outputs()

  xlsx_path <- find_first_existing(c(
    file.path(DIRS$tables_nqo1, "DGB_results_NQO1.xlsx"),
    file.path(ROOT, "DGB_results_NQO1.xlsx")
  ), label = "DGB_results_NQO1.xlsx")

  dgb <- read_nqo1_xlsx(xlsx_path)
  audit_log("INFO", "nqo1", paste("Linhas validas na planilha NQO1:", nrow(dgb)))
  audit_log("INFO", "nqo1", paste("Farmacos unicos (pos-unificacao case-insensitive):", dplyr::n_distinct(dgb$Drug_Name_Norm)))

  ranking <- rank_nqo1_drugs(dgb)
  top20 <- ranking |> dplyr::slice_head(n = 20)

  utils::write.csv(dgb, file.path(DIRS$tables_nqo1, "DGB_results_NQO1_from_xlsx_clean.csv"), row.names = FALSE)
  utils::write.csv(ranking, file.path(DIRS$tables_nqo1, "ranking_drugs_NQO1.csv"), row.names = FALSE)
  utils::write.csv(top20, file.path(DIRS$tables_nqo1, "top20_drugs_NQO1.csv"), row.names = FALSE)

  # Tabela 1: distribuicao por plataforma e direcao
  tab1_platform <- ranking |>
    dplyr::mutate(
      platform_group = dplyr::case_when(
        sources == "LINCS L1000" ~ "LINCS L1000",
        sources == "CREEDS" ~ "CREEDS",
        sources == "CMAP" ~ "CMAP",
        TRUE ~ "Multi-plataforma"
      )
    ) |>
    dplyr::count(platform_group, direction, name = "n_farmacos") |>
    tidyr::pivot_wider(
      names_from = direction,
      values_from = n_farmacos,
      values_fill = 0
    )
  utils::write.csv(tab1_platform, file.path(DIRS$tables_nqo1, "table1_platform_direction.csv"), row.names = FALSE)
  audit_log("INFO", "table1", paste("Tabela 1 gerada: plataformas x direcao"))

  doxo <- ranking |>
    dplyr::filter(stringr::str_detect(Drug_Name_Norm, "doxo|doxorubic|doxorubicin"))
  utils::write.csv(doxo, file.path(DIRS$tables_nqo1, "doxorubicin_NQO1_summary.csv"), row.names = FALSE)
  if (nrow(doxo) > 0) {
    audit_log("INFO", "nqo1", paste("Doxorubicin encontrada (pos-unificacao); melhor q-value:", signif(min(doxo$min_q_value), 4)))
    audit_log("INFO", "nqo1", paste("Direcao da doxorubicina:", doxo$direction[1], "| mean_fc:", signif(doxo$mean_fc[1], 4)))
  } else {
    audit_log("WARN", "nqo1", "Doxorubicin nao foi encontrada no ranking.")
  }

  save_plot(plot_volcano(dgb), file.path(DIRS$figures_nqo1, "volcano_NQO1_drugs.png"), 10, 7.5)
  save_plot(plot_dose_response(dgb, ranking),
            file.path(DIRS$figures_nqo1, "dose_response_NQO1_drugs.png"), 11.5, 8.5)

  pubchem <- load_pubchem()
  if (nrow(pubchem) > 0) {
    utils::write.csv(pubchem, file.path(DIRS$tables_pubchem, "pubchem_results_all.csv"), row.names = FALSE)
    pubchem_named <- pubchem |> dplyr::filter(!is.na(Name), Name != "")
    utils::write.csv(pubchem_named, file.path(DIRS$tables_pubchem, "pubchem_results_with_names.csv"), row.names = FALSE)
    p_pub <- plot_pubchem_counts(pubchem)
    if (!is.null(p_pub)) {
      save_plot(p_pub, file.path(DIRS$figures_pubchem, "pubchem_top_genes_compound_counts.png"), 9, 6.5)
    }
  }

  plot_ppi_network()

  audit_csv <- file.path(DIRS$audit, "audit_events_nqo1_pubchem_ppi.csv")
  utils::write.csv(AUDIT_EVENTS, audit_csv, row.names = FALSE)

  # Stats for report
  n_down <- sum(ranking$direction == "Downregula NQO1")
  n_up <- sum(ranking$direction == "Upregula NQO1")
  n_multi <- sum(stringr::str_detect(ranking$sources, ";"))
  n_lincs <- sum(stringr::str_detect(ranking$sources, "LINCS"))
  n_creeds <- sum(stringr::str_detect(ranking$sources, "CREEDS"))
  n_cmap <- sum(stringr::str_detect(ranking$sources, "CMAP"))

  report <- c(
    "# Auditoria NQO1/PubChem/PPI",
    "",
    paste("- Data:", format(Sys.time(), "%Y-%m-%d %H:%M:%S")),
    paste("- Raiz:", ROOT),
    paste("- Fonte unica dos farmacos:", xlsx_path),
    paste("- Abas processadas:", paste(openxlsx::getSheetNames(xlsx_path), collapse = ", ")),
    paste("- Linhas NQO1 validas:", nrow(dgb)),
    paste("- Farmacos unicos (pos-unificacao case-insensitive):", nrow(ranking)),
    paste("- Downreguladores:", n_down, "(", round(100*n_down/nrow(ranking), 1), "%)"),
    paste("- Upreguladores:", n_up, "(", round(100*n_up/nrow(ranking), 1), "%)"),
    paste("- LINCS L1000:", n_lincs, "(", round(100*n_lincs/nrow(ranking), 1), "%)"),
    paste("- CREEDS:", n_creeds, "(", round(100*n_creeds/nrow(ranking), 1), "%)"),
    paste("- CMAP:", n_cmap, "(", round(100*n_cmap/nrow(ranking), 1), "%)"),
    paste("- Multi-plataforma:", n_multi, "(", round(100*n_multi/nrow(ranking), 1), "%)"),
    paste("- Top 1:", top20$Drug_Name[1]),
    paste("- Doxorubicin no ranking:", ifelse(nrow(doxo) > 0, "sim", "nao")),
    if (nrow(doxo) > 0) paste("- Doxorubicin direcao:", doxo$direction[1], "| mean_fc:", signif(doxo$mean_fc[1], 4), "| min_q:", signif(doxo$min_q_value[1], 4)) else "- Doxorubicin: nao encontrada",
    "- Volcano revisado: q-values iguais a zero capados em -log10(q)=50.",
    "- Boxplot top drugs removido conforme solicitado.",
    "- PCA removido conforme solicitado.",
    paste("- Tabela top20 preservada como tabela:", file.path(DIRS$tables_nqo1, "top20_drugs_NQO1.csv")),
    "- Grafico top20_drugs_NQO1 removido conforme solicitado.",
    "- Unificacao case-insensitive aplicada no ranking (ex: Doxorubicin + doxorubicin).",
    paste("- Figuras NQO1:", DIRS$figures_nqo1),
    paste("- Figura PPI:", file.path(DIRS$figures_ppi, "PPI_AMPK_network.png")),
    paste("- Log estruturado:", audit_csv),
    "",
    "## Top 20",
    paste(capture.output(print(top20[, c("Drug_Name", "min_q_value", "mean_fc", "score", "direction", "sources")], row.names = FALSE)), collapse = "\n")
  )
  writeLines(report, file.path(DIRS$audit, "AUDIT_NQO1_PUBCHEM_PPI.md"), useBytes = TRUE)
  audit_log("INFO", "done", "Analise concluida.")
}

main()
