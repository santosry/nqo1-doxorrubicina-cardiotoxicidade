# Ranqueamento Farmacogenômico in silico de Fármacos Moduladores da Expressão do Gene NQO1

**Integração das bases CMAP, LINCS L1000 e CREEDS via Drug Gene Budger (DGB)**

Repositório: `nqo1-doxorrubicina-cardiotoxicidade` · Status: auditoria pré-commit concluída (ver [`AUDIT_REPORT.md`](AUDIT_REPORT.md))

---

## 1. Resumo Executivo

A **NQO1** (NAD(P)H quinona desidrogenase 1) é uma enzima citoprotetora que catalisa a redução bieletrônica de quinonas, prevenindo a formação de espécies reativas de oxigênio. Em uma análise prévia de expressão do miocárdio humano, a NQO1 mostrou-se **significativamente reprimida na cardiomiopatia dilatada (DCM)**, o que motivou sua investigação como alvo farmacogenômico.

O estudo integra três bases de assinaturas transcriptômicas de perturbação por fármacos (**CMAP**, **LINCS L1000**, **CREEDS**) por meio da plataforma **Drug Gene Budger (DGB)**, identificando e ranqueando compostos capazes de modular a transcrição de NQO1. Entre os moduladores, a **doxorrubicina** destacou-se como o **repressor transcricional mais significativo de NQO1** (q = 2,21 × 10⁻²⁸; Fold Change médio = −1,065), com evidência convergente em duas plataformas independentes (CREEDS e LINCS L1000).

---

## 2. Fluxo Analítico Definitivo (passo a passo)

> Esta seção registra **exatamente** o que foi feito, com as correções de auditoria aplicadas.

### Etapa 1 — Expressão diferencial do transcriptoma completo (GSE116250)
Foi realizada a análise de expressão diferencial do transcriptoma humano (GSE116250, **RNA-seq**, 64 amostras: 14 não-falha, 37 DCM, 13 ICM) com o pacote **limma**.

- O gene **NQO1** foi identificado como o **2º gene mais significativamente desregulado**, com **expressão REDUZIDA na DCM** (`log2FoldChange = −1,686`; `padj = 6,86 × 10⁻¹³`).
- ⚠️ **Correção importante:** NQO1 está **DOWNREGULADO** (e não upregulado). A direção negativa indica menor expressão no miocárdio patológico.

### Etapa 2 — Caracterização da via AMPK (hsa04152, KEGG)
Paralelamente, a via de sinalização **AMPK** (hsa04152) foi avaliada por **GSEA** e **GSVA**.

- A via AMPK **NÃO apresentou enriquecimento estatisticamente significativo** (`NES = −1,19`; `p = 0,125`; `padj = 0,419`).
- ⚠️ **Correção importante:** a NQO1 **NÃO pertence à via canônica AMPK (hsa04152)**. Ela foi selecionada a partir da análise do **transcriptoma completo** (Etapa 1), por sua relevância funcional na citoproteção redox — e não por ter sido "filtrada da via AMPK".

### Etapa 3 — Seleção do alvo NQO1
A NQO1 foi escolhida como alvo da triagem farmacogenômica com base em: (i) forte significância no DEG (padj = 6,86 × 10⁻¹³); e (ii) papel citoprotetor/redox, com relação direta à cardiotoxicidade por antraciclinas.

### Etapa 4 — Varredura farmacogenômica no Drug Gene Budger (DGB)
O gene NQO1 foi submetido à varredura de assinaturas de perturbação por xenobióticos no **DGB**, integrando **CMAP**, **LINCS L1000** e **CREEDS**.

- **3.523 registros válidos** · **1.266 fármacos únicos** (pós-unificação case-insensitive).
- Concentrações 0,1–1000 µM · tempos 6 h/24 h · linhagens MCF7, PC3, VCAP, entre outras.

### Etapa 5 — Ranqueamento e identificação da doxorrubicina
Um escore composto (`score = −log₁₀(q-value) × |Fold Change|`) integrou significância e magnitude. O ranking foi ordenado por q-value mínimo e escore.

- **Doxorrubicina** = repressor de NQO1 mais significativo com q-value real (não-zero) mais baixo (q = 2,21 × 10⁻²⁸), convergente em **CREEDS + LINCS L1000** (n = 15 registros).
- O topo do ranking inclui ainda trichostatin A (up; CMAP), parthenolide (up; LINCS) e panobinostat (down; LINCS) — ver `output/tables/nqo1/top20_drugs_NQO1.csv`.

> **Leitura correta do "mais significativo":** trichostatin A tem q = 0 (artefato de *underflow*, capado em 10⁻⁵⁰, fonte única CMAP). A doxorrubicina tem o **menor q-value real** (2,21 × 10⁻²⁸) e é o único do top 5 com evidência convergente em **duas bases independentes**.

---

## 3. Estrutura do Repositório

```
.
├── config.yaml                  # Thresholds + contexto do projeto centralizados
├── DESCRIPTION                  # Manifesto de dependências R (Imports/Suggests)
├── requirements.txt             # Dependências Python (auditoria/QC)
├── AUDIT_REPORT.md              # Auditoria de dados (pré-commit)
├── README.md                    # Este documento
├── .gitignore
├── script.R                     # Pipeline downstream (DGB → ranking → figuras → PPI)
├── scripts/
│   ├── ampk_de_analysis.R       # RECONSTRUÍDO: DE + AMPK + NQO1 + GSEA/GSVA (upstream)
│   └── audit_data.py            # Auditoria reprodutível (Python/pandas)
├── reports/QC_report.qmd        # Relatório HTML de QC (Quarto/RMarkdown)
├── data/raw/                    # DADOS BRUTOS (imutáveis)
│   ├── GSE116250_rpkm.txt.gz             # RNA-seq RPKM (Cufflinks)
│   ├── GSE116250_series_matrix.txt.gz    # Metadados GEO
│   ├── GPL16791.soft.gz                  # Anotação de plataforma
│   └── 9606.protein.*.v11.5.*.gz         # STRING v11.5 (PPI)
└── output/
    ├── tables/{nqo1,ampk,pubchem,ppi}/
    ├── figures/{nqo1_drug_regulation,ampk,pubchem,ppi}/
    ├── rdata/  audit/  archive/
```

| Pasta | Função |
|---|---|
| `data/raw/` | Dados imutáveis baixados de GEO/STRING. |
| `scripts/` | Scripts reprodutíveis (upstream DE + downstream DGB + auditoria). |
| `output/tables/` | Tabelas derivadas (DEG, ranking, GSEA, GSVA, PPI). |
| `output/figures/` | Figuras para publicação. |
| `output/audit/` | Logs (`.csv`, `.log`) e relatórios de auditoria. |
| `output/archive/` | Artefatos antigos preservados para rastreio. |

---

## 4. Fluxo Analítico Detalhado

### 4.1 Pré-processamento
1. **Leitura das 6 abas DGB** (`cmap/l1000/creeds` × `up/down`), com **remoção da coluna-índice residual `X1`** (melhoria aplicada).
2. **Conversão robusta de tipos** (`decimal_to_numeric`), incluindo formato europeu de decimal.
3. **Tratamento de `q = 0`**: cap de `−log₁₀(q)` em 50 e registro do valor bruto como `<1e-300` (`q_value_raw`).
4. **Imputação de missing**: não aplicada (dados analíticos têm 0% de NA); método configurável (`imputation.method: knn`).
5. **Outliers (|Z| > 3)**: sinalizados, não removidos; ranking usa estatística robusta.
6. **Unificação case-insensitive** dos nomes de fármacos.

### 4.2 Análise Estatística
- **DE (upstream):** limma com ajuste por **sexo + idade** (`~ 0 + group + sex + age`), contraste `DCM − CTRL`, FDR (Benjamini–Hochberg) α = 0,05.
- **GSEA/GSVA:** via AMPK (hsa04152) — **não significativa** (padj = 0,419).
- **Ranking DGB:** agregação por fármaco (`mean_fc`, `median_fc`, `min_q`, `score = −log₁₀(q) × |FC|`).

### 4.3 Análise de Sensibilidade
- Sexo (51♂ / 13♀) e idade (20–66) são **covariáveis** no modelo DE (melhoria aplicada).
- ICM é **excluída do contraste DCM vs. CTRL** (documentado no `ampk_de_analysis.R`).

### 4.4 Machine Learning
Não aplicável no estado atual (sem alvo rotulado de viabilidade). Plano (Random Forest/Regressão Logística + SHAP) configurável em `config.yaml`.

---

## 5. Instruções de Reprodução

### 5.1 Dependências
```r
# R (pipeline): ver DESCRIPTION. Instalação mínima:
install.packages(c("dplyr","ggplot2","stringr","openxlsx","igraph","tidyr","data.table","readr","yaml"))
# Bioconductor (DE/GSEA/GSVA):
if (!requireNamespace("BiocManager", quietly=TRUE)) install.packages("BiocManager")
BiocManager::install(c("limma","edgeR","GSVA","fgsea","GEOquery"))
```
```bash
# Python (auditoria/QC):
pip install -r requirements.txt
```

> **Reprodutibilidade total (renv):** para gerar o `renv.lock` com as versões exatas, rode
> `renv::init()` seguido de `renv::snapshot()` na raiz do projeto. (O `DESCRIPTION` já lista
> as dependências; o `renv.lock` deve ser gerado no ambiente R local.)

### 5.2 Ordem de execução
```bash
# 1) UPSTREAM — expressão diferencial + via AMPK + NQO1 + GSEA/GSVA
Rscript scripts/ampk_de_analysis.R

# 2) DOWNSTREAM — DGB: leitura do XLSX, ranking, figuras, PPI
Rscript script.R

# 3) Auditoria reprodutível
python scripts/audit_data.py

# 4) Relatório HTML de QC
quarto render reports/QC_report.qmd     # ou rmarkdown::render("reports/QC_report.Rmd")
```

### 5.3 Consulta DGB/PubChem (opcional)
Por padrão, o script usa `output/tables/nqo1/DGB_results_NQO1.xlsx` e os CSVs PubChem já salvos.
Para refazer a consulta, usar variáveis de ambiente (`RUN_PUBCHEM_QUERY`, `NCBI_API_KEY`) — **nunca** chaves no código.

---

## 6. Verificação de Scripts e Resultados (o que existe × o que foi reconstruído)

| Etapa | Script | Resultado | Status |
|---|---|---|---|
| Download GSE116250 | — (GEOquery) | `data/raw/GSE116250_*` | ✅ dado presente · script ausente |
| DE (limma, DCM vs CTRL) | `scripts/ampk_de_analysis.R` | `DEG_full_table.csv` | 🔁 **reconstruído** |
| Identificação da NQO1 | `scripts/ampk_de_analysis.R` | linha NQO1 no DEG | 🔁 **reconstruído** |
| GSEA/GSVA via AMPK (hsa04152) | `scripts/ampk_de_analysis.R` | `GSEA_AMPK_results.csv`, `GSVA_AMPK_scores_per_sample.csv` | 🔁 **reconstruído** |
| Consulta DGB (geração do XLSX) | ferramenta web DGB | `DGB_results_NQO1.xlsx` | ✅ dado presente · consulta externa |
| Leitura/processamento do XLSX | `script.R` → `read_nqo1_xlsx()` | `DGB_results_NQO1_from_xlsx_clean.csv` | ✅ presente |
| Ranking / top20 / doxo | `script.R` → `rank_nqo1_drugs()` | `ranking_drugs_NQO1.csv`, `top20_*`, `doxorubicin_NQO1_summary.csv` | ✅ presente |
| Figuras (volcano, dose-response) | `script.R` | `volcano_NQO1_drugs.png`, `dose_response_NQO1_drugs.png` | ✅ presente |
| PPI (STRING + igraph) | `script.R` → `plot_ppi_network()` | `PPI_AMPK_STRING_{nodes,edges}.csv`, `PPI_AMPK_network.png` | ✅ presente |
| Auditoria de dados | `scripts/audit_data.py` | `AUDIT_REPORT.md` | ✅ presente |

**Conclusão da verificação:** os **resultados** de todas as etapas estão presentes. Os **scripts upstream** (DE, GSEA/GSVA, consulta DGB) **não estavam versionados** — a etapa de DE/GSEA/GSVA foi **reconstruída** em `scripts/ampk_de_analysis.R`. A consulta DGB em si é uma ferramenta web externa (não é scriptável localmente), mas seu resultado (`DGB_results_NQO1.xlsx`) está preservado.

---

## 7. Melhorias implementadas (auditoria de código)

1. ✅ Remoção da coluna-índice residual `X1` ao ler as abas do XLSX.
2. ✅ Documentação de `Fold_Change` como **escore sinalizado** + colunas `effect_score` e `log2FC_derived`.
3. ✅ Registro de `q_value == 0` como `<1e-300` (`q_value_raw`).
4. ✅ Exclusão de ICM documentada no `ampk_de_analysis.R`.
5. ✅ Padronização do grupo controle como `CTRL` (mapa `non-failing = NF = CTRL` em `config.yaml`).
6. ✅ Proveniência da normalização documentada (`log2(RPKM+1)` + limma; nota sobre `voom` para contagens).
7. ✅ Ajuste por covariáveis `sex + age` no modelo DE.
8. ✅ Performance: `data.table::fread` nos `.gz` do STRING.
9. ✅ `DESCRIPTION` (dependências R) + comando `renv::snapshot()` documentado.
10. ✅ Ensaio de viabilidade (MTT/SRB) sinalizado como ausente; threshold pronto em `config.yaml`.

---

## 8. Controle de Qualidade Visual (rápido)

```bash
quarto render reports/QC_report.qmd
```
Gráficos: boxplot antes/depois da normalização · PCA (batch effect) · densidade por amostra · % de genes detectados.

---

## 9. Declaração de Transparência Algorítmica

> *"Este projeto utilizou a ferramenta de Inteligência Artificial generativa **deepseek-v4-pro** para auxiliar na automação de tarefas repetitivas de codificação (scripting), na estruturação e correção sintática da documentação, bem como na geração de rotinas de auditoria de dados. O uso desta ferramenta está em conformidade com a Portaria CNPq nº 2.664/2026, que regulamenta o uso de IA em pesquisas financiadas pela agência, garantindo a revisão crítica e validação final por um pesquisador titular (humano) responsável pelos resultados."*

---

## 10. Dados e Citação

- **GSE116250** (Sweet et al.): RNA-seq de miocárdio humano (DCM/ICM/não-falha) — GEO/NCBI.
- **STRING v11.5** (Szklarczyk et al.).
- **CMAP** (Lamb et al., 2006) · **LINCS L1000** (Subramanian et al., 2017) · **CREEDS**.
- **Drug Gene Budger (DGB)** (Wang et al., Bioinformatics, 2021).
