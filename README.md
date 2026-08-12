# NQO1 na Cardiomiopatia Dilatada e o Repressor Farmacológico Doxorrubicina

**Análise integrativa: via AMPK → transcriptoma completo → NQO1 → Drug Gene Budger (DGB) → Doxorrubicina**

Repositório: `nqo1-doxorrubicina-cardiotoxicidade` · Auditoria de dados concluída (ver [`AUDIT_REPORT.md`](AUDIT_REPORT.md)) · Narrativa oficial: [`NARRATIVA_FINAL_OFICIAL.md`](NARRATIVA_FINAL_OFICIAL.md)

---

## 1. Introdução

Este estudo investiga, por abordagem **in silico**, a regulação do gene **NQO1** (NAD(P)H quinona desidrogenase 1) no contexto da **cardiomiopatia dilatada (DCM)** e identifica fármacos capazes de modular sua expressão. O fluxo lógico segue **três passos encadeados**:

### PASSO 1 — Via AMPK (KEGG: hsa04152) como ponto de partida
A investigação iniciou-se pela análise de expressão diferencial da via de sinalização **AMPK**
(`hsa04152`, KEGG) na DCM, a partir do dataset **GSE116250** (RNA-seq de miocárdio humano:
14 não-falha, 37 DCM, 13 ICM). A via AMPK foi avaliada por **limma**, **GSEA** e **GSVA**, e
**não apresentou enriquecimento estatisticamente significativo** (**NES = −1,19; padj = 0,419**).
Esta foi a **primeira parte de um estudo exploratório maior**: partiu-se dela, mas o resultado
negativo foi devidamente registrado.

### PASSO 2 — Do transcriptoma completo ao gene NQO1
Em seguida, foram analisados **todos os genes diferencialmente expressos do transcriptoma
completo** (GSE116250). Dentre eles, selecionou-se o gene com **significância estatística mais
robusta e relevância biológica**: o **NQO1**, significativamente **downregulado na DCM**
(**log2 Fold Change = −1,686; padj = 6,86 × 10⁻¹³**). A seleção baseou-se na função
citoprotetora/redox do NQO1 — **não** por pertencer à via AMPK (à qual o gene **não** pertence).

### PASSO 3 — DGB e a identificação da Doxorrubicina
Com o alvo NQO1 definido, consultou-se a plataforma **Drug Gene Budger (DGB)**, integrando as
bases **CMAP**, **LINCS L1000** e **CREEDS**. O fármaco com **resultado mais significativo e
evidência convergente em duas bases independentes** foi a **Doxorrubicina**
(**q = 2,21 × 10⁻²⁸**; Fold Change médio = −1,065; downregulador; CREEDS + LINCS L1000, n = 15).

### Nota sobre o ranking: Doxorrubicina × Trichostatin A
O **trichostatin A** aparece com `q = 0`, valor decorrente de **underflow numérico** (fonte única,
CMAP), capado em 10⁻⁵⁰. Portanto, a **Doxorrubicina** é o fármaco com o **menor q-value real**
do ranking (2,21 × 10⁻²⁸) e o único do topo com **replicação em múltiplas plataformas
independentes**, sendo a resposta correta do estudo.

> ⚠️ **Escopo:** a associação doxorrubicina↔NQO1 provém de bases de assinatura farmacológica
> (DGB), **não** de um ensaio de viabilidade próprio. A hipótese de cardiotoxicidade é
> **indireta/integrativa** e requer validação experimental.

---

## Documentação

| Documento | Conteúdo |
|---|---|
| [`NARRATIVA_FINAL_OFICIAL.md`](NARRATIVA_FINAL_OFICIAL.md) | Narrativa oficial do estudo — os 3 passos (AMPK → NQO1 → Doxorrubicina) |
| [`JUSTIFICATIVA_CIENTIFICA.md`](JUSTIFICATIVA_CIENTIFICA.md) | Racional científico robusto do desenho do estudo |
| [`IMPLICACOES_CLINICAS_ENFERMAGEM.md`](IMPLICACOES_CLINICAS_ENFERMAGEM.md) | Implicações clínicas e de enfermagem (com limitações honestas) |
| [`AUDIT_REPORT.md`](AUDIT_REPORT.md) | Auditoria de integridade de dados (pré-commit) |
| [`ANALISES_BIOINFORMATICA_PASSO_A_PASSO.md`](ANALISES_BIOINFORMATICA_PASSO_A_PASSO.md) | Guia passo a passo das 7 análises bioinformáticas (dados secundários) |
| [`TODO_FUTURE_ANALYSES.md`](TODO_FUTURE_ANALYSES.md) | Roadmap priorizado de análises futuras |
| [`config.yaml`](config.yaml) · [`DESCRIPTION`](DESCRIPTION) · [`requirements.txt`](requirements.txt) | Configuração e dependências |

---

## 2. Estrutura do Repositório

```
.
├── NARRATIVA_FINAL_OFICIAL.md           # Narrativa oficial (3 passos)
├── JUSTIFICATIVA_CIENTIFICA.md          # Racional científico do desenho
├── IMPLICACOES_CLINICAS_ENFERMAGEM.md   # Implicações clínicas e de enfermagem
├── ANALISES_BIOINFORMATICA_PASSO_A_PASSO.md  # Guia passo a passo (7 análises)
├── TODO_FUTURE_ANALYSES.md              # Roadmap de análises futuras
├── README.md                            # Este documento
├── AUDIT_REPORT.md                      # Auditoria de dados (pré-commit)
├── config.yaml                          # Thresholds + contexto centralizados
├── DESCRIPTION                          # Dependências R (Imports/Suggests)
├── requirements.txt                     # Dependências Python (auditoria/QC)
├── .gitignore
│
├── script.R                     # PIPELINE DOWNSTREAM (PASSO 3): DGB → ranking → figuras → PPI
├── scripts/
│   ├── ampk_de_analysis.R       # PIPELINE UPSTREAM (PASSO 1 + PASSO 2): DE + AMPK + NQO1 + GSEA/GSVA
│   └── audit_data.py            # Auditoria reprodutível (Python/pandas)
├── reports/QC_report.qmd        # Relatório HTML de QC (Quarto/RMarkdown)
│
├── data/raw/                    # DADOS BRUTOS (imutáveis)
│   ├── GSE116250_rpkm.txt.gz             # RNA-seq RPKM (Cufflinks)
│   ├── GSE116250_series_matrix.txt.gz    # Metadados GEO
│   ├── GPL16791.soft.gz                  # Anotação de plataforma
│   └── 9606.protein.*.v11.5.*.gz         # STRING v11.5 (PPI)
│
└── output/
    ├── tables/{nqo1,ampk,pubchem,ppi}/   # Tabelas derivadas
    ├── figures/{...}/                    # Figuras
    ├── rdata/  audit/  archive/          # Objetos R, logs, arquivos antigos
```

---

## 3. Fluxo Analítico Detalhado

### 3.1 Pré-processamento (script.R, downstream)
1. Leitura das 6 abas do `DGB_results_NQO1.xlsx` (CMAP/LINCS/CREEDS × up/down), com **remoção da coluna-índice residual `X1`**.
2. Conversão robusta de tipos (`decimal_to_numeric`), incluindo formato europeu de decimal.
3. Tratamento de `q = 0` (underflow): cap de `−log₁₀(q)` em 50 e registro do valor bruto como `<1e-300` (`q_value_raw`).
4. **Documentação de `Fold_Change`** como **escore sinalizado** (não razão): colunas `effect_score` e `log2FC_derived`.
5. Unificação case-insensitive dos nomes de fármacos (`Doxorubicin` = `doxorubicin`).

### 3.2 Análise Estatística
- **Upstream (PASSO 1/2, `ampk_de_analysis.R`):** limma com ajuste por **sexo + idade**
  (`~ 0 + group + sex + age`), contraste `DCM − CTRL`, FDR (Benjamini–Hochberg) α = 0,05;
  GSEA (fgsea) e GSVA da via AMPK.
- **Downstream (PASSO 3, `script.R`):** ranking por `score = −log₁₀(q) × |Fold Change|`,
  ordenado por q-value mínimo e escore.

### 3.3 Resultados-chave
| Etapa | Métrica | Valor |
|---|---|---|
| PASSO 1 | GSEA AMPK | NES = −1,19 · padj = 0,419 (não significativo) |
| PASSO 2 | NQO1 (DEG) | log2FC = −1,686 · padj = 6,86 × 10⁻¹³ (down) |
| PASSO 3 | Doxorrubicina (DGB) | q = 2,21 × 10⁻²⁸ · FC médio = −1,065 (down, 2 bases) |

---

## 4. Instruções de Reprodução

### 4.1 Dependências
```r
# R — pipeline (ver DESCRIPTION)
install.packages(c("dplyr","ggplot2","stringr","openxlsx","igraph","tidyr","data.table","readr","yaml"))
# Bioconductor — DE/GSEA/GSVA
if (!requireNamespace("BiocManager", quietly=TRUE)) install.packages("BiocManager")
BiocManager::install(c("limma","edgeR","GSVA","fgsea","GEOquery"))
```
```bash
# Python — auditoria/QC
pip install -r requirements.txt
```
> **Reprodutibilidade total:** rode `renv::init()` + `renv::snapshot()` para gerar o `renv.lock`
> (o `DESCRIPTION` já lista as dependências).

### 4.2 Ordem de execução
```bash
# 1) UPSTREAM — PASSO 1 (AMPK) + PASSO 2 (transcriptoma -> NQO1)
Rscript scripts/ampk_de_analysis.R

# 2) DOWNSTREAM — PASSO 3 (DGB -> ranking -> Doxorrubicina -> figuras -> PPI)
Rscript script.R

# 3) Auditoria reprodutível
python scripts/audit_data.py

# 4) Relatório HTML de QC
quarto render reports/QC_report.qmd     # ou rmarkdown::render("reports/QC_report.Rmd")
```

### 4.3 Saídas geradas
| Script | Saída |
|---|---|
| `ampk_de_analysis.R` | `DEG_full_table.csv`, `NQO1_results.csv`, `DEG_AMPK_pathway.csv`, `GSEA_AMPK_results.csv`, `GSVA_AMPK_results.csv` (em `output/tables/ampk/`) |
| `script.R` | `ranking_drugs_NQO1.csv`, `top20_drugs_NQO1.csv`, `doxorubicin_NQO1_summary.csv`, figuras volcano/dose-response, rede PPI |

### 4.4 Consulta DGB/PubChem (opcional)
Por padrão, são usados o `DGB_results_NQO1.xlsx` e os CSVs PubChem já salvos. Para refazer,
usar variáveis de ambiente (`RUN_PUBCHEM_QUERY`, `NCBI_API_KEY`) — **nunca** chaves no código.

---

## 5. Controle de Qualidade e Auditoria

- **QC visual:** `quarto render reports/QC_report.qmd` (boxplot antes/depois, PCA, densidade, % genes detectados).
- **Auditoria de dados:** `python scripts/audit_data.py` reproduz as checagens do `AUDIT_REPORT.md`.
- **Logs:** `script.R` e `ampk_de_analysis.R` registram timestamp + nível (`INFO`/`WARN`/`ERROR`)
  em `output/audit/logs/*.log` e `output/audit/audit_events_*.csv`.

---

## 6. Declaração de Transparência Algorítmica

> *"Este projeto utilizou a ferramenta de Inteligência Artificial generativa **deepseek-v4-pro** para auxiliar na automação de tarefas repetitivas de codificação (scripting), na estruturação e correção sintática da documentação, bem como na geração de rotinas de auditoria de dados. O uso desta ferramenta está em conformidade com a Portaria CNPq nº 2.664/2026, que regulamenta o uso de IA em pesquisas financiadas pela agência, garantindo a revisão crítica e validação final por um pesquisador titular (humano) responsável pelos resultados."*

---

## 7. Dados e Citação

- **GSE116250** (Sweet et al.): RNA-seq de miocárdio humano (DCM/ICM/não-falha) — GEO/NCBI.
- **STRING v11.5** (Szklarczyk et al.).
- **CMAP** (Lamb et al., 2006) · **LINCS L1000** (Subramanian et al., 2017) · **CREEDS**.
- **Drug Gene Budger (DGB)** (Wang et al., Bioinformatics, 2021).
