# NQO1 como Modulador da Cardiotoxicidade por Doxorrubicina

**Repositório:** `nqo1-doxorrubicina-cardiotoxicidade`
**Status:** dados íntegros · auditoria pré-commit concluída (ver [`AUDIT_REPORT.md`](AUDIT_REPORT.md))

---

## 1. Resumo Executivo

A **NQO1** (NAD(P)H quinona desidrogenase 1) é uma flavoproteína com papel ambíguo na toxicidade de quinonas: em doses fisiológicas atua como **antioxidante/desintoxicante**, mas em contextos de estresse redox pode agir como **pró-fármaco redutor (bioativação)**. A **doxorrubicina**, antraciclina de primeira linha, tem sua **cardiotoxicidade** classicamente associada à geração de espécies reativas de oxigênio e à formação de metabólitos semiquinona — um processo em que a NQO1 participa diretamente da ciclagem redox.

Este repositório investiga, de forma **exploratória e integrativa**, se a modulação farmacológica e a expressão de NQO1 se correlacionam com o contexto cardiotóxico. O projeto combina três camadas de evidência:

1. **Expressão diferencial de NQO1** em miocárdio humano (GEO `GSE116250`, RNA-seq, n = 64) — a NQO1 está **significativamente reduzida** na cardiomiopatia dilatada (DCM): `log2FC = −1.69`, `padj = 6.86e-13`.
2. **Assinaturas farmacológicas** (CMAP, LINCS L1000, CREEDS) — a **doxorrubicina downregula NQO1** (`mean_fc = −1.07`, `q = 2.21e-28`, 15 registros).
3. **Rede de interação proteína–proteína (PPI)** da via AMPK (STRING v11.5) como contexto de estresse metabólico cardíaco.

> ⚠️ **Transparência de escopo (leia antes de prosseguir):** o dataset de expressão primário (`GSE116250`) é de **insuficiência cardíaca humana (DCM/ICM vs. não-falha)**, **não** um modelo experimental tratado com doxorrubicina. **Não há rótulo `Doxo` × `Controle`** no dado de expressão. A ligação doxorrubicina ↔ NQO1 provém de bancos de assinatura farmacológica. Portanto, a hipótese de cardiotoxicidade é sustentada de forma **indireta/integrativa** e deve ser lida como tal. Detalhes em [`AUDIT_REPORT.md`](AUDIT_REPORT.md) §2.1.

---

## 2. Estrutura do Repositório

```
.
├── config.yaml                  # Thresholds centralizados (viabilidade, FC, Z, FDR, ML)
├── requirements.txt             # Dependências Python (auditoria/QC)
├── AUDIT_REPORT.md              # Relatório de auditoria de dados (pré-commit)
├── README.md                    # Este documento
├── .gitignore                   # Lixo de SO/IDE/cache/segredos
├── script.R                     # Pipeline analítico principal (R)
│
├── data/
│   └── raw/                     # DADOS BRUTOS — nunca editar manualmente
│       ├── GSE116250_rpkm.txt.gz             # Matriz de expressão RPKM (RNA-seq)
│       ├── GSE116250_series_matrix.txt.gz    # Metadados GEO (série)
│       ├── GPL16791.soft.gz                  # Anotação da plataforma
│       └── 9606.protein.{info,links,aliases}.v11.5.*.gz  # STRING v11.5 (PPI)
│
├── output/
│   ├── tables/
│   │   ├── nqo1/                # Tabelas NQO1 × fármacos (ranking, top20, doxo)
│   │   ├── ampk/                # DEG, GSEA, GSVA e fenotipagem da via AMPK
│   │   ├── pubchem/             # Associação gene–composto (PubChem)
│   │   └── ppi/                 # Nós e arestas da rede PPI AMPK
│   ├── figures/
│   │   ├── nqo1_drug_regulation/  # Volcano + dose-response NQO1
│   │   ├── ampk/                  # Volcano, heatmaps, GSEA/GSVA
│   │   ├── pubchem/               # Lollipop/top genes por compostos
│   │   └── ppi/                   # Rede PPI AMPK
│   ├── rdata/                   # `preprocessed_microarray.RData`
│   ├── audit/                   # Logs e relatórios de auditoria
│   └── archive/                 # Resultados antigos, preservados p/ rastreio
│
└── (diretórios canônicos sugeridos para expandir o fluxo — ver §4)
    ├── data/processed/          # Derivados limpos (log2/TMM, imputados)
    ├── scripts/                 # Scripts auxiliares (auditoria, QC)
    ├── results/                 # Objetos intermediários/tabelas finais
    └── reports/                 # Relatórios HTML/PDF (Quarto/RMarkdown)
```

**Para que serve cada pasta:**

| Pasta | Função |
|---|---|
| `data/raw/` | Dados **imutáveis** baixados de GEO/STRING. Qualquer transformação gera um novo arquivo em `data/processed/` ou `output/`, nunca sobrescreve aqui. |
| `data/processed/` *(sugerido)* | Dados limpos/normalizados (log2/TMM, imputação, outliers tratados). |
| `scripts/` | Scripts auxiliares reprodutíveis (auditoria, QC, ML). |
| `output/tables/` | Tabelas derivadas (DEG, ranking, GSEA, GSVA, PPI). |
| `output/figures/` | Figuras prontas para publicação. |
| `output/rdata/` | Objetos R serializados (`.RData`). |
| `output/audit/` | Logs estruturados (`.csv`, `.log`) e relatórios de auditoria. |
| `output/archive/` | Artefatos obsoletos preservados para rastreabilidade (não apagar). |
| `results/` / `reports/` | Saídas finais e relatórios renderizados (HTML/PDF). |

---

## 3. Fluxo Analítico Detalhado (Passo a Passo)

O fluxo é **reproduzível por script** (`script.R`), sem nenhuma etapa manual em planilha.

### 3.1 Pré-processamento

1. **Leitura das abas da planilha NQO1** (`output/tables/nqo1/DGB_results_NQO1.xlsx`): 6 abas (`cmap/l1000/creeds` × `up/down`). Cada aba é lida, os nomes de coluna são padronizados (espaços → `_`, `-` → `_`, `.` → `_`) e recebem as colunas `source_database`, `source_sheet` e `expected_direction`.
2. **Conversão de tipo robusta** (`decimal_to_numeric`): valores com vírgula decimal europeia (`"1,23"`) são normalizados para ponto; `p_value`, `q_value` e `Fold_Change` tornam-se numéricos.
3. **Limpeza de `q-value = 0`**: valores de q exatamente zero (underflow numérico) são **capados** via `-log10(q)` limitado a 50 (`q_to_neglog10`), evitando distorção de escala no volcano plot.
4. **Imputação de missing — método:** **não é aplicada imputação** às variáveis contínuas, pois os dados analíticos têm **0% de NA** (ver auditoria). As colunas de ID específico de banco (`CREEDS_ID`, `GEO_ID`, `DrugBank_ID`) têm NA por "não se aplica", não por ausência. Caso futuras variáveis contínuas apresentem missing > 5%, o `config.yaml` define o método (`imputation.method: knn`, `k = 5`).
5. **Outliers — remoção/winsorização:** os outliers (|Z| > 3, `config.yaml` `outliers.z_score_threshold`) são **sinalizados, não removidos** por padrão. O ranking de fármacos usa estatística robusta (mínimo de q-value + escore), e os valores brutos permanecem em `data/raw/`. A ação configurável é `outliers.action: winsorize`.
6. **Unificação de nomes de fármacos (case-insensitive):** `Doxorubicin` e `doxorubicin` são unificados por `Drug_Name_Norm`, evitando contagem dupla no ranking.

### 3.2 Análise Estatística

- **Contraste de dois grupos (DCM vs. controle):** os testes de expressão diferencial foram conduzidos com o framework **limma** (colunas `t`, `pvalue`, `padj`, `B` em `DEG_full_table.csv`). A escolha de teste paramétrico/robusto segue a distribuição dos dados (limma aplica moderadores de variância empírica).
- **Correção de múltiplas comparações:** **Benjamini–Hochberg** (`statistics.multiple_testing_method: benjamini-hochberg`, FDR α = 0.05). Em `DEG_full_table.csv`, 485/41.842 genes são significativos (`padj < 0.05`).
- **Ranking NQO1 × fármacos:** cada fármaco é agregado por `mean_fc`, `median_fc`, `min_p`, `min_q`, `n_tests`, `n_sources` e um `score = -log10(q) × |FC|`. Ordenação por `min_q` decrescente de significância.
- **Tabela 1:** distribuição dos fármacos por plataforma (CMAP / LINCS L1000 / CREEDS / multi) × direção (up/down).

#### Resultados-chave
| Hipótese | Métrica | Valor | Veredito |
|---|---|---|---|
| NQO1 down em DCM | `log2FC` | −1.69 (padj 6.9e-13) | ✅ passa em FC>1.5 |
| Doxorrubicina → NQO1 | `mean_fc` | −1.07 (q 2.2e-28) | 🟡 significativo, efeito médio <1.5× |
| Via AMPK enriquecida | GSEA `NES` / `padj` | −1.19 / 0.419 | 🔴 **não significativo** |

> O enriquecimento GSEA da via AMPK **não** atinge FDR < 0.05. Qualquer afirmação sobre "ativação/inibição da AMPK" deve ser evitada ou claramente rotulada como exploratória.

### 3.3 Análise de Sensibilidade (estratificação)

- **Sexo:** 51 ♂ / 13 ♀ em `pheno_data_raw.csv` (desequilíbrio relevante). A análise de sensibilidade por sexo **ainda não está implementada**; recomenda-se um modelo ajustado por sexo/idade no DE (`~ 0 + group + sex + age`) para verificar se o efeito de NQO1 é robusto.
- **Idade:** faixa 20–66 anos (`age:ch1`). Recomenda-se análise de subgrupos por tercis de idade e por etiologia (DCM vs. ICM separadamente).
- **Subgrupo ICM:** o GSVA da via AMPK foi calculado **excluindo ICM** (51 = 37 DCM + 14 CTRL), enquanto o DEG cobre as 64 amostras. Documentar e justificar essa exclusão (atualmente não está no script).

### 3.4 Machine Learning (se aplicável)

- **Status: não aplicável ao estado atual** — não há variável-alvo rotulada de cardiotoxicidade/viabilidade (ver `config.yaml` `machine_learning.enabled: false`).
- **Plano:** se um ensaio de viabilidade (MTT/SRB) rotulado (`Doxo` × `Controle`) for incorporado, o fluxo previsto é:
  1. Modelo preditivo de cardiotoxicidade: **Random Forest** ou **Regressão Logística** (configurável em `config.yaml`).
  2. Importância de variáveis via **SHAP** (`shap`).
  3. Split treino/teste 75/25, seed 42.

---

## 4. Instruções de Reprodução

### 4.1 Pré-requisitos

**R ≥ 4.2** com os pacotes:

```r
install.packages(c("dplyr", "ggplot2", "stringr", "openxlsx", "igraph", "tidyr", "yaml"))
```

**Python ≥ 3.10** (auditoria/QC):

```bash
pip install -r requirements.txt
```

### 4.2 Ordem de execução

```bash
# 0) (opcional) clonar
git clone https://github.com/santosry/nqo1-doxorrubicina-cardiotoxicidade.git
cd nqo1-doxorrubicina-cardiotoxicidade

# 1) Pipeline principal (R) — gera tabelas e figuras
#    Windows (PowerShell):
& 'C:\Program Files\R\R-4.6.0\bin\Rscript.exe' script.R
#    Linux/macOS:
Rscript script.R

# 2) Auditoria de dados reprodutível (Python)
python scripts/audit_data.py

# 3) Relatório HTML de QC (rápido)
#    Opção A — RMarkdown:
Rscript -e 'rmarkdown::render("reports/QC_report.Rmd")'
#    Opção B — Quarto:
quarto render reports/QC_report.qmd
```

### 4.3 Consulta PubChem/NCBI (opcional, controlada por ambiente)

Por padrão, o script usa os resultados PubChem já salvos em `output/tables/pubchem/`. Para refazer a consulta:

```powershell
$env:RUN_PUBCHEM_QUERY='true'
$env:NCBI_API_KEY='sua_chave_ncbi_opcional'
```

> A chave NCBI **nunca** deve ser escrita no código — apenas no ambiente (e o `.gitignore` ignora arquivos de segredo).

### 4.4 Configuração centralizada

Todos os limiares estão em [`config.yaml`](config.yaml). Para alterar sem mexer no código:

```yaml
viability.doxo_alert_if_above: 0.8   # alerta se viabilidade Doxo > 0.8
expression.fold_change_significant: 1.5
outliers.z_score_threshold: 3.0
statistics.fdr_alpha: 0.05
```

Carregar no R: `cfg <- yaml::read_yaml("config.yaml")`; no Python: `cfg = yaml.safe_load(open("config.yaml"))`.

---

## 5. Logging e Auditoria

- O `script.R` registra cada etapa com **timestamp, nível e mensagem** em:
  - console (`message()`),
  - CSV estruturado `output/audit/audit_events_nqo1_pubchem_ppi.csv`,
  - arquivo `.log` diário em `output/audit/logs/run_YYYYMMDD.log` (registra timestamp, `INFO`/`WARN`/`ERROR`).
- O script auxiliar `scripts/audit_data.py` reproduz as checagens do [`AUDIT_REPORT.md`](AUDIT_REPORT.md).

---

## 6. Controle de Qualidade Visual (rápido)

Para gerar um relatório HTML com os principais gráficos de QC:

```bash
quarto render reports/QC_report.qmd   # ou rmarkdown::render("reports/QC_report.Rmd")
```

Gráficos de QC sugeridos:

1. **Boxplot antes/depois** da normalização (distribuição por amostra) — detecta escala/outliers.
2. **PCA** (PC1 × PC2 colorido por grupo) — detecta *batch effect* / agrupamento biológico.
3. **Densidade de expressão** por amostra — detecta amostras degradadas.
4. **Heatmap de correlação** entre amostras — detecta amostras trocadas/duplicadas.
5. **Barplot de % de genes detectados** por amostra.

---

## 7. Declaração de Transparência Algorítmica

> *"Este projeto utilizou a ferramenta de Inteligência Artificial generativa **deepseek-v4-pro** para auxiliar na automação de tarefas repetitivas de codificação (scripting), na estruturação e correção sintática da documentação, bem como na geração de rotinas de auditoria de dados. O uso desta ferramenta está em conformidade com a Portaria CNPq nº 2.664/2026, que regulamenta o uso de IA em pesquisas financiadas pela agência, garantindo a revisão crítica e validação final por um pesquisador titular (humano) responsável pelos resultados."*

---

## 8. Licença e Citação

- **Dados GEO:** `GSE116250` (Sweet et al., RNA-seq de miocárdio humano — DCM/ICM/não-falha).
- **STRING v11.5:** Szklarczyk et al.
- **Bancos de assinatura:** CMAP, LINCS L1000, CREEDS.
- Licença do código: definir conforme política institucional antes da publicação.
