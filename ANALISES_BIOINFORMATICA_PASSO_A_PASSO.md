# Análises Bioinformáticas — Passo a Passo (somente dados secundários)

> Guia prático para as 7 análises que podem ser executadas **apenas com programação e dados públicos** (sem bancada).
> A análise nº 3 do roadmap (validação in vitro MTT/SRB + qPCR) foi **excluída** por exigir dados experimentais primários.

---

## Análise 1 — Validação independente de NQO1 em outra coorte de DCM

**Objetivo biológico:** confirmar que a downregulação de NQO1 não é um artefato do GSE116250.

**Dados (secundários):**
- GEO: `GSE57338`, `GSE120895`, `GSE71613`, `GSE3585` (DCM vs controle).
- GTEx (Heart – Left Ventricle) como tecido normal de referência.

**Ferramentas:** R (`GEOquery`, `limma`, `dplyr`).

**Passo a passo:**
```r
# 1) Instalar (uma vez)
if (!requireNamespace("BiocManager", quietly=TRUE)) install.packages("BiocManager")
BiocManager::install(c("GEOquery","limma"))

# 2) Baixar uma coorte
library(GEOquery); library(limma)
gse <- getGEO("GSE57338", GSEMatrix = TRUE)[[1]]   # ExpressionSet
expr <- exprs(gse)                                  # matriz
pheno <- pData(gse)                                 # metadados

# 3) Identificar a coluna do grupo (DCM vs controle) — inspecionar pheno
#    (adaptar o nome da coluna conforme o dataset)
group <- pheno$`disease state:ch1`  # exemplo; ajustar
group <- ifelse(grepl("control|non-failing", group, ignore.case=TRUE), "CTRL", "DCM")

# 4) Expressão diferencial
design <- model.matrix(~ 0 + factor(group))
colnames(design) <- c("CTRL","DCM")
fit <- eBayes(contrasts.fit(lmFit(expr, design), makeContrasts(DCM - CTRL, levels=design)))
res <- topTable(fit, number=Inf, sort.by="none")

# 5) Verificar NQO1
res[rownames(res) == "NQO1", c("logFC","adj.P.Val")]
```

**Critério de sucesso:** NQO1 com `logFC < 0` e `adj.P.Val < 0.05` em ≥ 1 coorte independente.

---

## Análise 2 — Docking molecular NQO1 × Doxorrubicina

**Objetivo biológico:** avaliar se a doxorrubicina pode interagir fisicamente com a NQO1 (plausibilidade estrutural).

**Dados (secundários):**
- NQO1 (PDB): `1D4A`, `2F1O` ou `1QBG`.
- Doxorrubicina (PubChem): CID `31703` → `.sdf`.

**Ferramentas:** AutoDock Vina, OpenBabel, PyMOL, Python (`biopython`, `meeko`).

**Passo a passo:**
```bash
# 1) Instalar (conda recomendado)
conda create -n docking -c conda-forge autodock-vina openbabel pymol-open-source -y
conda activate docking
```
```python
# 2) Preparar receptor (remover água/ligantes, adicionar hidrogênios)
from Bio.PDB import PDBParser, PDBIO, Select
# baixar PDB: https://files.rcsb.org/download/1D4A.pdb
# remover HETATM (água/FAD do cristal, mantendo só a proteína) e salvar receptor.pdb

# 3) Preparar ligante (SDF -> PDBQT)
import subprocess
subprocess.run(["obabel", "doxorubicin.sdf", "-O", "doxorubicin.pdbqt",
                "--gen3d", "-p", "7.4"])   # protonar em pH fisiológico

# 4) Definir a caixa de docking (sítio ativo da NQO1 — coordenadas do FAD)
#    (obter centro/tamanho do cofator no PyMOL: selecionar FAD e "get_extent")
# 5) Executar o docking
subprocess.run([
    "vina", "--receptor", "nqo1_receptor.pdbqt",
    "--ligand", "doxorubicin.pdbqt",
    "--center_x", "10", "--center_y", "20", "--center_z", "30",  # ajustar
    "--size_x", "25", "--size_y", "25", "--size_z", "25",
    "--exhaustiveness", "32", "--num_modes", "9",
    "--out", "doxo_nqo1.pdbqt"
])
```
**Critério:** afinidade de ligação ≤ −7 kcal/mol sugere interação plausível; validar visualmente no PyMOL.

---

## Análise 4 — NQO1*2 (rs1800566) e eQTL

**Objetivo biológico:** ligar a variação genética de NQO1 (alelo *2 = atividade reduzida) à sua expressão cardíaca e à suscetibilidade à cardiotoxicidade.

**Dados (secundários):** GTEx Portal (eQTL), PharmGKB, dbSNP.

**Ferramentas:** Python (`requests`) ou R (`httr`).

**Passo a passo:**
```python
import requests
# 1) eQTL da NQO1 no ventrículo esquerdo (GTEx v8)
url = "https://gtexportal.org/api/v2/association/singleTissueEqtl"
params = {
    "gencodeId": "ENSG00000181019",           # NQO1
    "tissueSiteDetailId": "Heart_Left_Ventricle",
    "datasetId": "gtex_v8"
}
r = requests.get(url, params=params)
data = r.json()
# 2) Filtrar o SNP rs1800566 (NQO1*2) nos resultados de eQTL
snp = [x for x in data["data"] if x.get("snpId") == "rs1800566"]
print(snp)
```
**Critério:** verificar se `rs1800566` é um eQTL de NQO1 (beta/valor-p) e anotar a direção do efeito no PharmGKB (`https://www.pharmgkb.org/variant/PA166153558`).

---

## Análise 5 — NQO1 × fração de ejeção / sobrevida (Cox)

**Objetivo biológico:** correlacionar a expressão de NQO1 com a função cardíaca (FE) ou desfechos clínicos.

**Dados (secundários):** coortes GEO com metadados clínicos (ex.: FE) — **condicional** à disponibilidade desses campos.

**Ferramentas:** R (`survival`, `dplyr`).

**Passo a passo:**
```r
library(survival)
# 1) Obter expressão de NQO1 + metadados clínicos (FE, tempo, evento) de uma coorte GEO
# 2) Correlação com fração de ejeção
cor.test(nqo1_expr, ejection_fraction, method = "spearman")

# 3) Análise de sobrevida (Cox)
cox <- coxph(Surv(time, event) ~ nqo1_expr + age + sex, data = clin)
summary(cox)
```
**Critério:** correlação significativa (p < 0.05) e/ou HR ≠ 1 com p < 0.05 no Cox.

---

## Análise 6 — scRNA-seq cardíaco (localização celular de NQO1)

**Objetivo biológico:** identificar em qual tipo celular (cardiomiócito, fibroblasto, endotélio) a NQO1 é expressa/alterada.

**Dados (secundários):** atlas de coração humano (ex.: Litviňuková et al. 2020, GEO `GSE161157`; Human Heart Cell Atlas).

**Ferramentas:** R (`Seurat`) ou Python (`scanpy`).

**Passo a passo:**
```r
library(Seurat)
# 1) Baixar matriz (10x) e criar objeto
obj <- Read10X("path/to/filtered_feature_bc_matrix")
obj <- CreateSeuratObject(obj)

# 2) QC + normalização + clustering
obj <- subset(obj, subset = nFeature_RNA > 200 & percent.mt < 10)
obj <- NormalizeData(obj) |> FindVariableFeatures() |> ScaleData() |>
  RunPCA() |> FindNeighbors() |> FindClusters()

# 3) Anotar tipos celulares (marcadores: TNNI3/MYH7 = cardiomiócito, DCN = fibroblasto, PECAM1 = endotélio)
# 4) Visualizar NQO1
FeaturePlot(obj, features = "NQO1")
VlnPlot(obj, features = "NQO1", group.by = "celltype")
```
**Critério:** identificar o(s) tipo(s) celular(es) com maior expressão de NQO1.

---

## Análise 7 — ML preditivo de cardiotoxicidade + SHAP

**Objetivo biológico:** construir um classificador de fármacos cardiotóxicos a partir de assinaturas transcriptômicas.

**Dados (secundários):** assinaturas LINCS L1000 (clue.io/GEO) + rótulos de cardiotoxicidade (literatura, CTD).

**Ferramentas:** Python (`pandas`, `scikit-learn`, `shap`).

**Passo a passo:**
```python
import pandas as pd
from sklearn.ensemble import RandomForestClassifier
from sklearn.model_selection import train_test_split
import shap

# 1) Construir X (assinaturas L1000 por fármaco) e y (0 = não cardiotóxico, 1 = cardiotóxico)
X = pd.read_csv("l1000_signatures.csv")   # genes x fármacos (transposta)
y = pd.read_csv("labels.csv")["cardiotoxic"]

# 2) Treino/teste
X_train, X_test, y_train, y_test = train_test_split(X, y, test_size=0.25, random_state=42)
model = RandomForestClassifier(n_estimators=500, random_state=42).fit(X_train, y_train)

# 3) Importância (SHAP)
explainer = shap.TreeExplainer(model)
shap_values = explainer.shap_values(X_test)
shap.summary_plot(shap_values, X_test)   # genes mais relevantes
```
**Critério:** AUC > 0.8 em validação cruzada e identificar os genes top por SHAP.

---

## Análise 8 — Painel antioxidante NRF2/ARE (GCLC, HMOX1, TXNRD1, NQO1) no DGB

**Objetivo biológico:** estender a triagem DGB para todo o eixo NRF2/ARE, contextualizando a NQO1 no programa antioxidante.

**Dados (secundários):** Drug Gene Budger (CMAP, LINCS L1000, CREEDS).

**Ferramentas:** R (reuso das funções `read_nqo1_xlsx()`/`rank_nqo1_drugs()` do `script.R`).

**Passo a passo:**
```r
# 1) Para cada gene do painel, baixar o XLSX no DGB (https://drugbudger.com ou equivalente)
#    e salvar como DGB_results_<GENE>.xlsx
# 2) Reusar as funções do script.R (source com cuidado, ou copiar as funções)
genes <- c("NQO1","GCLC","HMOX1","TXNRD1","KEAP1","NFE2L2")
rankings <- lapply(genes, function(g) {
  d <- read_nqo1_xlsx(sprintf("DGB_results_%s.xlsx", g))
  rank_nqo1_drugs(d) |> dplyr::mutate(target_gene = g)
})
# 3) Agregar: fármacos que modulam MÚLTIPLOS genes do eixo antioxidante
dplyr::bind_rows(rankings) |>
  dplyr::group_by(Drug_Name_Norm) |>
  dplyr::summarise(n_genes = dplyr::n_distinct(target_gene),
                   genes = paste(sort(unique(target_gene)), collapse="; ")) |>
  dplyr::filter(n_genes >= 2) |> dplyr::arrange(dplyr::desc(n_genes))
```
**Critério:** fármacos que modulam ≥ 2 genes do painel NRF2/ARE são candidatos a moduladores do programa antioxidante como um todo.

---

## Ordem recomendada de execução

| Prioridade | Análise | Dificuldade | Tempo estimado |
|---|---|---|---|
| 1 | #1 Validação independente | Baixa | 1–2 h |
| 1 | #2 Docking molecular | Média | 1 dia |
| 1 | #4 eQTL (rs1800566) | Baixa | 1 h |
| 2 | #6 scRNA-seq | Média-alta | 1–2 dias |
| 2 | #8 Painel NRF2/ARE | Baixa-média | 2–3 h |
| 2 | #7 ML + SHAP | Média | 1 dia |
| 3 | #5 FE/sobrevida | Média (depende de dados) | variável |
