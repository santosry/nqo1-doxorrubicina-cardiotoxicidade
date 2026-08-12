# TODO — Roadmap de Análises Futuras

> Priorização: **P1** = validação imediata do achado central · **P2** = robustez e contexto · **P3** = expansão/exploração.

Cada item descreve: objetivo biológico, método, dados necessários (e fonte) e código-base para começar.

---

## P1 — Validação independente da downregulação de NQO1 em outra coorte de DCM

**Objetivo biológico:** confirmar que a redução de NQO1 na DCM não é um artefato do GSE116250.

**Método:** expressão diferencial (limma) em coortes independentes + meta-análise (efeito aleatório) dos `log2FC`.

**Dados:** GEO (GSE57338, GSE120895, GSE71613, GSE3585, GSE42955) ou GTEx (Heart — Left Ventricle).

**Código-base (R):**
```r
library(GEOquery); library(limma)
gse <- getGEO("GSE57338", GSEMatrix = TRUE)[[1]]
expr <- exprs(gse)
# (adaptar: construir pheno com grupo DCM vs controle e rodar limma)
fit <- lmFit(expr, model.matrix(~ group, data = pheno)) |> eBayes()
topTable(fit, coef = 2) |> dplyr::filter(gene == "NQO1")
```

---

## P1 — Docking molecular NQO1 × Doxorrubicina

**Objetivo biológico:** avaliar a plausibilidade estrutural de interação direta entre doxorrubicina e NQO1.

**Método:** docking molecular (AutoDock Vina), energia de ligação (kcal/mol), visualização (PyMOL).

**Dados:** estrutura da NQO1 (PDB `1D4A`, `2F1O` ou `1QBG`) · doxorrubicina (PubChem CID `31703` → SDF).

**Código-base (Python):**
```python
# 1) baixar PDB + ligante e preparar (remove água/ligantes, adiciona H)
# 2) definir caixa de docking no sítio ativo (coordenadas do FAD)
# 3) executar o Vina
import subprocess
subprocess.run([
    "vina", "--receptor", "nqo1_prep.pdbqt",
    "--ligand", "doxorubicin.pdbqt",
    "--center_x", "10", "--center_y", "20", "--center_z", "30",
    "--size_x", "25", "--size_y", "25", "--size_z", "25",
    "--out", "doxo_nqo1_out.pdbqt", "--exhaustiveness", "32"
])
```

---

## P1 — Validação funcional in vitro (viabilidade MTT/SRB + qPCR de NQO1)

**Objetivo biológico:** validar experimentalmente que a doxorrubicina (i) suprime NQO1 e (ii) reduz viabilidade — acionando o threshold já pronto em `config.yaml`.

**Método:** curva dose-resposta (ajuste 4PL), qRT-PCR (ΔΔCt), teste t/Wilcoxon entre Doxo × Controle.

**Dados:** a serem gerados em laboratório (cardiomiócitos H9c2/AC16).

**Código-base (R):**
```r
library(drc)
fit <- drm(viability ~ dose, data = assay, fct = LL.4())
ED50 <- ED(fit, 50)
# NQO1 por qPCR (delta-delta-Ct)
t.test(nqo1_dct ~ group, data = qpcr)   # ou wilcox.test
```

---

## P2 — Associação genética: NQO1*2 (rs1800566) e eQTL

**Objetivo biológico:** ligar a variação genética de NQO1 (alelo *2, atividade reduzida) à sua expressão e à suscetibilidade à cardiotoxicidade.

**Método:** consulta de eQTL (GTEx), colocalização, anotação farmacogenética (PharmGKB).

**Dados:** GTEx V8 (Heart — Left Ventricle), PharmGKB (rs1800566).

**Código-base (R):**
```r
# GTEx Portal API (exemplo simplificado)
httr::GET("https://gtexportal.org/api/v2/association/singleTissueEqtl",
          query = list(geneId = "NQO1", tissueSiteDetailId = "Heart_Left_Ventricle"))
```

---

## P2 — Correlação clínica: NQO1 × fração de ejeção / desfechos

**Objetivo biológico:** verificar se a expressão de NQO1 correlaciona com a função cardíaca (FE) ou com desfechos clínicos.

**Método:** correlação de Spearman/Pearson, regressão linear e análise de sobrevida (Cox).

**Dados:** coortes GEO com metadados clínicos (FE) ou uma coorte de cardiomiopatia.

**Código-base (R):**
```r
library(survival)
cor.test(nqo1_expr, ejection_fraction, method = "spearman")
coxph(Surv(time, event) ~ nqo1_expr + age + sex, data = clin)
```

---

## P2 — Localização célula-específica de NQO1 (scRNA-seq cardíaco)

**Objetivo biológico:** identificar em qual tipo celular (cardiomiócito, fibroblasto, endotélio) ocorre a mudança de NQO1.

**Método:** análise de scRNA-seq (Seurat/Scanpy), DE por tipo celular.

**Dados:** atlas de coração humano (ex.: Litviňuková et al. 2020), Human Cell Atlas.

**Código-base (R):**
```r
library(Seurat)
obj <- Seurat::NormalizeData(obj) |> Seurat::FindVariableFeatures() |>
  Seurat::ScaleData() |> Seurat::RunPCA() |> Seurat::FindNeighbors() |> Seurat::FindClusters()
Seurat::FeaturePlot(obj, features = "NQO1")
```

---

## P3 — Modelo preditivo de cardiotoxicidade + SHAP

**Objetivo biológico:** construir um classificador de fármacos cardiotóxicos a partir de assinaturas transcriptômicas.

**Método:** Random Forest / Regressão Logística + importância via SHAP.

**Dados:** assinaturas LINCS L1000 + rótulos de cardiotoxicidade (literatura/CTD).

**Código-base (Python):**
```python
from sklearn.ensemble import RandomForestClassifier
import shap
model = RandomForestClassifier(n_estimators=500, random_state=42).fit(X_train, y_train)
explainer = shap.TreeExplainer(model)
shap_values = explainer.shap_values(X_test)
shap.summary_plot(shap_values, X_test)
```

---

## P3 — Expansão para o painel antioxidante NRF2/ARE (GCLC, HMOX1, TXNRD1, NQO1)

**Objetivo biológico:** estender a triagem DGB para todo o eixo NRF2/ARE, contextualizando NQO1 no programa antioxidante.

**Método:** consulta DGB multi-gene + agregação de ranks (reuso das funções de `script.R`).

**Dados:** DGB (CMAP/LINCS/CREEDS) para os genes do painel.

**Código-base (R):**
```r
# Reusar read_nqo1_xlsx()/rank_nqo1_drugs() de script.R para cada gene-alvo,
# depois agregar: bind_rows(lista_rankings) |> group_by(Drug_Name_Norm) |> summarise(n_genes)
```

---

## Resumo de prioridades

| # | Análise | Prioridade |
|---|---|---|
| 1 | Validação NQO1 em coorte independente de DCM | P1 |
| 2 | Docking molecular NQO1 × Doxorrubicina | P1 |
| 3 | Validação in vitro (viabilidade + qPCR) | P1 |
| 4 | NQO1*2 / eQTL (rs1800566) | P2 |
| 5 | NQO1 × fração de ejeção / sobrevida | P2 |
| 6 | scRNA-seq (localização celular de NQO1) | P2 |
| 7 | ML preditivo + SHAP | P3 |
| 8 | Painel antioxidante NRF2/ARE | P3 |
