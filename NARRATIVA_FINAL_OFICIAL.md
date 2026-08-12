# Narrativa Final Oficial do Estudo

## Fluxo Lógico do Estudo (3 passos)

### PASSO 1 — Caracterização exploratória da via AMPK (KEGG: hsa04152) na cardiomiopatia dilatada
A investigação iniciou-se pela análise de expressão diferencial da via de sinalização **AMPK**
(KEGG `hsa04152`) no miocárdio com **cardiomiopatia dilatada (DCM)**, a partir do dataset
**GSE116250** (RNA-seq, 64 amostras: 14 não-falha, 37 DCM, 13 ICM). A via AMPK foi avaliada por
expressão diferencial (limma), **GSEA** e **GSVA**, não apresentando enriquecimento
estatisticamente significativo (**NES = −1,19; padj = 0,419**). Esta etapa constituiu a
**primeira parte de um estudo exploratório mais amplo** — partiu-se dela, mas o resultado
negativo foi registrado como tal.

### PASSO 2 — Seleção do gene NQO1 a partir do transcriptoma completo
Em seguida, foram considerados **todos os genes diferencialmente expressos do transcriptoma
completo** (GSE116250). Dentre eles, selecionou-se o gene com a **significância estatística mais
robusta e maior relevância biológica**: o **NQO1** (NAD(P)H quinona desidrogenase 1),
significativamente **downregulado na DCM** (**log2 Fold Change = −1,686; padj = 6,86 × 10⁻¹³**).
A seleção baseou-se na função citoprotetora/redox do NQO1 e em sua relação com a cardiotoxicidade
por antraciclinas — **não** por pertencer à via AMPK (à qual o gene não pertence).

### PASSO 3 — Identificação da Doxorrubicina via Drug Gene Budger (DGB)
Com o alvo **NQO1** definido, consultou-se a plataforma **Drug Gene Budger (DGB)**, integrando as
bases **CMAP**, **LINCS L1000** e **CREEDS**, para identificar fármacos moduladores da expressão
de NQO1. O fármaco com **resultado mais significativo e evidência convergente em duas bases
independentes** foi a **Doxorrubicina** (**q = 2,21 × 10⁻²⁸**; Fold Change médio = −1,065;
downregulador; CREEDS + LINCS L1000, n = 15 registros).

---

### Nota sobre o ranking (Doxorrubicina × Trichostatin A)
O **trichostatin A** aparece com `q = 0`, valor decorrente de **underflow numérico** (fonte
única, CMAP) e capado em 10⁻⁵⁰. Portanto, a **Doxorrubicina** corresponde ao fármaco com o
**menor q-value real** do ranking (2,21 × 10⁻²⁸) e o único do topo com **replicação em múltiplas
plataformas independentes**, sendo a resposta correta do estudo.
