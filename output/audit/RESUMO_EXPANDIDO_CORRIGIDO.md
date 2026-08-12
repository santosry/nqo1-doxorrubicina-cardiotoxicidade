# Análise Farmacogenômica in silico da Modulação do Gene NQO1 por Doxorrubicina na Cardiomiopatia Dilatada

**L.P. Souza¹\*; R.P. Santos¹; P.A. Borges¹**

¹Instituto Federal de Educação, Ciência e Tecnologia Fluminense Campus Campos Guarus

\*luciele.s@gsuite.iff.edu.br

---

## Resumo

Este estudo realizou uma análise transcriptômica in silico a partir do dataset GSE116250 (microarranjos de miocárdio humano com cardiomiopatia dilatada — CMD — versus controles sadios), complementada pela caracterização da via de sinalização AMPK e subsequente triagem farmacogenômica focada no gene NQO1 (NAD(P)H Quinona Desidrogenase 1). A análise de expressão diferencial global revelou que o gene NQO1 se encontra significativamente reprimido no miocárdio com CMD (log₂ Fold Change = −1,686; padj = 6,86 × 10⁻¹³). Paralelamente, a via AMPK não apresentou enriquecimento estatisticamente significativo no GSEA (p = 0,125; padj = 0,419), embora diversos genes constituintes da via tenham mostrado desregulação individual. O gene NQO1 foi selecionado para a etapa farmacogenômica com base em sua relevância funcional na citoproteção redox, independentemente de sua não filiação à via canônica AMPK do KEGG. A modelagem analítica integrativa por meio da plataforma Drug Gene Budger (DGB), integrando as bases CMAP, LINCS L1000 e CREEDS, totalizando 3523 registros válidos de 1281 fármacos únicos, demonstrou que a doxorrubicina atua como um significativo repressor transcricional de NQO1 (Fold Change médio = −2,565 no LINCS L1000 e −0,315 no CREEDS; q-value mínimo = 2,21 × 10⁻²⁸). Este achado indica que a antraciclina suprime a expressão de um gene antioxidante já deficitário na CMD, potencialmente exacerbando a vulnerabilidade ao estresse oxidativo. O mapeamento farmacogenômico identificou ainda outros fármacos com capacidade de downregulação de NQO1, como alvocidib (q-value = 1,14 × 10⁻¹¹) e panobinostat (q-value = 6,60 × 10⁻¹⁸), e upreguladores como trichostatin A (q-value < 10⁻⁵⁰) e parthenolide (q-value = 4,42 × 10⁻¹⁹). Conclui-se que a caracterização computacional deste circuito genético-farmacológico fornece subsídios para a compreensão dos mecanismos moleculares de cardiotoxicidade por antraciclinas, destacando a supressão de NQO1 como um fator que pode agravar o desequilíbrio redox preexistente na CMD. Estes achados in silico requerem validação experimental em modelos de cardiomiócitos primários antes de qualquer extrapolação clínica.

**Palavras-chave:** Biologia computacional; Cardiomiopatias; Farmacogenética; NQO1; Doxorrubicina.

---

## 1. Introdução

A Cardiomiopatia Dilatada (CMD) consiste em uma das principais etiologias causadoras de insuficiência cardíaca congestiva progressiva, distinguindo-se clinicamente pela dilatação do ventrículo esquerdo e grave disfunção sistólica. No plano molecular, a patogênese da enfermidade envolve um remodelamento metabólico complexo e crônico do miocárdio, caracterizado por disfunção mitocondrial e marcada exacerbação de estresse oxidativo tecidual [1,2].

A sinalização mediada pela Proteína Quinase Ativada por AMP (AMPK) funciona como um sensor energético celular mestre, coordenando respostas adaptativas microambientais ao estresse no miocárdio [4]. Em virtude disso, a caracterização do estado transcricional de genes da via AMPK no miocárdio patológico oferece uma perspectiva molecular relevante para a compreensão das alterações metabólicas associadas à CMD.

Paralelamente, a análise de expressão diferencial global do transcriptoma permite a identificação de genes significativamente desregulados, independentemente de sua filiação a vias canônicas específicas. Neste contexto, o gene NQO1 (NAD(P)H Quinona Desidrogenase 1) destacou-se como um dos genes mais significativamente alterados no miocárdio com CMD (padj = 6,86 × 10⁻¹³; log₂ Fold Change = −1,686). A seleção deste alvo para investigação farmacogenômica justifica-se pelo papel citoprotetor da enzima na regulação redox, atuando na redução bieletrônica de quinonas e hidroquinonas, o que previne a geração de espécies reativas de oxigênio (ERO) [5].

O gene NQO1 guarda estrita relação com mecanismos de cardiotoxicidade induzida por agentes quimioterápicos como a doxorrubicina, uma antraciclina eficaz no tratamento oncológico, mas cujo emprego clínico encontra forte limitação no risco dose-dependente de indução de cardiomiopatia severa provocada por danos mitocondriais e estresse oxidativo [7]. Investigar a arquitetura farmacogenômica celular e a modulação transcriptômica do gene NQO1 por xenobióticos contribui para a compreensão dos perfis moleculares de vulnerabilidade à cardiotoxicidade [6,8].

---

## 2. Materiais e Métodos

### 2.1. Materiais

Utilizaram-se matrizes de dados de expressão gênica públicas depositadas no repositório Gene Expression Omnibus (GEO, NCBI), com foco no conjunto de microarranjos indexado sob o código GSE116250, correspondente a amostras teciduais de miocárdio humano acometido por CMD (n = 37) em paralelo a controles sadios (n = 14). Para a análise de interações proteína-proteína, empregaram-se os dados da plataforma STRING versão 11.5 (Search Tool for the Retrieval of Interacting Genes, Homo sapiens), com escore de confiança mínimo de 700. Para a triagem farmacogenômica, utilizou-se o sistema gerenciador de assinaturas de perturbação Drug Gene Budger (DGB), integrando os bancos de dados químicos e biológicos do CMAP, LINCS L1000 e CREEDS [8].

### 2.2. Metodologia

A análise estatística de expressão diferencial dos genes foi realizada no ambiente R (versão 4.6.0) de computação por intermédio dos pacotes bioinformáticos GEOquery e limma, aplicando o ajuste de p-valor pela taxa de falsa descoberta (FDR, significância para padj < 0,05). A via de sinalização AMPK (hsa04152, KEGG) foi avaliada por GSEA (Gene Set Enrichment Analysis) e GSVA (Gene Set Variation Analysis). Concomitantemente, os dados de interações físicas e funcionais foram filtrados no STRING utilizando um escore de confiança restritivo mínimo de 700. A modelagem topológica da rede proteína-proteína (PPI) foi executada via pacote igraph, permitindo isolar a conectividade (degree) e projetar o perfil de log₂ Fold Change (logFC) dos nós.

Na etapa de mapeamento fármaco-gene, o gene NQO1 foi submetido à varredura de perturbações transcricionais no DGB, integrando assinaturas de exposição a xenobióticos das plataformas CMAP, LINCS L1000 e CREEDS. Os dados fármaco-gene foram padronizados limpando ruídos de formatação de strings e convertendo as dosagens químicas para a escala micromolar (µM). Para quantificar e ranquear a relevância estatística combinada à magnitude biológica da resposta induzida pelos compostos, computou-se um escore analítico customizado definido como: **score = −log₁₀(q-value) × |Fold Change|**.

---

## 3. Resultados e Discussão

### 3.1. Expressão diferencial e via AMPK

A análise de expressão diferencial global do transcriptoma (GSE116250) identificou alterações significativas em múltiplos genes no miocárdio com CMD. O gene NQO1 figurou como o segundo gene mais significativamente desregulado (log₂FC = −1,686; padj = 6,86 × 10⁻¹³), apresentando expressão reduzida em relação aos controles sadios — achado que sugere comprometimento das defesas antioxidantes dependentes de NQO1 no miocárdio patológico.

Paralelamente, genes constituintes da via AMPK (hsa04152, KEGG) apresentaram alterações individuais de expressão, com 46 genes contribuindo para o leading edge do GSEA. Contudo, a via AMPK como um todo **não apresentou enriquecimento estatisticamente significativo** (p = 0,125; padj = 0,419; NES = −1,19). A análise de GSVA confirmou uma tendência de redução dos escores da via AMPK no grupo CMD (média = −0,081) em comparação aos controles (média = +0,042), embora com sobreposição considerável entre os grupos. A rede PPI da via AMPK, construída a partir de 119 genes com 1874 interações de alta confiança no STRING, revelou alterações topológicas que envolvem hubs como PRKAB1 (degree = 74), PRKAB2 (degree = 72) e RPS6KB1 (degree = 64).

### 3.2. Triagem farmacogenômica do gene NQO1

O mapeamento fármaco-transcriptômico integrativo gerou um banco consolidado contendo 3523 assinaturas de perturbação biológica provocadas por xenobióticos, abrangendo 1281 fármacos únicos oriundos de três plataformas independentes (CMAP, LINCS L1000 e CREEDS), em condições experimentais variadas de concentração (0,1 µM a 1000 µM), tempo de exposição (6 h a 24 h) e linhagens celulares.

Os resultados numéricos derivados do modelo de triagem analítica para os cinco primeiros colocados no ranking molecular encontram-se sumarizados na **Tabela 1**.

**Tabela 1.** Top 5 fármacos moduladores de NQO1 segundo o ranking integrado.

| Posição | Fármaco | q-value mínimo | Fold Change médio | Score | Direção | Fonte |
|---------|---------|---------------|-------------------|-------|---------|-------|
| 1º | trichostatin A | < 10⁻⁵⁰* | +0,681 | 60,7 | Upregula NQO1 | CMAP |
| 2º | Doxorubicina | 2,21 × 10⁻²⁸ | −0,315 | 41,4 | Downregula NQO1 | CREEDS; LINCS L1000 |
| 3º | manumycin-a | 3,62 × 10⁻²⁰ | −0,010 | 68,1 | Upregula NQO1† | LINCS L1000 |
| 4º | parthenolide | 4,42 × 10⁻¹⁹ | +1,569 | 94,0 | Upregula NQO1 | LINCS L1000 |
| 5º | panobinostat | 6,60 × 10⁻¹⁸ | −0,539 | 24,9 | Downregula NQO1 | LINCS L1000 |

\* q-value igual a zero foi capado em 10⁻⁵⁰ para evitar distorção de escala (−log₁₀(0) = ∞).
† manumycin-a apresenta mean_fc ligeiramente negativo (−0,010) mas median_fc positivo (+0,746), indicando distribuição assimétrica dos efeitos; a direção foi determinada pela mediana.

A modelagem gráfica por Volcano plot (Figura 1) evidenciou a distribuição global dos fármacos, com trichostatin A e parthenolide destacando-se como upreguladores de NQO1, enquanto doxorrubicina, alvocidib (q-value = 1,14 × 10⁻¹¹; Fold Change médio = −1,485; posição 18ª no ranking), panobinostat e wortmannin (q-value = 6,74 × 10⁻⁸; Fold Change médio = −0,159; posição 58ª) figuram como downreguladores significativos.

A análise de dose-resposta (Figura 2), restrita a fármacos com pelo menos duas concentrações numéricas distintas, demonstrou que a doxorrubicina exibe perfil de downregulação consistente na concentração de 10 µM nas plataformas LINCS L1000 e CREEDS.

### 3.3. Implicações da downregulação de NQO1 pela doxorrubicina

A doxorrubicina demonstrou associação estatisticamente significante como **repressor transcricional** de NQO1, com q-value mínimo de 2,21 × 10⁻²⁸ (CREEDS) e Fold Change médio de −2,565 (LINCS L1000) a −0,315 (CREEDS). Sob a ótica da plausibilidade biológica, essa supressão transcriptômica sugere que a doxorrubicina reduz os níveis de NQO1 — enzima chave na defesa antioxidante celular — em um contexto onde o gene já se encontra significativamente reprimido pela própria CMD (log₂FC = −1,686). A formação intracelular de radicais semiquinona decorrentes do metabolismo da doxorrubicina, combinada à capacidade reduzida de detoxificação pela via NQO1, pode contribuir para um ciclo de exacerbação do estresse oxidativo e consequente disfunção mitocondrial [3,7].

É relevante notar que a direção da regulação (downregulation) é consistente entre duas plataformas independentes (CREEDS e LINCS L1000), conferindo maior confiabilidade ao achado.

Como limitação metodológica manifesta, destaca-se o caráter estritamente in silico das matrizes do DGB, coletadas majoritariamente em modelos celulares tumorais (MCF7, PC3, VCAP), cujas dinâmicas regulatórias podem diferir substancialmente daquelas de cardiomiócitos humanos primários in vivo submetidos a estresse mecânico e hemodinâmico. Desta forma, os achados aqui apresentados não permitem inferência causal direta ou transposição imediata de magnitudes farmacológicas para contextos clínicos sem a prévia validação experimental em modelos cardíacos.

Não obstante tais limitações, os achados assumem relevância translacional potencial. O desvelamento destas interações sugere que a supressão de NQO1 pela doxorrubicina pode representar um mecanismo adicional de cardiotoxicidade, particularmente em pacientes com CMD cujos níveis basais de NQO1 já se encontram reduzidos. Estudos futuros de validação funcional em cardiomiócitos primários e a investigação de polimorfismos de NQO1 (como o alelo NQO1\*2, de atividade enzimática reduzida) em coortes clínicas são necessários para validar e expandir estes achados.

---

## 4. Conclusões

A análise transcriptômica em ambiente computacional revelou que o gene NQO1 se encontra significativamente reprimido no miocárdio de pacientes com CMD (log₂FC = −1,686; padj = 6,86 × 10⁻¹³). A via de sinalização AMPK, embora apresente alterações em genes individuais, não demonstrou enriquecimento global estatisticamente significativo (p = 0,125), indicando que as alterações metabólicas na CMD envolvem mecanismos adicionais além desta via.

A modelagem farmacogenômica integrativa revelou que a doxorrubicina atua como um **repressor transcricional** da expressão de NQO1 (q-value mínimo = 2,21 × 10⁻²⁸), indicando que esta antraciclina suprime um gene antioxidante já deficitário no miocárdio com CMD. Este achado sugere um potencial mecanismo de exacerbação do estresse oxidativo mediado pela combinação entre a repressão patológica de NQO1 na CMD e sua supressão farmacológica pela doxorrubicina. O ranking farmacogenômico completo (1281 fármacos) e as análises de rede PPI (119 genes, 1874 interações) constituem um recurso para investigações futuras de fármacos moduladores de NQO1. Embora requeiram corroboração em sistemas celulares miocárdicos primários, os dados mapeados elucidam interfaces farmacogenômicas relevantes que fornecem subsídios conceituais para a compreensão dos mecanismos moleculares de cardiotoxicidade por antraciclinas.

---

## Agradecimentos

Ao Instituto Federal Fluminense (IFFluminense) e ao Núcleo de Pesquisa e Estudos em Saúde (NUPES) pelo suporte institucional e técnico para a execução desta pesquisa.

---

## Referências

[1] BOZKURT, B. Dilated cardiomyopathy: a review of current and future treatment strategies. **Nature Reviews Cardiology**, v. 18, p. 673–694, 2021.

[2] RAMACCINI, D.; MONTOYA-URIBE, V.; AAN, F. J. et al. Mitochondrial Function and Dysfunction in Dilated Cardiomyopathy. **Frontiers in Cell and Developmental Biology**, v. 8, 2021.

[3] MARTINS, D. et al. Oxidative Stress as a Therapeutic Target of Cardiac Remodeling. **Antioxidants**, v. 11, n. 12, 2022.

[4] HERZIG, S.; SHAW, R. J. AMPK: guardian of metabolism and mitochondrial homeostasis. **Nature Reviews Molecular Cell Biology**, v. 19, p. 121–135, 2018.

[5] ROSS, D.; SIEGEL, D. Functions of NQO1 in cellular protection and CoQ10 metabolism and its potential role as a redox sensitive molecular switch. **Frontiers in Physiology**, v. 8, 2017.

[6] HODOS, R. A.; KIDD, B. A.; SHAMEER, K. et al. In silico methods for drug repurposing and pharmacology. **Wiley Interdisciplinary Reviews: Systems Biology and Medicine**, v. 8, n. 3, p. 186–210, 2016.

[7] WALLACE, K. B. Doxorubicin-induced cardiac mitochondrionopathy. **Pharmacology & Toxicology**, v. 93, n. 3, p. 105–115, 2003.

[8] WANG, Z. Y. et al. Drug Gene Budger (DGB): an application for ranking drugs to modulate a specific gene based on transcriptomic signatures. **Bioinformatics**, v. 37, n. 8, p. 1247–1248, 2021.

[9] EKINS, S.; PUHL, A. C.; ZORN, K. M. et al. Exploiting machine learning for end-to-end drug discovery and development. **Nature Materials**, v. 18, n. 5, p. 435–441, 2019.

[10] ZANGIROLAMI-RAIMUNDO, J.; ECHEIMBERG, J. O.; LEONE, C. Tópicos de metodologia de pesquisa: estudos de corte transversal. **Journal of Human Growth and Development**, v. 28, n. 3, p. 356–360, 2018.

[11] COSTA, M. A. M. S. Resultados e discussões integradas de bioinformática e dados transcriptômicos. In: Resultados Esperados de Projeto de Iniciação Científica, IFFluminense, 2026.

---

## Nota de auditoria (2026-06-26)

Esta versão do manuscrito foi integralmente corrigida após auditoria forense dos dados e do código-fonte (script.R). As principais correções aplicadas foram:

1. **Título** alterado para refletir os achados reais (sem subtítulo)
2. **q-value da doxorrubicina** corrigido de 6,32 × 10⁻⁶ (inexistente) para 2,21 × 10⁻²⁸ (CREEDS) e 2,50 × 10⁻⁹ (LINCS L1000)
3. **Direção da regulação** corrigida: doxorrubicina, geldanamicina e wortmannin são DOWNREGULADORES de NQO1 (anteriormente descritos incorretamente como upreguladores)
4. **Monensina** removida — composto não consta nos 1281 fármacos do ranking
5. **Relação AMPK → NQO1** corrigida: NQO1 foi identificado na análise do transcriptoma completo, não como parte da via AMPK (à qual não pertence); GSEA da via AMPK não significativo agora reportado
6. **Referência [13]** removida — era de engenharia de materiais (carbono pirolítico), sem relação com o tema
7. **Discussão e Conclusão** reescritas com base nos dados reais
8. **Ranking** unificado com normalização case-insensitive ("Doxorubicin" e "doxorubicin" fundidos)
9. **Código-fonte** (script.R) corrigido para normalização case-insensitive, remoção de código morto e documentação de distribuições assimétricas
10. **Tipografia** corrigida ("Institudo" → "Instituto")
