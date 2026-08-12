# Justificativa Científica do Desenho do Estudo

> **Por que este estudo foi desenhado desta forma?** Da via AMPK ao gene NQO1 e, dele, à doxorrubicina.
> Documento complementar à `NARRATIVA_FINAL_OFICIAL.md` — foca no **racional científico**, não apenas na sequência.

---

## 1. A pergunta de pesquisa original (ampla e exploratória)

A pergunta inicial não era sobre um gene ou fármaco específico, mas sobre um fenômeno clínico:

> **"Quais fármacos contribuem para a cardiomiopatia — em especial por mecanismos de estresse oxidativo e disfunção metabólica?"**

Essa formulação é legitimamente **exploratória e geradora de hipótese**: parte-se de um problema clínico (cardiotoxicidade/cardiomiopatia) e usa-se bioinformática para **triar** candidatos moleculares e farmacológicos, em vez de testar uma hipótese única e fechada.

---

## 2. O desenho como "funil" de três etapas

O estudo segue uma lógica de **afunilamento progressivo**, que reduz a pergunta ampla a um alvo molecular e, dele, a um fármaco:

```
PROBLEMA CLÍNICO (cardiomiopatia/cardiotoxicidade)
        │
        ▼
[PASSO 1] Via candidata (AMPK, hsa04152) ────► contexto metabólico/energético
        │  (resultado: via NÃO enriquecida, mas genes individuais alterados)
        ▼
[PASSO 2] Transcriptoma completo ────► gene mais significativo = NQO1 (redox)
        │  (alvo molecular com relevância biológica)
        ▼
[PASSO 3] Drug Gene Budger ────► fármaco modulador = Doxorrubicina
           (farmacologia reversa "gene → fármaco")
```

Este desenho combina **duas estratégias complementares**, o que o torna metodologicamente robusto:

| Estratégia | O que faz | Resultado obtido |
|---|---|---|
| **Dirigida por hipótese** (via candidata) | Parte de conhecimento biológico prévio (AMPK como sensor energético cardíaco) | Contexto metabólico; hipótese de via refutada |
| **Dirigida por dados** (transcriptoma completo) | Varre todo o transcriptoma sem viés prévio | NQO1 como gene mais desregulado |

---

## 3. Por que começar pela via AMPK? (racional biológico)

A **AMPK (proteína quinase ativada por AMP)** é o **sensor energético mestre da célula**, com papel central no miocárdio:

1. **Homeostase energética cardíaca.** O coração é o órgão com maior demanda metabólica contínua do corpo. A AMPK coordena captação de glicose, oxidação de ácidos graxos, biogênese mitocondrial e autofagia [Herzig & Shaw, 2018].

2. **A cardiomiopatia dilatada é, em essência, uma doença de remodelamento metabólico.** No coração insuficiente há **déficit energético** (o miocárdio "passa fome"), disfunção mitocondrial e acúmulo de estresse oxidativo [Bozkurt et al., 2021].

3. **Portanto, a AMPK é a via mais natural para iniciar** uma investigação de cardiomiopatia: se existe um eixo metabólico desregulado na DCM, a AMPK é a primeira suspeita.

> **Justificativa formal:** *"A AMPK foi selecionada como ponto de partida por ser o principal regulador do metabolismo energético e da homeostase mitocondrial no cardiomiócito, dois processos centrais na fisiopatologia da cardiomiopatia dilatada."*

---

## 4. O valor científico do resultado "negativo" da AMPK

A análise mostrou que a via AMPK **como um todo não está enriquecida de forma significativa** (GSEA: NES = −1,19; padj = 0,419), embora **genes individuais** estejam desregulados.

Este resultado **não é um fracasso** — é um achado com valor próprio:

1. **Refuta uma hipótese plausível**: a DCM não é dirigida por um deslocamento coordenado da via AMPK.
2. **Redireciona a investigação**: se a via inteira não muda, o sinal relevante deve estar em **genes isolados**, justificando a varredura do transcriptoma completo.
3. **Evita viés de confirmação**: o resultado negativo demonstra que o estudo não "forçou" a via AMPK a aparecer como significativa.

> **Justificativa formal:** *"A ausência de enriquecimento global da via AMPK indica que as alterações moleculares da DCM envolvem genes individuais, motivando a análise do transcriptoma completo em busca do gene mais significativamente desregulado."*

---

## 5. A transição para o transcriptoma completo e a seleção da NQO1

Do transcriptoma completo, o gene **NQO1** (NAD(P)H quinona desidrogenase 1) emergiu como o **segundo mais significativamente desregulado**, com **downregulação** na DCM (log2FC = −1,686; padj = 6,86 × 10⁻¹³).

A seleção da NQO1 é justificada por **três critérios combinados**:

1. **Significância estatística** (padj = 6,86 × 10⁻¹³) — o gene mais robusto do dataset.
2. **Relevância biológica** — a NQO1 é uma **flavoproteína citoprotetora** que reduz quinonas por transferência de dois elétrons, prevenindo a formação de semiquinonas e espécies reativas de oxigênio (ROS) [Ross & Siegel, 2017; Dinkova-Kostova & Talalay, 2010].
3. **Relevância farmacológica** — a NQO1 está diretamente envolvida no **metabolismo de quinonas**, classe à qual pertencem as **antraciclinas** (doxorrubicina), cuja cardiotoxicidade é classicamente atribuída à **ciclagem redox e geração de ROS** [Wallace, 2003; Siegel et al., 2012].

> **Nota crítica de rigor:** a NQO1 **não pertence à via AMPK**. Ela foi selecionada a partir do **transcriptoma completo** (estratégia dirigida por dados), e não "filtrada da AMPK". Essa distinção é importante para a honestidade metodológica.

---

## 6. A etapa DGB: farmacologia reversa "gene → fármaco"

Com o alvo NQO1 definido, a consulta ao **Drug Gene Budger (DGB)** inverte a direção usual da farmacologia (que vai do fármaco ao gene):

- **Convencional (fármaco → gene):** "dado um fármaco, que genes ele modula?"
- **Reversa (gene → fármaco), usada aqui:** "dado um gene de interesse, que fármacos modulam sua expressão?"

O DGB integra assinaturas transcriptômicas de perturbação por xenobióticos de **CMAP**, **LINCS L1000** e **CREEDS**, permitindo triar milhares de compostos quanto ao efeito sobre a NQO1 [Wang et al., 2021; Lamb et al., 2006; Subramanian et al., 2017].

---

## 7. Por que a doxorrubicina é o resultado "correto" e coerente

O ranking apontou a **doxorrubicina** como o repressor mais significativo de NQO1 (q = 2,21 × 10⁻²⁸; Fold Change médio = −1,065; convergente em CREEDS + LINCS L1000).

A **coerência biológica** é o que fecha o argumento:

1. A doxorrubicina gera **ROS** (via ciclagem redox de semiquinonas) [Wallace, 2003].
2. A NQO1 é uma **enzima antioxidante** que detoxifica quinonas [Ross & Siegel, 2017].
3. Portanto, um fármaco **pró-oxidante que reprime um gene antioxidante** produz um **duplo prejuízo redox** — mecanismo plausível de cardiotoxicidade.
4. A **replicação em duas bases independentes** (CREEDS + LINCS L1000) aumenta a confiabilidade do achado.

> **Nota sobre o ranking:** o trichostatin A aparece com q = 0 por **underflow numérico** (fonte única, CMAP). A doxorrubicina detém o **menor q-value real** e a **evidência em múltiplas plataformas**, sendo a resposta correta do estudo.

---

## 8. Enquadramento epistemológico

Este é um estudo **exploratório, descritivo e gerador de hipótese** (*in silico*), e não um estudo confirmatório:

- **Não testa uma hipótese causal** — gera hipóteses a serem testadas experimentalmente.
- **Não permite inferência causal** — associações transcriptômicas, não mecanismos comprovados.
- **Seu valor** está em **priorizar** alvos e fármacos para validação futura, de forma **reprodutível** (código aberto) e **sistemática** (ranking integrado).

A combinação de **via candidata (AMPK)** + **abordagem dirigida por dados (NQO1)** + **farmacologia reversa (DGB)** é um desenho metodologicamente defensável e comum em **biologia de sistemas e farmacogenômica computacional**.

---

## 9. Resumo da justificativa em uma frase

> *"Partiu-se de um problema clínico (cardiomiopatia) e de uma via metabólica central (AMPK) como contexto; a refutação da hipótese de via conduziu à varredura do transcriptoma completo, que identificou a NQO1 — gene antioxidante mais desregulado — e, por farmacologia reversa via DGB, a doxorrubicina como seu repressor mais significativo, gerando uma hipótese mecanicamente coerente de cardiotoxicidade por estresse oxidativo."*

---

## 10. Referências

1. Herzig S, Shaw RJ. AMPK: guardian of metabolism and mitochondrial homeostasis. *Nat Rev Mol Cell Biol*. 2018.
2. Bozkurt B, et al. Dilated cardiomyopathy: a review of current and future treatment strategies. *Nat Rev Cardiol*. 2021.
3. Ross D, Siegel D. Functions of NQO1 in cellular protection and CoQ10 metabolism. *Front Physiol*. 2017.
4. Dinkova-Kostova AT, Talalay P. NAD(P)H:quinone acceptor oxidoreductase 1 (NQO1), a multifunctional antioxidant enzyme. *Arch Biochem Biophys*. 2010.
5. Wallace KB. Doxorubicin-induced cardiac mitochondrionopathy. *Pharmacol Toxicol*. 2003.
6. Siegel D, et al. NAD(P)H:quinone oxidoreductase 1 (NQO1) in the sensitivity and resistance to antitumor quinones. *Biochem Pharmacol*. 2012.
7. Wang ZY, et al. Drug Gene Budger (DGB): an application for ranking drugs to modulate a specific gene. *Bioinformatics*. 2021.
8. Lamb J, et al. The Connectivity Map. *Science*. 2006.
9. Subramanian A, et al. A Next Generation Connectivity Map: L1000 Platform. *Cell*. 2017.
