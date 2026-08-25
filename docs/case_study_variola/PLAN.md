# Planejamento de escrita — estudo de caso Orthopoxvirus / PhyloTreeMiner

**Status:** rascunho de trabalho, 12 de agosto de 2026
**Texto acadêmico associado:** [`PAPER.md`](PAPER.md)
**Código que produz todos os números citados:** `workflow/stability/` (executar `python -m workflow.stability.case_study`)

---

## 1. Tese

O artigo defende **duas** afirmações, uma metodológica e uma biológica, ligadas por um único mecanismo.

> **Metodológica.** O suporte de um clado *entre pipelines* — a fração de combinações alinhador × método de inferência que o recuperam — é uma medida de robustez ortogonal ao bootstrap, e mede a única fonte de erro que o bootstrap é estruturalmente incapaz de ver: a escolha do método. Medi-la corretamente exige uma **identidade de clado canônica**; sem ela, a medida colapsa.
>
> **Biológica.** Em um painel de 121 genomas completos de Orthopoxvirus, essa medida se comporta como um detector de resolução filogenética: o sinal que sobrevive a todos os pipelines é exatamente a estrutura de nível de espécie, enquanto o *backbone* profundo entre espécies é inteiramente dependente do método. O alinhamento genômico completo, na prática, não sustenta a topologia profunda que rotineiramente se publica a partir dele.

O elo entre as duas: a segunda afirmação só é mensurável porque a primeira foi corrigida. Esse é o argumento de coesão do artigo — não é "ferramenta + aplicação", é "a medida estava quebrada, consertamos, e o conserto revelou onde o sinal genômico realmente termina".

---

## 2. O que já está estabelecido

Todos os números abaixo saem de `projects/*/out/outputs/stability/summary.json`, gerados a partir das árvores já existentes no repositório. Nada aqui depende de rodar novos experimentos.

### 2.1 Três modos independentes de falha de identidade

O esquema anterior (`treeUtils.encode_list_to_int`) identifica um clado pelo MD5 truncado em 16 bits da *representação textual da lista ordenada* dos hashes de seus terminais. Isso falha de três maneiras, e as três foram quantificadas:

| Falha | Mecanismo | Efeito na medida | VARV‑121 |
|---|---|---|---|
| Dependência de ordem | o hash muda se os filhos forem percorridos em outra ordem | **subestima** suporte (fragmenta um clado em vários itens) | 95 de 269 clados (35,3%) fragmentados |
| Colisão de 16 bits | 65 536 valores para milhares de clados | **superestima** suporte (funde clados distintos) | 1 item em colisão (0,2%) |
| Adulteração de rótulo | IQ‑TREE e RAxML‑NG truncam o dígito de versão de acessos RefSeq (`NC_008030.1` → `NC_008030.`) | **subestima** suporte; táxons afetados deixam de casar entre pipelines | 4 táxons, 2 pipelines |

Efeito combinado, em uma frase que serve de gancho para o abstract: **o esquema anterior reporta zero clados recuperados por todos os 8 pipelines; a identidade canônica encontra 31.** Em k ≥ 6 a razão é 50 contra 5.

### 2.2 O aligner não importa; o método de inferência é tudo

Distância Robinson–Foulds normalizada média entre pares de pipelines que diferem em **um único fator**:

| Experimento | trocando o alinhador | trocando a inferência | razão |
|---|---|---|---|
| VARV‑121 (121 táxons) | 0,004 (máx 0,017) | 0,548 (máx 0,676) | ~130× |
| VARV‑49 (49 táxons) | 0,000 | 0,488 (máx 0,606) | ∞ |
| VARV‑6‑noITR (6 táxons) | 0,000 | 0,625 (máx 0,875) | ∞ |

MAFFT e Clustal Omega produzem topologias **idênticas** (RF = 0) para todo método de inferência, nas três escalas; a única exceção é IQ‑TREE em 121 táxons, com RF = 0,017. Este é o resultado mais limpo do conjunto e o mais fácil de comunicar.

### 2.3 A estabilidade rastreia a taxonomia

Dos 14 clados com ≥ 3 táxons recuperados por todos os 8 pipelines, **13 são de espécie única**. A única exceção é `{VACV, VACV, BPXV}` — buffalopox aninhado em vaccinia, que é a biologia correta, não um erro. Os clados universais incluem: MPXV (23/23 táxons), CMLV (11/11), e um clado VARV de 8.

No sentido oposto: em suporte 4/8 aparecem clados de 116, 65, 64, 62 e 50 táxons com pureza taxonômica de 0,35 a 0,62 — o backbone profundo é grande, method-dependent e taxonomicamente incoerente.

### 2.4 A mineração de padrões maximais recupera as famílias algorítmicas sem ser informada delas

Enumerando exatamente o reticulado de conjuntos fechados (barato: 2^M com M = número de pipelines, não de táxons), os blocos de 4 pipelines que mais compartilham clados são:

- `{FastTree, IQ-TREE} × {MAFFT, Clustal}` → 95 clados em comum (família ML)
- `{NJ, UPGMA} × {MAFFT, Clustal}` → 57 clados (família de distância)
- qualquer par cruzando as famílias → 44 clados

O maior padrão maximal frequente tem 38 clados suportados por 5 dos 8 pipelines. Isto é uma validação atraente do método: a estrutura de agrupamento dos algoritmos emerge dos dados.

---

## 3. O que falta para o alvo

Avaliação honesta: **o pilar metodológico está pronto; o pilar biológico ainda é confirmatório.** Recuperar monofilia de MPXV e CMLV valida o método, mas não é achado novo. Para que a afirmação dupla se sustente em veículo de alto impacto, é preciso converter "o backbone profundo é instável" de observação em resultado. Em ordem de prioridade:

| # | Experimento | Por que é necessário | Esforço |
|---|---|---|---|
| **A** | **Bootstrap × estabilidade entre pipelines no mesmo dataset.** Para cada clado, o par (suporte bootstrap dentro do pipeline, suporte entre pipelines). | É *a* comparação que a tese exige. Sem ela, "ortogonal ao bootstrap" é asserção. A previsão testável e forte: existem clados com bootstrap ≥ 95 e suporte entre pipelines ≤ 0,5 — alta confiança em um artefato. | médio (IQ‑TREE com `-B 1000` já está no ambiente) |
| **B** | **Ampliar o delineamento de alinhadores.** MUSCLE, PRANK e, se viável, Cactus/progressiveMauve para genomas completos. | O achado "o alinhador não importa" é atualmente baseado em dois alinhadores progressivos que compartilham heurísticas. Um revisor vai apontar isso na primeira leitura. Com PRANK (consciente de filogenia) e um alinhador genômico real, o achado ou fica muito mais forte ou é corretamente qualificado. | médio‑alto (custo de CPU no alinhamento de 121 genomas) |
| **C** | **Mascaramento de ITRs no painel completo de 121.** Hoje o experimento sem ITRs tem só 6 táxons e uma composição diferente — não é controle do principal. | As ITRs são repetições terminais invertidas de dezenas de kb; sua presença em um alinhamento genômico completo é uma explicação candidata para a instabilidade do backbone. Testar isso transforma a observação em mecanismo. | médio |
| **D** | **Validação em dados simulados.** Simular sob topologia conhecida (INDELible/AliSim), aplicar os 8 pipelines, verificar se estabilidade entre pipelines prediz clados corretos melhor que bootstrap. | Converte a medida de descritiva em validada. É o que separa "métrica interessante" de "métrica que você deveria usar". | alto |
| **E** | **Um segundo táxon-alvo.** Ex.: o dataset de Zika já presente no repositório (`projects/Zika_*`, até 480 sequências). | Generalidade. Um revisor perguntará se o efeito é peculiar a genomas de poxvírus, que são grandes, ricos em AT e com ITRs. | baixo (dados já existem) |

**Recomendação.** A + E são baratos e, juntos, já sustentam uma submissão forte. B + C elevam o teto. D é o que seria necessário para um veículo da família *Nature*; sem D, o enquadramento honesto é "medida nova, bem validada empiricamente", não "medida provadamente superior".

---

## 4. Estratégia de veículo

| Veículo | Encaixe | Avaliação |
|---|---|---|
| *Nature Methods* | alto, **se** D for feito | O formato deles é exatamente "métrica nova + validação + caso real". Sem validação em simulação, desk reject é provável. |
| *Nature Communications* | médio‑alto | Aceita a combinação metodológico + aplicado com menos exigência de validação formal. Alvo realista com A + B + C + E. |
| *Molecular Biology and Evolution* | alto | Público que se importa com a diferença entre robustez amostral e metodológica. Provavelmente o melhor encaixe real. |
| *Bioinformatics* (Original Paper) | alto, teto mais baixo | Aceitaria essencialmente o que já existe hoje, com A e E. |
| *Systematic Biology* | médio | Vai exigir tratamento teórico da medida (propriedades, consistência), não apenas empírico. |

**Recomendação:** escrever para *Nature Communications* / *MBE* como alvo primário, mantendo o texto compatível com o formato Article da família *Nature* (é o que `PAPER.md` faz). Rebaixar para *Bioinformatics* é reformatação, não reescrita.

---

## 5. Plano de figuras

| Fig | Conteúdo | Estado |
|---|---|---|
| **1** | Esquema: delineamento fatorial (2 alinhadores × 4 métodos), transação = pipeline, item = clado; ao lado, os três modos de falha de identidade com um exemplo mínimo de cada | **a desenhar** |
| **2** | Perfil de suporte: número de clados recuperados por ≥ k pipelines, identidade canônica vs. legada | **pronta** — `fig_support_profile.png` |
| **3** | Matriz RF entre os 8 pipelines, com a estrutura de blocos ML / distância visível | **pronta** — `fig_rf_heatmap.png` |
| **4** | Tamanho do clado × suporte entre pipelines, colorido por coerência taxonômica; mostra a separação entre a copa estável e o backbone instável | **pronta** — `fig_size_vs_support.png` |
| **5** | Bootstrap × suporte entre pipelines (dispersão), com o quadrante "alto bootstrap, baixa estabilidade" destacado | **bloqueada pelo experimento A** |
| **ED 1** | Réplicas em VARV‑49 e VARV‑6‑noITR (as três tabelas de fatores) | **pronta** |
| **ED 2** | Padrões maximais e os blocos algorítmicos que eles recuperam | dados prontos, figura a desenhar |
| **ED 3** | Efeito do mascaramento de ITRs | **bloqueada pelo experimento C** |

Sistema visual já fixado em `workflow/stability/report.py`: paleta categórica validada para deficiência de visão de cores (azul `#2a78d6` / laranja `#eb6834`), rampa sequencial de matiz único, eixos recessivos, rótulos diretos.

---

## 6. Estrutura do texto e alocação

Formato Article da família *Nature*, ~3 000 palavras no corpo.

| Seção | Palavras | Conteúdo |
|---|---|---|
| Abstract | 150 | Os dois achados, com os números 0→31 e 130×. |
| Introduction | 500 | Bootstrap mede reamostragem, não escolha de método; a lacuna; por que a identidade de clado é o obstáculo técnico. |
| Results §1 | 400 | Identidade canônica e os três modos de falha, quantificados. |
| Results §2 | 400 | Alinhador vs. inferência; RF = 0 entre alinhadores. |
| Results §3 | 500 | O perfil de estabilidade e sua correspondência com a taxonomia. |
| Results §4 | 400 | Padrões maximais recuperam famílias algorítmicas. |
| Results §5 | 350 | Bootstrap vs. estabilidade *(pendente do experimento A)*. |
| Discussion | 500 | O que isso significa para filogenômica de genoma completo; limitações; o que a medida não é. |
| Methods | sem limite | Dataset, pipelines, identidade canônica, enumeração exata do reticulado, disponibilidade. |

---

## 7. Objeções de revisor previstas

| Objeção | Resposta |
|---|---|
| "UPGMA e NJ são métodos ruins; incluí-los infla artificialmente a discordância." | Justo, e é por isso que o resultado é reportado **por par de métodos**, não como um número agregado. O ponto se sustenta mesmo dentro da família ML: FastTree × IQ‑TREE dá RF = 0,18–0,19, ordens de grandeza acima do efeito do alinhador. Reportar isso explicitamente no corpo do texto, não em suplementar. |
| "Dois alinhadores progressivos não são amostra de espaço de alinhamento." | Reconhecer como limitação e executar o experimento B. Sem B, restringir a afirmação a "entre alinhadores progressivos de uso corrente". |
| "Isto é apenas um consenso de métodos com nome novo." | Não: o consenso descarta as árvores discordantes; a mineração de padrões maximais preserva *quais* pipelines concordam em *quais* clados, e é isso que recupera a estrutura de famílias algorítmicas (§2.4). O consenso não tem esse produto. |
| "A monofilia de MPXV/CMLV já é conhecida." | Correto — é validação, não descoberta, e o texto deve dizer isso literalmente. O achado biológico é o **negativo**: o backbone profundo não é sustentado. |
| "Crocodilepox é um outgroup distante demais e distorce a raiz." | Verdadeiro e precisa ser tratado (ver §8). Rodar novamente com outgroup adequado (Taterapox/Camelpox) e reportar ambas as versões. |
| "16 bits era obviamente insuficiente; corrigir isso não é resultado." | O resultado não é o conserto, é a magnitude do viés que ele revela: 0 versus 31. Esquemas de hash truncado são comuns em pipelines de bioinformática; a lição é transferível. |

---

## 8. Problemas do dataset a resolver antes de submeter

Levantados na inspeção; nenhum invalida a análise metodológica, mas todos são munição de revisor.

1. **O nome não descreve o conteúdo.** O diretório é `..._li_et_al_2007_replication`, mas o painel de 121 táxons não é uma réplica de Li et al. 2007: são 77 VARV, 23 MPXV, 11 CMLV, 3 crocodilepox, 2 CPXV, 2 VACV, 1 BPXV, 1 YOKA, 1 TATV. Renomear e descrever honestamente como painel Orthopoxvirus.
2. **Acesso duplicado.** `NC_008291.1` (Taterapox) aparece duas vezes no FASTA: 125 registros, 124 acessos únicos.
3. **Três táxons desaparecem entre o FASTA e as árvores** (`DQ437594.1`, `HQ849551.1`, `NC_003391.1`): 124 → 121. Descobrir onde e por quê (provável remoção de sequências idênticas), e documentar a etapa.
4. **Crocodilepox como outgroup.** Três genomas de crocodilepox pertencem a uma subfamília distinta de Poxviridae. Como outgroup de Orthopoxvirus é distante demais e é candidato natural a explicar a instabilidade do backbone. Refazer com outgroup próximo.
5. **`Variola_Yu_li_2007_noITRs_6seqs` não é controle do experimento principal:** táxons e escala diferentes. Substituir pelo experimento C.

---

## 9. Sequência de trabalho sugerida

1. Corrigir a proveniência do dataset (§8, itens 1–4) — sem isso os números finais mudarão depois de escritos.
2. Rodar o experimento A (bootstrap). É a evidência de maior retorno por unidade de esforço.
3. Rodar o experimento E (Zika, dados já no repositório) para generalidade.
4. Escrever o primeiro rascunho completo com Figuras 1–5 (`PAPER.md` já contém tudo que não depende de A).
5. Rodar B e C em paralelo à revisão interna.
6. Decidir o veículo com base em se D será feito.

---

## 10. Reprodutibilidade

- Toda a análise: `python -m workflow.stability.case_study` (a partir de `BioComp_UFF/`), sem dependências além de Biopython, NumPy, pandas e Matplotlib — **`mlxtend` não é mais necessário**, porque a enumeração do reticulado é exata.
- Testes: `python -m unittest workflow.tests.test_stability` (16 testes; cobrem invariância à ordem, reconciliação de rótulos, maximalidade dos padrões e separação de fatores).
- Saídas por experimento em `projects/<projeto>/out/outputs/stability/`: `clade_support.csv`, `maximal_patterns.csv`, `rf_matrix.csv`, `summary.json` e as três figuras.
- Para a submissão: fixar a imagem Docker por digest, depositar as árvores e alinhamentos no Zenodo e citar o DOI em Data Availability.
