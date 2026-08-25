# Cross-pipeline clade stability separates resolved from unresolved structure in whole-genome orthopoxvirus phylogenies

**Draft manuscript — Article format.** Every number in this text is produced by
`python -m workflow.stability.case_study` from the trees deposited in this
repository. Claims that depend on experiments not yet run are marked
**[pending]** and must be removed or completed before submission.

---

## Abstract

Phylogenomic conclusions are routinely reported with bootstrap support, a measure
of how a single analytical pipeline responds to resampling its own alignment
columns. It is silent about a larger source of variation: the choice of pipeline
itself. Here we measure clade support *across* pipelines — the fraction of
aligner × inference-method combinations that recover a given clade — in a panel
of 121 complete orthopoxvirus genomes. Measuring this quantity at all requires a
clade identity that is invariant to traversal order and free of hash collisions;
we show that the ordered, 16-bit identity in common use fragments 35.3% of clades
and reports **zero** clades shared by all eight pipelines where a canonical
identity finds 31. With identity fixed, two results follow. First, topological
variation is governed almost entirely by the inference method and essentially not
at all by the aligner: swapping MAFFT for Clustal Omega leaves the topology
unchanged (mean normalised Robinson–Foulds distance 0.004) while swapping the
inference method changes it profoundly (0.548), a ratio of roughly 130. Second,
cross-pipeline support behaves as a detector of phylogenetic resolution: 13 of
the 14 clades of three or more taxa recovered by every pipeline are
single-species groups, whereas the deep inter-species backbone — clades of 50 to
116 taxa — is recovered by at most half of the pipelines and is taxonomically
incoherent. Whole-genome orthopoxvirus alignments support species-level
structure robustly and the deep backbone not at all, a distinction invisible to
within-pipeline support.

---

## Main

Bootstrap resampling asks a precise question: if the alignment columns had been
drawn differently, would this clade survive? It is a question about sampling
error, posed *inside* a fixed analytical pipeline — one aligner, one substitution
model, one tree-search algorithm. It cannot, by construction, report on the
consequences of having chosen that pipeline rather than another. Yet in
whole-genome phylogenomics the pipeline is often the dominant free parameter, and
different reasonable choices are known to yield conflicting topologies.

The natural complement is to hold the data fixed and vary the method: build the
tree under every combination of aligner and inference algorithm and ask, for each
clade, how many combinations recover it. We call this quantity **cross-pipeline
support**. It is cheap to compute — the pipelines are usually run anyway during
method selection — and it measures exactly the error term the bootstrap omits.

Computing it, however, exposes a technical obstacle that turns out to dominate
the measurement. Cross-pipeline support requires deciding when a clade in one
tree *is the same clade* as one in another. That decision is trivial in principle
(two clades are the same when they contain the same taxa) and is routinely
implemented in a way that is not. Below we quantify the damage, correct it, and
then use the corrected measure on a panel of 121 orthopoxvirus genomes.

### Clade identity determines what can be measured

The identity scheme we audit — used by the frequent-subtree miner in the
PhyloTreeMiner workflow, and representative of a common pattern — hashes the
textual representation of the *ordered list* of terminal hashes and truncates the
digest to 16 bits. Three failure modes follow, and they push the measurement in
opposite directions.

**Order dependence.** The same clade, traversed with its children in a different
order, hashes differently. Two trees that agree perfectly on a clade are recorded
as disagreeing. In the 121-genome panel, 95 of 269 distinct clades (35.3%) were
split across more than one identifier; the effect is comparable at the other
scales tested (40.0% at 49 taxa, 40.0% at 6 taxa). This *deflates* measured
agreement.

**Hash collision.** Sixteen bits provide 65,536 identifiers. With hundreds of
clades observed across eight trees, collisions are expected rather than rare: one
identifier in the 121-genome panel was shared by two genuinely distinct clades,
merging them into a single apparently better-supported item. This *inflates*
measured agreement.

**Label mangling by upstream tools.** IQ-TREE rewrote four RefSeq accessions when
emitting Nexus, truncating the version digit (`NC_008030.1` → `NC_008030.`);
RAxML-NG did the same in the six-taxon experiment. Under any exact-match
identity, every clade containing an affected taxon fails to match its counterpart
in the other pipelines. Left unreconciled, this alone inflated the apparent
number of distinct clades from 269 to 284 and concealed one of the 31 universal
clades.

The combined effect is not marginal. Under the legacy identity, **no clade
whatsoever** appears in all eight pipelines, and only five appear in at least
six. Under a canonical identity — the unordered set of terminal names, digested
at 128 bits, with accession versions reconciled — the same trees yield 31 and 50
respectively (Fig. 2). An analysis using the legacy scheme would have concluded
that these eight pipelines share no common phylogenetic signal at all. They share
a great deal.

We therefore treat canonical clade identity not as an implementation detail but
as a precondition: cross-pipeline support is unmeasurable without it.

### The aligner does not matter; the inference method is everything

The panel comprises 121 complete orthopoxvirus genomes (77 variola, 23 monkeypox,
11 camelpox, 3 crocodilepox, 2 cowpox, 2 vaccinia, and one each of buffalopox,
Yoka and taterapox), aligned to 283,874 columns. Eight pipelines were run: MAFFT
and Clustal Omega crossed with neighbour-joining, UPGMA, FastTree and IQ-TREE.
The factorial design lets each factor's contribution be isolated by comparing
pipeline pairs that differ in one factor only.

The separation is stark (Fig. 3). Across pairs differing only in the aligner, the
mean normalised Robinson–Foulds distance is **0.004**, with a maximum of 0.017;
MAFFT and Clustal Omega produced *topologically identical* trees for
neighbour-joining, UPGMA and FastTree, and differed for IQ-TREE by four
bipartitions out of 238. Across pairs differing only in the inference method, the
mean is **0.548**, with a maximum of 0.676 — a ratio of approximately 130.

This is not an artefact of including distance methods. Within the maximum-
likelihood family, FastTree and IQ-TREE differ by 0.18–0.19 — still one to two
orders of magnitude above the aligner effect. The pattern replicates at both
other scales tested: 0.000 versus 0.488 in a 49-genome variola panel, and 0.000
versus 0.625 in a six-genome panel aligned after inverted terminal repeat
removal.

The practical implication is direct. Effort spent comparing progressive aligners
for this class of data buys almost nothing in topological terms, while the choice
of inference method silently determines the result. We note the scope: both
aligners tested are progressive and share heuristics, and the claim should be
read as bounded by that.

### Cross-pipeline support tracks taxonomic resolution

Of 269 distinct clades observed, 31 (11.5%) were recovered by all eight
pipelines, 50 by at least six, and 123 by at least four. Restricting to clades of
three or more taxa — those carrying real topological content — 14 are universal,
and **13 of these are single-species groups** (Fig. 4). Monkeypox virus is
recovered as monophyletic by every pipeline (all 23 isolates), as is camelpox
(all 11), alongside a stable eight-genome variola clade and nested subclades of
each.

The single exception is instructive rather than troubling: the clade
`{vaccinia, vaccinia, buffalopox}`. Buffalopox is a vaccinia-derived lineage, so
the one universal clade that is not species-pure is the one where species labels
and evolutionary history genuinely disagree. The measure recovered the biology,
not an error.

The complementary observation carries the weight of the paper. At support 4/8,
clades of 116, 65, 64, 62 and 50 taxa appear, with taxonomic purity between 0.35
and 0.62 — large, deep groupings that mix every species in the panel and that
half the pipelines reject. The deep inter-species backbone of this whole-genome
alignment is not merely weakly supported; it is *method-determined*. Which
species groups with which depends on whether one ran UPGMA or IQ-TREE.

Cross-pipeline support thus partitions the tree into a stable canopy and an
unstable backbone, and does so without external information. The boundary it
draws — at roughly the species level — is a statement about how deep the
phylogenetic signal in these alignments actually reaches.

### Maximal pattern mining recovers the algorithmic family structure

Because the number of pipelines M is small, the lattice of closed clade sets can
be enumerated exactly in O(2^M·|C|) time; no heuristic frequent-itemset mining is
required, and the maximal patterns returned are exact rather than approximate.

Applied to the panel, the four-pipeline groups sharing the most clades are
{FastTree, IQ-TREE} × {MAFFT, Clustal} with 95 clades in common, and
{NJ, UPGMA} × {MAFFT, Clustal} with 57; any group spanning the two families
shares 44. The largest maximal frequent pattern comprises 38 clades supported by
five of the eight pipelines.

The mining therefore reconstructs the division between maximum-likelihood and
distance methods from topological agreement alone, having been given no
information about the algorithms. This is a property a majority-rule consensus
tree cannot provide: a consensus discards the identity of the disagreeing trees,
whereas the pattern lattice preserves which pipelines agree on which structure.

---

## Discussion

Three conclusions follow, of decreasing generality.

**Identity is a measurement instrument, not plumbing.** Truncated, order-sensitive
hashes are widespread in bioinformatics pipelines, where they are adopted for
compactness and treated as neutral. Here such a scheme did not add noise — it
inverted the result, reporting zero cross-pipeline agreement where substantial
agreement exists. Any pipeline that identifies biological objects by hashing an
ordered representation should be audited before its outputs are interpreted; the
audit is inexpensive and, in our case, decisive.

**For this class of data, pipeline variance lives in the tree search.** The
~130-fold asymmetry between the inference and alignment factors argues that
method-comparison effort in whole-genome viral phylogenomics is currently
misallocated. We emphasise the boundaries of this claim: two progressive
aligners, one viral family, genomes with unusual composition and long terminal
repeats. Extending the aligner panel — in particular to a phylogeny-aware
aligner and to a true genome aligner — is the first experiment that could
overturn it.

**Whole-genome alignment does not license deep topology here.** The stable
structure stops at species boundaries. This does not mean the deep relationships
are unknowable; it means these alignments, analysed this way, do not determine
them, and that a bootstrap value computed inside any one pipeline will not
reveal this. We expect — but have not yet shown — that clades exist with high
bootstrap support and low cross-pipeline support, and that these are precisely
the ones where confidence is misplaced. **[pending: experiment A]**

Cross-pipeline support is not a replacement for the bootstrap; the two measure
different errors and can be reported together at no additional cost, since the
pipelines are typically run regardless. Nor is it a consensus method: it retains
the attribution of agreement to specific pipelines, which is what enables the
algorithmic structure of the analysis to be read off the results.

**Limitations.** The measure is defined relative to the panel of pipelines
chosen, and a panel weighted toward one algorithmic family will inflate apparent
stability; reporting the design explicitly, as we do, is therefore mandatory. The
biological findings here are confirmatory — monkeypox and camelpox monophyly are
established — and serve as validation that the measure recovers known structure;
the novel biological content is the negative result about the backbone. The
panel includes three crocodilepox genomes as an outgroup, which is distant enough
that its influence on the deep topology should be tested directly before the
backbone instability is attributed to the data rather than to the rooting.
Finally, validation against a known simulated topology remains to be done, and
until it is, cross-pipeline support should be read as a well-characterised
descriptive measure rather than a demonstrated improvement over existing support
values. **[pending: experiments B–E]**

---

## Methods

### Sequence data

Complete genomes were retrieved from GenBank and assembled into a panel of 121
orthopoxvirus genomes (`data/workflow_dataAcquisition_li_et_al_2007_replication-RetMax200/dataset_final.fasta`).
Two replicate panels were analysed: a 49-genome variola panel and a six-genome
panel from which inverted terminal repeats had been removed. Accession lists,
including the three records present in the FASTA but absent from the trees, are
given in Supplementary Table 1.

### Alignment and tree inference

Each panel was aligned with MAFFT and with Clustal Omega under the workflow's
default parameters, yielding alignments of 283,874, 235,955 and 250,517 columns
respectively. From each alignment, trees were inferred by neighbour-joining,
UPGMA, FastTree and IQ-TREE (RAxML-NG additionally for the six-genome panel),
giving eight pipelines per panel (ten for the six-genome panel). All software
versions are pinned by the workflow's Docker image.

### Canonical clade identity

A clade is identified by the unordered set of names of its terminal descendants.
For storage and comparison this set is reduced to a 128-bit MD5 digest computed
over the lexicographically sorted, newline-joined names, making the identifier
invariant to traversal order, to node rotation, and to the arbitrary rooting of
subtrees. Leaves and the universal clade are excluded as uninformative.

Terminal labels are reconciled before identification by stripping the GenBank
accession version suffix, which resolves the truncation introduced by IQ-TREE and
RAxML-NG. A label is recorded as mangled when its raw spelling differs from the
modal spelling across pipelines for the same taxon, with ties broken toward the
longer spelling, since mangling here is truncation. Reconciliation is reported
per pipeline so that its effect on any result can be audited.

The legacy identity, reimplemented for the audit, hashes `str(list_of_terminal_hashes)`
with MD5 and retains the first four hexadecimal digits, where each terminal hash
is itself the first four hexadecimal digits of the MD5 of the terminal name.

### Cross-pipeline support and pattern mining

Each pipeline contributes one transaction whose items are its canonical clades.
The support of a clade is the fraction of pipelines containing it. For each
non-empty subset S of the M pipelines, the intersection I(S) of their clade sets
is closed; its carrier set is the set of all pipelines containing I(S).
Enumerating all 2^M − 1 subsets yields every closed pattern exactly, and maximal
frequent patterns are obtained by discarding closed patterns properly contained
in another pattern above the support threshold. With M ≤ 10, this is exhaustive
and takes well under a second, so no heuristic miner is used.

### Topological distance and factor effects

Robinson–Foulds distance is the size of the symmetric difference between the two
clade sets, normalised by 2(n − 2), the maximum for rooted binary trees on n taxa
with the universal clade excluded. Factor effects are computed by partitioning
the pipeline pairs into those differing only in the aligner, only in the
inference method, and in both, and reporting the mean, minimum and maximum
normalised distance within each partition.

### Taxonomic coherence

Species labels were assigned from the FASTA description of each accession by
keyword match. A clade's purity is the fraction of its taxa carrying the majority
label; a clade is single-species when purity is 1.

### Code and data availability

The analysis is implemented in `workflow/stability/` and reproduced end to end by
`python -m workflow.stability.case_study`. It depends only on Biopython, NumPy,
pandas and Matplotlib. A unit-test suite covering order invariance, label
reconciliation, pattern maximality and factor separation runs with
`python -m unittest workflow.tests.test_stability`. Per-experiment outputs —
`clade_support.csv`, `maximal_patterns.csv`, `rf_matrix.csv`, `summary.json` and
all figures — are written to `projects/<project>/out/outputs/stability/`.

---

## Figures

**Fig. 1 | Design and the three identity failure modes.** *(to be drawn)*
The factorial design: two aligners crossed with four inference methods; each
pipeline is one transaction, each clade one item. Insets give a minimal example
of order dependence, hash collision and label mangling.

**Fig. 2 | Cross-pipeline clade support under canonical and legacy identity.**
`fig_support_profile.png` — number of clades recovered by at least k of eight
pipelines. The legacy identity reports zero clades at k = 8 and five at k ≥ 6,
against 31 and 50 under canonical identity.

**Fig. 3 | Topological distance between pipelines.**
`fig_rf_heatmap.png` — normalised Robinson–Foulds distances among the eight
pipelines. Aligner-only swaps lie on the zero diagonal blocks; the
maximum-likelihood and distance families separate cleanly.

**Fig. 4 | Taxonomic coherence versus methodological support.**
`fig_size_vs_support.png` — clade size against cross-pipeline support,
distinguishing single-species from mixed-species clades. Universal clades are
small and species-pure; large clades are supported by at most half the pipelines
and are taxonomically mixed.

**Fig. 5 | Bootstrap support versus cross-pipeline support.** **[pending:
experiment A]**

**Extended Data Fig. 1 | Replication at 49 and 6 genomes.**
**Extended Data Fig. 2 | Maximal patterns and the algorithmic family blocks.**
**Extended Data Fig. 3 | Effect of inverted terminal repeat masking.** **[pending:
experiment C]**
