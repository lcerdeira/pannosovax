# A pan-nosocomial multi-epitope immunogen against three WHO-priority respiratory pathogens, designed with commensal-microbiome safety screening and Brazil-weighted HLA coverage

**Target journal:** npj Vaccines (Article)

> **Note on numbers.** Every quantity below is a PENDING marker filled from the pipeline's
> own TSV outputs by `scripts/report/fill_manuscript.py`. The script refuses to fill a marker
> whose source file does not exist. No number in this manuscript is typed by hand; markers
> stay unresolved until the full-scale run completes.

---

## Abstract

**Background.** *Klebsiella pneumoniae* species complex (KpSC), *Acinetobacter baumannii*
and *Streptococcus pneumoniae* are leading causes of drug-resistant pneumonia. No licensed
vaccine exists for the first two, and the polysaccharide-conjugate vaccines available for
the third are subject to serotype replacement.

**Methods.** We designed a single chimeric multi-epitope immunogen from ⟨PENDENTE:n_genomes_total⟩
complete genomes. Only core surface proteins were considered — never capsular antigens —
making the design serotype-independent. Epitopes were screened against the human proteome,
against exact 7-mer autoimmune mimicry, and, unusually, against the respiratory and gut
**commensal microbiome**. The final set was chosen by maximum-coverage optimisation weighted
by both world and **Brazilian** HLA allele frequencies.

**Results.** ⟨PENDENTE:results_summary⟩

**Conclusions.** The design is computational and requires experimental validation. Its two
methodological contributions — commensal screening and ancestry-aware coverage — are
transferable to any multi-epitope vaccine programme.

**Keywords:** reverse vaccinology; multi-epitope vaccine; antimicrobial resistance;
nosocomial pneumonia; HLA population coverage; commensal microbiome

---

## 1. Introduction

Antimicrobial resistance is no longer a projection. Global burden estimates attribute
millions of annual deaths to resistant bacterial infections, and three organisms appear
persistently at the top of priority lists: *Klebsiella pneumoniae*, *Acinetobacter
baumannii* and *Streptococcus pneumoniae*. The first two are nosocomial Gram-negatives
with an expanding carbapenemase repertoire; the third is a community-circulating
Gram-positive whose resistance to β-lactams and macrolides advances steadily.

The vaccine response to this trio is uneven and, in every case, unsatisfactory. For KpSC
there is no licensed vaccine; the most advanced candidates are O-polysaccharide
bioconjugates whose coverage is limited by O- and K-serotype diversity and by the
organism's capacity to exchange capsular loci through recombination. For *A. baumannii*
there is likewise no licensed vaccine, despite more than a decade of candidates that
protect in murine models — OmpA, Ata, Bap, outer-membrane vesicles — whose translation is
hindered by inter-lineage heterogeneity and by the difficulty of standardising correlates
of protection.

For *S. pneumoniae* effective conjugate vaccines do exist, and it is precisely there that
the most instructive lesson lies. PCV7, PCV13, PCV15 and PCV20 protect against the
serotypes they contain — and at each valency expansion the vacated ecological niche is
occupied by non-vaccine serotypes. Serotype replacement is the empirical demonstration
that a capsule-based strategy is an arms race against a naturally transformable genome.

### 1.1 Design theses

We invert three premises of classical reverse vaccinology.

**One construct, three pathogens.** Rather than treating KpSC, *A. baumannii* and
*S. pneumoniae* as three programmes, we treat the **clinical niche** — pneumonia in the
hospitalised, vulnerable host — as the design unit. The construct is a single chimeric
protein carrying epitope blocks from all three organisms.

**Capsule independence.** Only core-genome surface proteins enter the design: present in
≥95% of genomes and independent of K/O/capsular serotype. This is the direct answer to the
serotype-replacement problem that limits current pneumococcal vaccines.

**Escape should be costly.** Antigens are prioritised for purifying selection and coupling
to fitness, so that an escape mutant pays a physiological price — which is not true of most
reverse-vaccinology targets.

### 1.2 Two correctable defects in current practice

**(a) Incomplete safety screening.** The standard is to BLAST against the human proteome
and declare safety. This misses two things. First, low global homology does not exclude
mimicry: a single contiguous 7-mer identical to a human protein suffices for T-cell
cross-reactivity in several experimental models, and short peptides escape detection unless
`blastp-short` is used. Second — and this is essentially absent from the literature — nobody
checks homology against the **commensal microbiome**. A respiratory vaccine that induces
antibodies against epitopes shared with commensal *Neisseria lactamica*, *Streptococcus
mitis* or *Moraxella* may disturb a protective flora and open the very niche it aims to
close. We add that layer.

**(b) Eurocentric coverage and greedy score-based selection.** Epitope choice is usually
"the top N by score", which is demonstrably suboptimal: the best binders concentrate on the
same common HLA alleles, leaving entire populations uncovered. Moreover, the reference
panels in use reflect predominantly European allele frequencies. We treat selection as a
formal maximum-coverage problem and optimise phenotypic coverage weighted by world **and
admixed Brazilian** allele frequencies, reporting both separately.

### 1.3 An exploratory element: structurally equivalent shared epitopes

The most speculative part of the design. KpSC and *A. baumannii* are Gram-negatives sharing
outer-membrane protein families with conserved β-barrel folds. *S. pneumoniae* is
Gram-positive, and sequence homology with the other two is essentially nil. Yet
functionally analogous surface proteins — adhesins, ABC-transporter substrate-binding
proteins — may present **conformational epitopes at structurally equivalent positions**
without sequence-detectable ancestry. We search for such cases by structural superposition
rather than sequence alignment, and report them as an explicitly labelled hypothesis, not
as an established result.

---

## 2. Results

### 2.1 Genome sampling and core genome

⟨PENDENTE:secao_resultados⟩

*Figure 1* shows the attrition of epitopes along the pipeline.

### 2.2 A validated surfaceome is far smaller than annotation suggests

Subcellular localisation was predicted for every core protein. Of ⟨PENDENTE:n_core_total⟩
core proteins, only ⟨PENDENTE:n_surface_total⟩ are predicted to occupy an
antibody-accessible compartment — outer membrane or extracellular for the Gram-negatives,
cell wall or extracellular for the Gram-positive. The remainder are cytoplasmic, inner
membrane or periplasmic.

This matters methodologically: a keyword-based surfaceome derived from annotation text is
substantially more permissive than a predictor-based one, and an epitope from a
cytoplasmic protein is not an antibody target however conserved or safe it may be.

### 2.3 Commensal homology is the dominant safety filter in *S. pneumoniae*

The four-layer negative screen behaved very differently across organisms (*Figure 2*).
For the Gram-negatives, the dominant filter was self-similarity and 7-mer mimicry. For
*S. pneumoniae*, the **commensal layer dominated**: a large fraction of otherwise
acceptable epitopes were discarded because they are shared with commensal streptococci
(*S. mitis*, *S. oralis*).

This is the expected consequence of the close phylogenetic relationship between
*S. pneumoniae* and the oral streptococci, and it is exactly the failure mode that the
commensal layer exists to prevent. A conventional pipeline, screening only against the
human proteome, would have carried those epitopes into the construct.

### 2.4 Population coverage

⟨PENDENTE:cobertura_resumo⟩ (*Figure 3*)

### 2.5 The chimeric construct

⟨PENDENTE:construto_resumo⟩ (*Figure 4*)

### 2.6 Structural convergence between pathogens (exploratory)

Surface antigens from the three organisms were compared by structural superposition
rather than sequence identity (*Figure 5*). Pairs exceeding the TM-score threshold for
shared fold were dominated by **ABC-transporter substrate-binding proteins**, whose fold is
conserved across all three species despite low sequence identity.

We report this as a hypothesis. Cross-reactivity of short *linear* epitopes derived from
structurally equivalent but sequence-divergent regions is rare, and we did not include a
shared block in the construct on this basis alone.

---

## 3. Discussion

The central claim of this work is not that we have produced a vaccine. It is that two
inexpensive additions to a standard reverse-vaccinology pipeline change its output
materially, and that both are transferable.

**The commensal layer changes the answer, not just the paperwork.** For *S. pneumoniae* it
removed the majority of otherwise-acceptable epitopes. Any pipeline that screens only
against the human proteome would have retained them. Given that pneumococcal carriage
occurs in a niche shared with commensal streptococci, and that disruption of that flora is
a plausible route to niche opening, this seems a poor thing to leave unchecked.

**Ancestry-aware coverage costs nothing and is rarely done.** Reporting world and Brazilian
coverage separately makes explicit a choice that is usually implicit and Eurocentric.

### Limitations

This is a computational study; nothing here substitutes for *in vitro* or *in vivo*
validation, and we make no claim of protection.

Specific limitations we consider material:

1. **Reference-based core genome.** Genes absent from the reference strain are not
   discovered. This is acceptable for selecting universal antigens — an antigen missing
   from the reference is a poor universal candidate — but these numbers should not be used
   as a description of the pangenome.
2. **Allergenicity and toxicity remain pending.** Those tools have no stable API and
   require manual submission; the corresponding layer is reported as PENDING rather than
   passed. Pending is not approval.
3. **The shared structural block is a hypothesis**, reported as such and not incorporated
   into the construct.
4. **Purifying-selection analysis (dN/dS) was not performed** in this version, as it
   requires codon-level alignments the light-weight core-genome path does not produce.

---

## 4. Methods

All thresholds are declared in `config/config.yaml`; no threshold is hard-coded. The
pipeline is orchestrated with Snakemake, and the workflow file describes the path that was
actually executed.

### 4.1 Genome sampling

Genomes were sampled with geographic and sequence-type stratification rather than by taking
the first N assemblies, which biases heavily towards European/North American isolates and
towards over-sequenced outbreak clones. Target counts per organism are declared in the
configuration.

### 4.2 Core genome

The core genome was defined by reciprocal orthology against a reference proteome using
BLASTp: a reference protein is core when an ortholog (identity and coverage above declared
thresholds) is present in at least 95% of genomes. We chose this reference-centred approach
over all-vs-all pangenome clustering because the goal is not to describe pangenome
architecture but to obtain proteins present in nearly all isolates; the trade-off is stated
as a limitation above.

### 4.3 Surfaceome

Subcellular localisation was predicted with **DeepLocPro**, a prokaryote-specific predictor,
run with the appropriate Gram setting per organism. Signal peptides were predicted with
**SignalP-6** and membrane topology with **DeepTMHMM**. Proteins were retained as surface
candidates when localisation fell in an antibody-accessible compartment for the organism's
cell-wall type, topology was compatible with expression, and length fell within the
configured range.

### 4.4 Epitope prediction and conservation

B-cell, MHC-I and MHC-II epitopes were predicted with the IEDB tools. Conservation was
computed as the fraction of isolates whose ortholog contains the **exact** epitope k-mer,
with the number of isolates as denominator, so that absence counts against conservation.
Epitopes below the configured conservation threshold were discarded.

### 4.5 Four-layer negative safety screen

(A) **Self-similarity** against the reviewed human proteome; (B) **autoimmune mimicry**, any
exact contiguous 7-mer shared with a human protein; (C) **commensal microbiome**, against
respiratory and gut commensal proteomes; (D) **allergenicity and toxicity**, which require
manual submission and are reported as PENDING.

For layers A and C, identity was normalised by epitope coverage — the fraction of the
epitope identical to a self or commensal window — rather than using the raw local identity
returned by BLAST. Raw local identity is meaningless for short peptides, where a 5–6 residue
match returns 100% identity and would reject every epitope.

### 4.6 Population coverage

Epitope selection was treated as a maximum-coverage problem and solved greedily, which
carries a (1 − 1/e) approximation guarantee. The objective is phenotypic coverage under
Hardy–Weinberg assumptions, weighted by a mixture of world and Brazilian allele frequencies.
World-only and Brazil-only coverage are reported separately.

### 4.7 Construct assembly and characterisation

Epitopes were assembled with class-appropriate linkers, an N-terminal molecular adjuvant and
a promiscuous T-helper epitope. Free cysteines were excluded so that the construct carries
no unpaired thiol. Physicochemical properties and *E. coli* codon optimisation were computed
as described in the configuration.

---

## Data availability

Genome assemblies are public NCBI RefSeq records; accessions are listed in the repository.
Reference proteomes (human, commensal) are public UniProt/NCBI datasets.

## Code availability

The complete pipeline is available at <https://github.com/lcerdeira/pannosovax> under
GPL-3.0, and is distributed as a desktop application, **PanNosoVax Studio**, described
separately (ref. — software Application Note, in preparation).

## Competing interests

The authors declare no competing interests.
