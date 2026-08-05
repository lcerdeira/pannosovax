# PanNosoVax Studio: reproducible multi-epitope vaccine design without the command line

**Target journal:** Bioinformatics (OUP) — Applications Note
**Length budget:** ~1300–2000 words, one figure, ≤20 references

> **Status.** This manuscript is *not submittable yet*. Applications Notes require the
> software to be publicly available and working. Outstanding items are tracked at the end
> of this file; the text below is drafted so that only those items, plus the figure, remain.

---

## Abstract

**Summary.** Reverse vaccinology pipelines are written for people who use a terminal. The
clinicians, epidemiologists and public-health researchers who most need them generally do
not. PanNosoVax Studio is a desktop application that runs a complete multi-epitope vaccine
design pipeline — core genome, surfaceome, epitope prediction, multi-layer safety screening,
HLA population coverage and construct assembly — behind a graphical interface, while
delegating execution to Snakemake so that reproducibility, provenance and resumability are
preserved rather than traded away. Long-running stages checkpoint continuously, so the
application can be closed and reopened without losing work.

**Availability and implementation.** Source code, documentation and a demonstration dataset
are available at <https://github.com/lcerdeira/pannosovax> under GPL-3.0. Packaged builds for
macOS and Windows are produced by continuous integration. Implemented in Python 3.11
(FastAPI, Snakemake, Biopython).

**Contact.** ⟨author contact⟩

---

## 1 Introduction

Reverse vaccinology — deriving vaccine candidates from genome sequences rather than from
cultured organisms — is now a standard computational discipline, and the individual tools
are mature and mostly free. What is not solved is *who can run them*.

Existing options fall into two groups, and both exclude part of the intended audience. Web
servers are approachable but impose submission limits, cannot express a multi-stage
dependent workflow, and give the user no provenance record. Command-line pipelines are
reproducible and scalable but assume comfort with a shell, package managers and cluster
schedulers. The people best positioned to ask good questions about a nosocomial pathogen —
infectious-disease clinicians, hospital epidemiologists, public-health laboratories — are
frequently in neither group.

The obvious answer, a graphical wrapper, is usually the wrong one: wrapping tool invocations
in a GUI typically discards exactly the properties that make a pipeline trustworthy. Our
design point is different. The interface never executes analysis itself; it writes
configuration and drives a **Snakemake** workflow, which remains the single source of truth
for the dependency graph, resumption and provenance. The graphical layer buys accessibility
and costs nothing in reproducibility.

## 2 Implementation

### 2.1 Architecture

The application has three layers:

1. **Interface** — a local web UI (project setup, stage selection, live progress, results
   browser with figures), rendered in a native window.
2. **Control** — a local FastAPI service that translates user choices into Snakemake
   invocations, streams logs, and exposes results.
3. **Execution** — the Snakemake workflow and the underlying analysis scripts.

Nothing is duplicated between layers: the stage list shown to the user is derived from the
workflow itself by a dry run, so the interface cannot drift out of step with the pipeline.
The application runs entirely on the user's machine; no data is uploaded except to the public
prediction services the pipeline already depends on.

### 2.2 Stages exposed

The default target covers acquisition, reference-based core genome, surfaceome prediction,
epitope prediction with conservation filtering, four-layer safety screening, HLA coverage
optimisation, construct assembly, physicochemical characterisation and codon optimisation.

Stages requiring an external service without a stable programmatic interface — 3D structure
prediction, molecular dynamics, immune simulation, allergenicity — are deliberately **not**
in the default path. Instead, the pipeline prepares a submission package and reports the
stage as pending. This is a design decision: a pipeline that appears to complete while
silently skipping a safety layer is worse than one that says it did not run.

### 2.3 Long-running work

Epitope prediction against a public API takes hours at realistic scale. The application
treats this as a product requirement rather than an inconvenience: every prediction is cached
per protein, so an interrupted run resumes where it stopped; the interface states the
expected duration before the user commits; and the workflow's own dependency tracking means
re-running recomputes only what is missing.

### 2.4 Distribution

Builds for macOS and Windows are produced by continuous integration on every release tag, so
the packaged artefacts are always built from the tested tree. A demonstration dataset of nine
genomes ships with the repository and completes in minutes, allowing installation to be
verified without committing to a full-scale run.

## 3 Use case

The application was developed alongside, and validated by, a full-scale design study against
three WHO-priority causes of nosocomial pneumonia — *Klebsiella pneumoniae* species complex,
*Acinetobacter baumannii* and *Streptococcus pneumoniae* — described separately (ref. Paper A).
That study exercises every stage in the default path at production scale and produced a
single chimeric immunogen. Two of its methodological features are implemented as
first-class pipeline stages rather than post-hoc scripts: negative screening against the
commensal microbiome, and epitope selection by maximum population coverage weighted by
non-European HLA allele frequencies.

**Figure 1.** (a) The stage view during execution, showing per-stage state and streaming log.
(b) The workflow dependency graph rendered by Snakemake, illustrating that the graphical layer
drives — and does not replace — the reproducible workflow.

## 4 Conclusion

PanNosoVax Studio makes an established computational method usable by researchers who are
excluded from it by tooling rather than by expertise, without weakening the guarantees that
make the method defensible.

Its principal limitation is inherited from the field: several analyses that a complete design
should include depend on services with no stable programmatic access, and remain manual
steps. The application makes those steps explicit and prepares their inputs, but cannot
remove them.

---

## Blocking items before submission

| item | status |
|---|---|
| Public repository, OSI licence | ✅ GPL-3.0 |
| CI passing on macOS + Windows | ✅ |
| Automated tests | ✅ API contract tests |
| Demonstration dataset | ✅ 9 genomes, minutes |
| Packaged binaries downloadable | ⬜ built by CI; need a tagged release |
| Archived DOI (Zenodo) | ⬜ |
| End-to-end demo through all default stages | ⬜ pending full-scale run settling |
| Figure 1 (interface + DAG) | ⬜ |
| Author list, contact, funding | ⬜ |
| Paper A citation (for cross-reference) | ⬜ depends on Paper A status |
